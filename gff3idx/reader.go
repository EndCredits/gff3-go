package gff3idx

import (
	"fmt"
	"os"
	"unsafe"

	"github.com/zeebo/xxh3"
	"golang.org/x/sys/unix"
)

type Reader struct {
	data []byte
	f    *os.File
	hdr  Header
}

func Open(path string) (*Reader, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("gff3idx open: %w", err)
	}

	fi, err := f.Stat()
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("gff3idx stat: %w", err)
	}
	sz := fi.Size()
	if sz < int64(HeaderSize) {
		f.Close()
		return nil, fmt.Errorf("gff3idx too small: %d bytes", sz)
	}

	data, err := unix.Mmap(int(f.Fd()), 0, int(sz), unix.PROT_READ, unix.MAP_SHARED)
	if err != nil {
		f.Close()
		return nil, fmt.Errorf("gff3idx mmap: %w", err)
	}

	r := &Reader{data: data, f: f}
	hdr := (*Header)(unsafe.Pointer(&data[0]))
	if string(hdr.Magic[:]) != Magic {
		r.Close()
		return nil, fmt.Errorf("gff3idx bad magic: %q", string(hdr.Magic[:]))
	}
	if hdr.Version != Version {
		r.Close()
		return nil, fmt.Errorf("gff3idx unsupported version: %d", hdr.Version)
	}
	r.hdr = *hdr

	if err := r.validateBounds(uint64(sz)); err != nil {
		r.Close()
		return nil, err
	}

	return r, nil
}

func (r *Reader) validateBounds(fileSize uint64) error {
	h := r.hdr
	if h.EntryCount > 100_000_000 || h.GeneCount > 100_000_000 {
		return fmt.Errorf("gff3idx implausible entry/gene count")
	}

	check := func(off uint64, size uint64, name string) error {
		if off+size > fileSize {
			return fmt.Errorf("gff3idx %s out of bounds: offset=%d size=%d file=%d", name, off, size, fileSize)
		}
		return nil
	}

	eslots := uint64(nextPow2(h.EntryCount * 2))
	if err := check(h.EntriesOffset, eslots*HashSlotSize+uint64(h.EntryCount)*EntryRecordSize, "entries"); err != nil {
		return err
	}

	gslots := uint64(nextPow2(h.GeneCount * 2))
	if err := check(h.GenesOffset, gslots*HashSlotSize+uint64(h.GeneCount)*GeneRecordSize, "genes"); err != nil {
		return err
	}

	geneDataStart := h.GenesOffset + gslots*HashSlotSize + uint64(h.GeneCount)*GeneRecordSize
	if geneDataStart < h.SpatialOffset {
		if err := check(geneDataStart, h.SpatialOffset-geneDataStart, "gene data"); err != nil {
			return err
		}
	}

	if err := check(h.StringPoolOff, h.StringPoolSize, "string pool"); err != nil {
		return err
	}

	if err := check(h.SpatialOffset, 4, "spatial header"); err != nil {
		return err
	}

	return nil
}

func (r *Reader) Close() error {
	if r.data != nil {
		unix.Munmap(r.data)
		r.data = nil
	}
	if r.f != nil {
		r.f.Close()
		r.f = nil
	}
	return nil
}

func (r *Reader) stringAt(off uint32) string {
	poolStart := r.hdr.StringPoolOff
	poolEnd := poolStart + r.hdr.StringPoolSize
	if poolEnd > uint64(len(r.data)) {
		return ""
	}
	pool := r.data[poolStart:poolEnd]
	if int(off) >= len(pool) {
		return ""
	}
	end := off
	for end < uint32(len(pool)) && pool[end] != 0 {
		end++
	}
	return string(pool[off:end])
}

func (r *Reader) u32At(off uint64) uint32 {
	if int(off)+4 > len(r.data) {
		return 0
	}
	return byteOrder.Uint32(r.data[off : off+4])
}

func (r *Reader) lookupHash(hashOff uint64, count uint32, id string) (int, bool) {
	h := hashNonZero(xxh3.HashString(id))
	slots := unsafe.Slice((*HashSlot)(unsafe.Pointer(&r.data[hashOff])), count)
	pos := h % uint64(count)
	probed := uint32(0)
	for {
		slot := slots[pos]
		if slot.Hash == 0 {
			return 0, false
		}
		if slot.Hash == h {
			storedOff := uint32(slot.Val & 0xFFFFFFFF)
			if r.stringAt(storedOff) == id {
				return int(slot.Val >> 32), true
			}
		}
		pos = (pos + 1) % uint64(count)
		probed++
		if probed >= count {
			return 0, false
		}
	}
}

func (r *Reader) entrySlotCount() uint32 {
	if r.hdr.EntryCount > 0x7FFFFFFF {
		return 1
	}
	return nextPow2(r.hdr.EntryCount * 2)
}
func (r *Reader) geneSlotCount() uint32 {
	if r.hdr.GeneCount > 0x7FFFFFFF {
		return 1
	}
	return nextPow2(r.hdr.GeneCount * 2)
}

func (r *Reader) ByID(id string) (*Feature, bool) {
	idx, ok := r.lookupHash(r.hdr.EntriesOffset, r.entrySlotCount(), id)
	if !ok {
		return nil, false
	}
	if uint32(idx) >= r.hdr.EntryCount {
		return nil, false
	}

	recOff := r.hdr.EntriesOffset + uint64(r.entrySlotCount())*HashSlotSize + uint64(idx)*EntryRecordSize
	if int(recOff)+EntryRecordSize > len(r.data) {
		return nil, false
	}
	rec := (*EntryRecord)(unsafe.Pointer(&r.data[recOff]))

	return &Feature{
		SeqID:  r.stringAt(rec.ChrOffset),
		Source: r.stringAt(rec.SourceOffset),
		Type:   r.stringAt(rec.TypeOffset),
		Start:  int(rec.Start),
		End:    int(rec.End),
		Score:  r.stringAt(rec.ScoreOffset),
		Strand: r.stringAt(rec.StrandOffset),
		Phase:  int(rec.Phase),
	}, true
}

func (r *Reader) ChildrenOf(geneID string) (*GeneChildren, bool) {
	idx, ok := r.lookupHash(r.hdr.GenesOffset, r.geneSlotCount(), geneID)
	if !ok {
		return nil, false
	}
	if uint32(idx) >= r.hdr.GeneCount {
		return nil, false
	}

	recOff := r.hdr.GenesOffset + uint64(r.geneSlotCount())*HashSlotSize + uint64(idx)*GeneRecordSize
	if int(recOff)+GeneRecordSize > len(r.data) {
		return nil, false
	}
	rec := (*GeneRecord)(unsafe.Pointer(&r.data[recOff]))

	gc := &GeneChildren{
		Transcripts: make([]string, rec.TranscriptCount),
		CDSs:        make([]string, rec.CDSCount),
		Exons:       make([]string, rec.ExonCount),
	}

	off := rec.DataOffset
	n := rec.TranscriptCount
	for i := uint32(0); i < n; i++ {
		gc.Transcripts[i] = r.stringAt(r.u32At(off))
		off += 4
	}
	n = rec.CDSCount
	for i := uint32(0); i < n; i++ {
		gc.CDSs[i] = r.stringAt(r.u32At(off))
		off += 4
	}
	n = rec.ExonCount
	for i := uint32(0); i < n; i++ {
		gc.Exons[i] = r.stringAt(r.u32At(off))
		off += 4
	}
	return gc, true
}

func (r *Reader) InRange(chr string, minStart, maxEnd int) []SpatialFeat {
	off := r.hdr.SpatialOffset
	n := r.u32At(off)
	off += 4
	for i := uint32(0); i < n; i++ {
		if int(off)+SpatialHeaderSize > len(r.data) {
			return nil
		}
		h := (*SpatialHeader)(unsafe.Pointer(&r.data[off]))
		if r.stringAt(h.ChrOffset) == chr {
			featBytes := uint64(h.FeatureCount) * SpatialFeatureRecSize
			if h.DataOffset+featBytes > uint64(len(r.data)) {
				return nil
			}
			recs := unsafe.Slice((*SpatialFeatureRec)(unsafe.Pointer(&r.data[h.DataOffset])), h.FeatureCount)

			lo := sortSearchSpatial(recs, minStart)
			result := make([]SpatialFeat, 0)
			for j := lo; j < int(h.FeatureCount); j++ {
				rec := recs[j]
				if int(rec.Start) > maxEnd {
					break
				}
				result = append(result, SpatialFeat{
					Start: int(rec.Start),
					End:   int(rec.End),
					ID:    r.stringAt(rec.IDOffset),
					Type:  r.stringAt(rec.TypeOffset),
				})
			}
			return result
		}
		off += SpatialHeaderSize
	}
	return nil
}

func sortSearchSpatial(recs []SpatialFeatureRec, target int) int {
	lo, hi := 0, len(recs)-1
	for lo <= hi {
		mid := (lo + hi) / 2
		if int(recs[mid].Start) < target {
			lo = mid + 1
		} else {
			hi = mid - 1
		}
	}
	return lo
}

func (r *Reader) EntryCount() uint32 { return r.hdr.EntryCount }
func (r *Reader) GeneCount() uint32  { return r.hdr.GeneCount }
func (r *Reader) ChrCount() uint32   { return r.hdr.SpatialChrs }
