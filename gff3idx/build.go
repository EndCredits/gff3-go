package gff3idx

import (
	"encoding/binary"
	"fmt"
	"os"
	"sort"
	"strings"

	"github.com/EndCredits/gff3-go"

	"github.com/zeebo/xxh3"
)

type builder struct {
	f       *os.File
	strings strings.Builder
	poolCur uint64
	pool    map[string]uint32
}

// Build constructs a binary index file from parsed GFF3 records.
//
// The resulting file can be opened with Open() for repeated mmap-based
// queries. Index size is approximately 150 MB per million records.
func Build(records []*gff3.Record, outPath string) error {
	b := &builder{pool: make(map[string]uint32)}

	f, err := os.Create(outPath)
	if err != nil {
		return fmt.Errorf("create: %w", err)
	}
	defer f.Close()
	b.f = f

	reserve := int64(HeaderSize)
	if err := f.Truncate(reserve); err != nil {
		return fmt.Errorf("truncate: %w", err)
	}
	f.Seek(reserve, 0)

	entryRecs, geneRecs, spatial := b.collect(records)
	entrySlots := b.buildHashSlots(entryRecs)
	geneSlots := b.buildGeneHashSlots(geneRecs)

	entryOff := uint64(reserve)
	slotN := b.writeHashTable(entrySlots)
	entryRecOff := entryOff + uint64(slotN)*HashSlotSize
	b.writeEntryRecords(entryRecs)

	geneOff := entryRecOff + uint64(len(entryRecs))*EntryRecordSize
	b.writeHashTable(geneSlots)

	geneOffsets := make([]uint64, len(geneRecs))
	for i, gr := range geneRecs {
		off, _ := f.Seek(0, 1)
		geneOffsets[i] = uint64(off)
		binary.Write(f, byteOrder, &gr.rec)
	}
	geneDataOff := uint64(b.size())

	geneDataEnd := b.writeGeneData(geneRecs, geneOffsets, geneDataOff)

	spatialOff := geneDataEnd
	spatialEnd := b.writeSpatial(spatial, spatialOff)

	poolOff := spatialEnd
	poolSize := b.writeStringPool()

	hdr := Header{
		EntryCount:     uint32(len(entryRecs)),
		SpatialChrs:    uint32(len(spatial)),
		GeneCount:      uint32(len(geneRecs)),
		StringPoolSize: poolSize,
		EntriesOffset:  entryOff,
		SpatialOffset:  spatialOff,
		GenesOffset:    geneOff,
		StringPoolOff:  poolOff,
	}
	copy(hdr.Magic[:], Magic)
	hdr.Version = Version

	f.Seek(0, 0)
	binary.Write(f, byteOrder, &hdr)
	f.Sync()
	return nil
}

type indexedEntry struct {
	id  string
	rec EntryRecord
}

type indexedGene struct {
	id          string
	rec         GeneRecord
	transcripts []string
	cdss        []string
	exons       []string
}

type indexedSpatialChr struct {
	chr      string
	features []SpatialFeatureRec
}

type idRecord struct {
	r        *gff3.Record
	minStart int64
	maxEnd   int64
}

func (b *builder) collect(records []*gff3.Record) ([]indexedEntry, []indexedGene, []indexedSpatialChr) {
	seen := make(map[string]*idRecord)
	parentToChildren := make(map[string][]string)

	for _, r := range records {
		id := r.Attributes.Get("ID")
		if id == "" {
			continue
		}
		start := int64(r.Start)
		end := int64(r.End)
		if prev, ok := seen[id]; ok {
			if start < prev.minStart {
				prev.minStart = start
			}
			if end > prev.maxEnd {
				prev.maxEnd = end
			}
		} else {
			seen[id] = &idRecord{r: r, minStart: start, maxEnd: end}
		}
		for _, pid := range r.Attributes["Parent"] {
			parentToChildren[pid] = append(parentToChildren[pid], id)
		}
	}

	entries := make([]indexedEntry, 0, len(seen))
	for id, ir := range seen {
		r := ir.r
		entries = append(entries, indexedEntry{
			id: id,
			rec: EntryRecord{
				Start:        ir.minStart,
				End:          ir.maxEnd,
				ChrOffset:    b.stringRef(r.SeqID),
				SourceOffset: b.stringRef(r.Source),
				TypeOffset:   b.stringRef(r.Type),
				ScoreOffset:  b.stringRef(r.Score),
				StrandOffset: b.stringRef(r.Strand),
				Phase:        int32(r.Phase),
			},
		})
	}
	sort.Slice(entries, func(i, j int) bool { return entries[i].id < entries[j].id })

	genes := make([]indexedGene, 0)
	for id, ir := range seen {
		if ir.r.Type != "gene" {
			continue
		}
		txSet := make(map[string]bool)
		cdsSet := make(map[string]bool)
		exonSet := make(map[string]bool)

		for _, childID := range parentToChildren[id] {
			childIR, ok := seen[childID]
			if !ok {
				continue
			}
			switch childIR.r.Type {
			case "mRNA":
				txSet[childID] = true
				for _, grandchildID := range parentToChildren[childID] {
					gcIR, ok2 := seen[grandchildID]
					if !ok2 {
						continue
					}
					switch gcIR.r.Type {
					case "CDS":
						cdsSet[grandchildID] = true
					case "exon":
						exonSet[grandchildID] = true
					}
				}
			case "CDS":
				cdsSet[childID] = true
			case "exon":
				exonSet[childID] = true
			}
		}

		transcripts := sortedKeys(txSet)
		cdss := sortedKeys(cdsSet)
		exons := sortedKeys(exonSet)
		genes = append(genes, indexedGene{
			id: id,
			rec: GeneRecord{
				TranscriptCount: uint32(len(transcripts)),
				CDSCount:        uint32(len(cdss)),
				ExonCount:       uint32(len(exons)),
			},
			transcripts: transcripts,
			cdss:        cdss,
			exons:       exons,
		})
	}
	sort.Slice(genes, func(i, j int) bool { return genes[i].id < genes[j].id })

	chrMap := make(map[string][]SpatialFeatureRec)
	for _, r := range records {
		id := r.Attributes.Get("ID")
		if id == "" {
			continue
		}
		chrMap[r.SeqID] = append(chrMap[r.SeqID], SpatialFeatureRec{
			Start:      int64(r.Start),
			End:        int64(r.End),
			IDOffset:   b.stringRef(id),
			TypeOffset: b.stringRef(r.Type),
		})
	}
	chrKeys := make([]string, 0, len(chrMap))
	for k := range chrMap {
		chrKeys = append(chrKeys, k)
	}
	sort.Strings(chrKeys)
	spatial := make([]indexedSpatialChr, len(chrKeys))
	for i, chr := range chrKeys {
		feats := chrMap[chr]
		sort.Slice(feats, func(a, b int) bool { return feats[a].Start < feats[b].Start })
		spatial[i] = indexedSpatialChr{chr: chr, features: feats}
	}

	return entries, genes, spatial
}

func (b *builder) buildHashSlots(entries []indexedEntry) []HashSlot {
	cap := nextPow2(uint32(len(entries)) * 2)
	slots := make([]HashSlot, cap)
	for i, e := range entries {
		h := hashNonZero(xxh3.HashString(e.id))
		pos := h % uint64(cap)
		for slots[pos].Hash != 0 {
			pos = (pos + 1) % uint64(cap)
		}
		slots[pos].Hash = h
		slots[pos].Val = (uint64(i) << 32) | uint64(b.stringRef(e.id))
	}
	return slots
}

func (b *builder) buildGeneHashSlots(genes []indexedGene) []HashSlot {
	cap := nextPow2(uint32(len(genes)) * 2)
	slots := make([]HashSlot, cap)
	for i, g := range genes {
		h := hashNonZero(xxh3.HashString(g.id))
		pos := h % uint64(cap)
		for slots[pos].Hash != 0 {
			pos = (pos + 1) % uint64(cap)
		}
		slots[pos].Hash = h
		slots[pos].Val = (uint64(i) << 32) | uint64(b.stringRef(g.id))
	}
	return slots
}

func (b *builder) writeHashTable(slots []HashSlot) int {
	for _, s := range slots {
		binary.Write(b.f, byteOrder, &s)
	}
	return len(slots)
}

func (b *builder) writeEntryRecords(entries []indexedEntry) {
	for _, e := range entries {
		binary.Write(b.f, byteOrder, &e.rec)
	}
}

func (b *builder) writeGeneData(genes []indexedGene, offsets []uint64, dataStart uint64) uint64 {
	cur := dataStart
	for i, g := range genes {
		g.rec.DataOffset = cur
		total := g.rec.TranscriptCount + g.rec.CDSCount + g.rec.ExonCount
		buf := make([]uint32, 0, total)
		for _, t := range g.transcripts {
			buf = append(buf, b.stringRef(t))
		}
		for _, c := range g.cdss {
			buf = append(buf, b.stringRef(c))
		}
		for _, e := range g.exons {
			buf = append(buf, b.stringRef(e))
		}
		for _, v := range buf {
			binary.Write(b.f, byteOrder, v)
		}
		cur += uint64(len(buf) * 4)

		off, _ := b.f.Seek(0, 1)
		b.f.Seek(int64(offsets[i]), 0)
		binary.Write(b.f, byteOrder, g.rec)
		b.f.Seek(off, 0)
	}
	return cur
}

func (b *builder) writeSpatial(chrs []indexedSpatialChr, start uint64) uint64 {
	binary.Write(b.f, byteOrder, uint32(len(chrs)))
	cur := start + 4 + uint64(len(chrs))*SpatialHeaderSize
	headers := make([]SpatialHeader, len(chrs))
	for i, c := range chrs {
		headers[i] = SpatialHeader{
			ChrOffset:    b.stringRef(c.chr),
			FeatureCount: uint32(len(c.features)),
			DataOffset:   cur,
		}
		cur += uint64(len(c.features)) * SpatialFeatureRecSize
	}
	for _, h := range headers {
		binary.Write(b.f, byteOrder, &h)
	}
	for _, c := range chrs {
		for _, f := range c.features {
			binary.Write(b.f, byteOrder, &f)
		}
	}
	return uint64(b.size())
}

func (b *builder) writeStringPool() uint64 {
	b.f.Write([]byte(b.strings.String()))
	return uint64(b.strings.Len())
}

func (b *builder) size() int {
	b.f.Sync()
	off, _ := b.f.Seek(0, 1)
	return int(off)
}

func (b *builder) stringRef(s string) uint32 {
	if s == "" {
		return 0xFFFFFFFF
	}
	if off, ok := b.pool[s]; ok {
		return off
	}
	off := b.poolCur
	b.pool[s] = uint32(off)
	b.strings.WriteString(s)
	b.strings.WriteByte(0)
	b.poolCur += uint64(len(s)) + 1
	return uint32(off)
}

func sortedKeys(m map[string]bool) []string {
	result := make([]string, 0, len(m))
	for k := range m {
		result = append(result, k)
	}
	sort.Strings(result)
	return result
}
