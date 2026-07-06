package gff3idx

import (
	"encoding/binary"
	"math/bits"
)

const (
	Magic   = "GFFI"
	Version = 1
)

var byteOrder = binary.LittleEndian

func nextPow2(n uint32) uint32 {
	if n == 0 {
		return 1
	}
	return 1 << (32 - bits.LeadingZeros32(n-1))
}

type Header struct {
	Magic          [4]byte
	Version        uint32
	_              uint32
	EntryCount     uint32
	SpatialChrs    uint32
	GeneCount      uint32
	StringPoolSize uint64
	EntriesOffset  uint64
	SpatialOffset  uint64
	GenesOffset    uint64
	StringPoolOff  uint64
}

const HeaderSize = 64

type EntryRecord struct {
	Start        int64
	End          int64
	ChrOffset    uint32
	SourceOffset uint32
	TypeOffset   uint32
	ScoreOffset  uint32
	StrandOffset uint32
	Phase        int32
}

const EntryRecordSize = 40

type GeneRecord struct {
	TranscriptCount uint32
	CDSCount        uint32
	ExonCount       uint32
	_               uint32
	DataOffset      uint64
}

const GeneRecordSize = 24

type SpatialHeader struct {
	ChrOffset    uint32
	FeatureCount uint32
	DataOffset   uint64
}

const SpatialHeaderSize = 16

type SpatialFeatureRec struct {
	Start      int64
	End        int64
	IDOffset   uint32
	TypeOffset uint32
}

const SpatialFeatureRecSize = 24

type HashSlot struct {
	Hash uint64
	Val  uint64
}

const HashSlotSize = 16

func ByteOrder() binary.ByteOrder { return byteOrder }
func NextPow2(n uint32) uint32   { return nextPow2(n) }

func hashNonZero(h uint64) uint64 {
	if h == 0 {
		return 1
	}
	return h
}
