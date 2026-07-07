// Package gff3idx provides two backends for querying GFF3 features:
// an in-memory index (MemQuerier) and a persistent mmap-based binary index.
// Both implement the Querier interface.
//
// In-memory usage:
//
//	q := gff3idx.Wrap(records)
//	feat, _ := q.ByID("gene00001")
//
// Binary index usage:
//
//	gff3idx.Build(records, "genes.gff3idx")
//	idx, _ := gff3idx.Open("genes.gff3idx")
//	feat, _ := idx.ByID("gene00001")
package gff3idx

// Querier is the common interface for both in-memory and binary index backends.
type Querier interface {
	// ByID looks up a feature by its ID attribute.
	// Returns nil, false if not found.
	ByID(id string) (*Feature, bool)

	// ChildrenOf returns the gene's child features grouped by type:
	// transcripts (mRNA), CDSs, and exons.
	// Traverses the two-level hierarchy gene → mRNA → (CDS, exon).
	ChildrenOf(geneID string) (*GeneChildren, bool)

	// InRange returns all features that overlap the given genomic interval.
	// Returns nil if chromosome not found.
	InRange(chr string, minStart, maxEnd int) []SpatialFeat
}

// Feature represents a single GFF3 feature returned from a query.
type Feature struct {
	SeqID  string
	Source string
	Type   string
	Start  int
	End    int
	Score  string
	Strand string
	Phase  int
}

// GeneChildren holds the IDs of a gene's child features, grouped by type.
type GeneChildren struct {
	Transcripts []string // mRNA IDs
	CDSs        []string // CDS feature IDs
	Exons       []string // exon feature IDs
}

// SpatialFeat is a lightweight feature record returned by InRange queries.
type SpatialFeat struct {
	Start int
	End   int
	ID    string
	Type  string
}
