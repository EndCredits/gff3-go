package gff3idx

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

type GeneChildren struct {
	Transcripts []string
	CDSs        []string
	Exons       []string
}

type SpatialFeat struct {
	Start int
	End   int
	ID    string
	Type  string
}

type Querier interface {
	ByID(id string) (*Feature, bool)
	ChildrenOf(geneID string) (*GeneChildren, bool)
	InRange(chr string, minStart, maxEnd int) []SpatialFeat
}
