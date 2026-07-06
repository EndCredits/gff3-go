package gff3

import (
	"fmt"
	"strings"
)

// PhaseUndefined is the sentinel value for features that have no phase (not CDS).
const PhaseUndefined = -1

// Strand values.
const (
	StrandPlus  = "+"
	StrandMinus = "-"
	StrandNone  = "."
	StrandUnknown = "?"
)

// LineType distinguishes between kinds of GFF3 lines.
type LineType int

const (
	LineDirective LineType = iota
	LineComment
	LineFeature
	LineBlank
	LineFASTA
)

// DirectiveKind identifies the type of a ## directive line.
type DirectiveKind int

const (
	DirGFFVersion DirectiveKind = iota
	DirSequenceRegion
	DirFeatureOntology
	DirAttributeOntology
	DirSourceOntology
	DirSpecies
	DirGenomeBuild
	DirTerminator
	DirFASTA
	DirUnknown
)

// Directive represents a ##-prefixed meta-data line.
type Directive struct {
	Kind DirectiveKind
	Args []string
}

func (k DirectiveKind) String() string {
	switch k {
	case DirGFFVersion:
		return "gff-version"
	case DirSequenceRegion:
		return "sequence-region"
	case DirFeatureOntology:
		return "feature-ontology"
	case DirAttributeOntology:
		return "attribute-ontology"
	case DirSourceOntology:
		return "source-ontology"
	case DirSpecies:
		return "species"
	case DirGenomeBuild:
		return "genome-build"
	case DirTerminator:
		return "###"
	case DirFASTA:
		return "FASTA"
	default:
		return "unknown"
	}
}

// Record represents a single GFF3 feature line (9 columns).
type Record struct {
	SeqID      string
	Source     string
	Type       string
	Start      int
	End        int
	Score      string
	Strand     string
	Phase      int
	Attributes Attributes
}

func (r *Record) Validate() error {
	if r.SeqID == "" {
		return fmt.Errorf("gff3: seqid must not be empty")
	}
	if strings.Contains(r.SeqID, " ") {
		return fmt.Errorf("gff3: seqid contains unescaped whitespace")
	}
	if strings.HasPrefix(r.SeqID, ">") {
		return fmt.Errorf("gff3: seqid must not begin with '>'")
	}
	if r.Start > r.End {
		return fmt.Errorf("gff3: start (%d) must be <= end (%d)", r.Start, r.End)
	}
	if r.Start < 1 || r.End < 1 {
		return fmt.Errorf("gff3: coordinates must be positive integers")
	}
	if !validStrand(r.Strand) {
		return fmt.Errorf("gff3: invalid strand %q", r.Strand)
	}
	if r.Phase != PhaseUndefined && (r.Phase < 0 || r.Phase > 2) {
		return fmt.Errorf("gff3: invalid phase %d", r.Phase)
	}
	if r.Type == "CDS" && r.Phase == PhaseUndefined {
		return fmt.Errorf("gff3: CDS features require a defined phase (0, 1, or 2)")
	}
	if r.Source == "" {
		return fmt.Errorf("gff3: source must not be empty")
	}
	if r.Type == "" {
		return fmt.Errorf("gff3: type must not be empty")
	}
	return nil
}

// Clone returns a deep copy of the Record.
func (r *Record) Clone() *Record {
	c := *r
	c.Attributes = r.Attributes.Clone()
	return &c
}

// Unmarshal parses a single GFF3 feature line into a Record.
func (r *Record) Unmarshal(line string) error {
	columns := splitColumns(line)
	if len(columns) != 9 {
		return fmt.Errorf("gff3: expected 9 tab-separated columns, got %d", len(columns))
	}

	r.SeqID = Unescape(columns[0])
	r.Source = Unescape(columns[1])
	r.Type = Unescape(columns[2])

	start, err := parseInt(columns[3])
	if err != nil {
		return fmt.Errorf("gff3: invalid start coordinate: %w", err)
	}
	r.Start = start

	end, err := parseInt(columns[4])
	if err != nil {
		return fmt.Errorf("gff3: invalid end coordinate: %w", err)
	}
	r.End = end

	r.Score = Unescape(columns[5])

	strand := Unescape(columns[6])
	if !validStrand(strand) {
		return fmt.Errorf("gff3: invalid strand %q", strand)
	}
	r.Strand = strand

	phase, err := parsePhase(columns[7])
	if err != nil {
		return fmt.Errorf("gff3: invalid phase: %w", err)
	}
	r.Phase = phase

	r.Attributes, err = ParseAttributes(columns[8])
	if err != nil {
		return fmt.Errorf("gff3: invalid attributes: %w", err)
	}

	return nil
}

func splitColumns(line string) []string {
	result := make([]string, 0, 9)
	start := 0
	for i := 0; i < len(line); i++ {
		if line[i] == '\t' {
			result = append(result, line[start:i])
			start = i + 1
		}
	}
	result = append(result, line[start:])
	return result
}

func parseInt(s string) (int, error) {
	if s == "." {
		return 0, fmt.Errorf("gff3: undefined coordinate")
	}
	n := 0
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c < '0' || c > '9' {
			return 0, fmt.Errorf("gff3: invalid integer %q", s)
		}
		n = n*10 + int(c-'0')
	}
	return n, nil
}

func parsePhase(s string) (int, error) {
	if s == "." {
		return PhaseUndefined, nil
	}
	switch s {
	case "0":
		return 0, nil
	case "1":
		return 1, nil
	case "2":
		return 2, nil
	default:
		return 0, fmt.Errorf("gff3: invalid phase %q", s)
	}
}

func validStrand(s string) bool {
	return s == StrandPlus || s == StrandMinus || s == StrandNone || s == StrandUnknown
}
