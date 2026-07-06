package gff3

import (
	"bufio"
	"fmt"
	"io"
	"strings"
)

func ParseDirective(line string) (Directive, error) {
	if !strings.HasPrefix(line, "##") {
		return Directive{}, fmt.Errorf("gff3: not a directive line")
	}
	if line == "###" {
		return Directive{Kind: DirTerminator}, nil
	}
	if line == "##FASTA" {
		return Directive{Kind: DirFASTA}, nil
	}

	parts := strings.SplitN(line[2:], " ", 2)
	keyword := parts[0]
	var args []string
	if len(parts) > 1 {
		args = strings.Fields(parts[1])
	}

	switch keyword {
	case "gff-version":
		return Directive{Kind: DirGFFVersion, Args: args}, nil
	case "sequence-region":
		return Directive{Kind: DirSequenceRegion, Args: args}, nil
	case "feature-ontology":
		return Directive{Kind: DirFeatureOntology, Args: args}, nil
	case "attribute-ontology":
		return Directive{Kind: DirAttributeOntology, Args: args}, nil
	case "source-ontology":
		return Directive{Kind: DirSourceOntology, Args: args}, nil
	case "species":
		return Directive{Kind: DirSpecies, Args: args}, nil
	case "genome-build":
		return Directive{Kind: DirGenomeBuild, Args: args}, nil
	default:
		return Directive{Kind: DirUnknown, Args: append([]string{keyword}, args...)}, nil
	}
}

type SequenceRegion struct {
	SeqID string
	Start int
	End   int
}

func ParseSequenceRegion(d Directive) (SequenceRegion, error) {
	if d.Kind != DirSequenceRegion {
		return SequenceRegion{}, fmt.Errorf("gff3: not a sequence-region directive")
	}
	if len(d.Args) < 3 {
		return SequenceRegion{}, fmt.Errorf("gff3: sequence-region requires seqid, start, end")
	}
	start, err := parseInt(d.Args[1])
	if err != nil {
		return SequenceRegion{}, fmt.Errorf("gff3: invalid sequence-region start: %w", err)
	}
	end, err := parseInt(d.Args[2])
	if err != nil {
		return SequenceRegion{}, fmt.Errorf("gff3: invalid sequence-region end: %w", err)
	}
	return SequenceRegion{SeqID: d.Args[0], Start: start, End: end}, nil
}

type Reader struct {
	sc          *bufio.Scanner
	lineNum     int
	directives  []Directive
	inFASTA     bool
	terminated  bool
	sawVersion  bool
	firstLine   bool
}

func NewReader(r io.Reader) *Reader {
	return &Reader{sc: bufio.NewScanner(r), firstLine: true}
}

func (r *Reader) Read() (*Record, error) {
	if r.inFASTA || r.terminated {
		return nil, io.EOF
	}

	for r.sc.Scan() {
		r.lineNum++
		line := r.sc.Text()

		if line == "" {
			continue
		}

		if line[0] == '#' {
			if err := r.handleMeta(line); err != nil {
				return nil, r.wrapErr(err)
			}
			if r.inFASTA || r.terminated {
				return nil, io.EOF
			}
			continue
		}

		if r.firstLine && !r.sawVersion {
			return nil, r.wrapErr(fmt.Errorf("gff3: missing ##gff-version directive on first line"))
		}
		r.firstLine = false

		if line[0] == '>' {
			r.inFASTA = true
			return nil, io.EOF
		}

		var rec Record
		if err := rec.Unmarshal(line); err != nil {
			return nil, r.wrapErr(err)
		}
		return &rec, nil
	}

	if err := r.sc.Err(); err != nil {
		return nil, err
	}
	return nil, io.EOF
}

func (r *Reader) Directives() []Directive {
	return r.directives
}

func (r *Reader) handleMeta(line string) error {
	if len(line) > 1 && line[1] == '#' {
		d, err := ParseDirective(line)
		if err != nil {
			return err
		}
		r.directives = append(r.directives, d)
		if d.Kind == DirFASTA {
			r.inFASTA = true
		}
		if d.Kind == DirTerminator {
			r.terminated = true
		}
		if d.Kind == DirGFFVersion {
			r.sawVersion = true
		}
	}
	return nil
}

func (r *Reader) wrapErr(err error) error {
	return fmt.Errorf("gff3: line %d: %w", r.lineNum, err)
}
