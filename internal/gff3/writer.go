package gff3

import (
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// Writer writes GFF3 formatted records and directives to an io.Writer.
type Writer struct {
	w io.Writer
}

func NewWriter(w io.Writer) *Writer {
	return &Writer{w: w}
}

func (w *Writer) WriteDirective(d Directive) error {
	var line string
	switch d.Kind {
	case DirGFFVersion:
		line = "##gff-version " + strings.Join(d.Args, " ")
	case DirSequenceRegion:
		line = "##sequence-region " + strings.Join(d.Args, " ")
	case DirFeatureOntology:
		line = "##feature-ontology " + strings.Join(d.Args, " ")
	case DirAttributeOntology:
		line = "##attribute-ontology " + strings.Join(d.Args, " ")
	case DirSourceOntology:
		line = "##source-ontology " + strings.Join(d.Args, " ")
	case DirSpecies:
		line = "##species " + strings.Join(d.Args, " ")
	case DirGenomeBuild:
		line = "##genome-build " + strings.Join(d.Args, " ")
	case DirTerminator:
		line = "###"
	case DirFASTA:
		line = "##FASTA"
	case DirUnknown:
		line = "##" + strings.Join(d.Args, " ")
	default:
		return fmt.Errorf("gff3: unknown directive kind %v", d.Kind)
	}
	_, err := fmt.Fprintln(w.w, line)
	return err
}

func (w *Writer) WriteRecord(r *Record) error {
	line, err := r.Marshal()
	if err != nil {
		return err
	}
	_, err = fmt.Fprintln(w.w, line)
	return err
}

// Marshal serializes a Record back into a GFF3 feature line.
func (r *Record) Marshal() (string, error) {
	if err := r.Validate(); err != nil {
		return "", err
	}

	var b strings.Builder

	b.WriteString(Escape(r.SeqID))
	b.WriteByte('\t')
	b.WriteString(Escape(r.Source))
	b.WriteByte('\t')
	b.WriteString(Escape(r.Type))
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(r.Start))
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(r.End))
	b.WriteByte('\t')
	b.WriteString(Escape(r.Score))
	b.WriteByte('\t')
	b.WriteString(r.Strand)
	b.WriteByte('\t')
	if r.Phase == PhaseUndefined {
		b.WriteByte('.')
	} else {
		b.WriteByte('0' + byte(r.Phase))
	}
	b.WriteByte('\t')
	b.WriteString(marshalAttributes(r.Attributes))

	return b.String(), nil
}

func marshalAttributes(attrs Attributes) string {
	if len(attrs) == 0 {
		return "."
	}
	tags := make([]string, 0, len(attrs))
	for tag := range attrs {
		tags = append(tags, tag)
	}
	sort.Strings(tags)

	var parts []string
	for _, tag := range tags {
		values := attrs[tag]
		if multiValueTag(tag) {
			parts = append(parts, Escape(tag)+"="+escapeMultiValues(values))
		} else {
			for _, v := range values {
				parts = append(parts, Escape(tag)+"="+EscapeAttr(v))
			}
		}
	}
	return strings.Join(parts, ";")
}

func escapeMultiValues(values []string) string {
	escaped := make([]string, len(values))
	for i, v := range values {
		escaped[i] = EscapeAttr(v)
	}
	return strings.Join(escaped, ",")
}
