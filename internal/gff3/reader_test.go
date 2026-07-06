package gff3

import (
	"io"
	"reflect"
	"strings"
	"testing"
)

func TestParseDirective(t *testing.T) {
	tests := []struct {
		line string
		kind DirectiveKind
		args []string
	}{
		{"##gff-version 3.1.26", DirGFFVersion, []string{"3.1.26"}},
		{"##gff-version 3", DirGFFVersion, []string{"3"}},
		{"##sequence-region ctg123 1 1497228", DirSequenceRegion, []string{"ctg123", "1", "1497228"}},
		{"##feature-ontology http://example.com/sofa.obo", DirFeatureOntology, []string{"http://example.com/sofa.obo"}},
		{"##attribute-ontology http://example.com/attr.obo", DirAttributeOntology, []string{"http://example.com/attr.obo"}},
		{"##source-ontology http://example.com/src.obo", DirSourceOntology, []string{"http://example.com/src.obo"}},
		{"##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=6239", DirSpecies, []string{"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=6239"}},
		{"##genome-build WormBase ws110", DirGenomeBuild, []string{"WormBase", "ws110"}},
		{"###", DirTerminator, nil},
		{"##FASTA", DirFASTA, nil},
		{"##custom-directive arg1 arg2", DirUnknown, []string{"custom-directive", "arg1", "arg2"}},
	}

	for _, tt := range tests {
		d, err := ParseDirective(tt.line)
		if err != nil {
			t.Errorf("ParseDirective(%q) error: %v", tt.line, err)
			continue
		}
		if d.Kind != tt.kind {
			t.Errorf("ParseDirective(%q).Kind = %v, want %v", tt.line, d.Kind, tt.kind)
		}
		if !reflect.DeepEqual(d.Args, tt.args) {
			t.Errorf("ParseDirective(%q).Args = %v, want %v", tt.line, d.Args, tt.args)
		}
	}
}

func TestParseDirectiveInvalid(t *testing.T) {
	_, err := ParseDirective("ctg123\t.\tgene\t1000\t9000")
	if err == nil {
		t.Error("expected error for non-directive line")
	}
}

func TestParseSequenceRegion(t *testing.T) {
	d := Directive{Kind: DirSequenceRegion, Args: []string{"ctg123", "1", "1497228"}}
	sr, err := ParseSequenceRegion(d)
	if err != nil {
		t.Fatalf("ParseSequenceRegion error: %v", err)
	}
	if sr.SeqID != "ctg123" {
		t.Errorf("SeqID = %q, want ctg123", sr.SeqID)
	}
	if sr.Start != 1 {
		t.Errorf("Start = %d, want 1", sr.Start)
	}
	if sr.End != 1497228 {
		t.Errorf("End = %d, want 1497228", sr.End)
	}
}

func TestParseSequenceRegionInvalid(t *testing.T) {
	d := Directive{Kind: DirGFFVersion, Args: []string{"3"}}
	_, err := ParseSequenceRegion(d)
	if err == nil {
		t.Error("expected error for non-sequence-region directive")
	}
}

func TestReaderCanonicalGene(t *testing.T) {
	r := NewReader(strings.NewReader(canonicalGene))

	var records []*Record
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read error: %v", err)
		}
		records = append(records, rec)
	}

	if len(records) != 23 {
		t.Errorf("expected 23 records, got %d", len(records))
	}

	dirs := r.Directives()
	if len(dirs) != 2 {
		t.Errorf("expected 2 directives, got %d", len(dirs))
	}
	if dirs[0].Kind != DirGFFVersion {
		t.Errorf("first directive should be gff-version, got %v", dirs[0].Kind)
	}
	if dirs[1].Kind != DirSequenceRegion {
		t.Errorf("second directive should be sequence-region, got %v", dirs[1].Kind)
	}

	sr, err := ParseSequenceRegion(dirs[1])
	if err != nil {
		t.Fatalf("ParseSequenceRegion error: %v", err)
	}
	if sr.SeqID != "ctg123" {
		t.Errorf("SeqID = %q, want ctg123", sr.SeqID)
	}
}

func TestReaderGeneRecord(t *testing.T) {
	r := NewReader(strings.NewReader(canonicalGene))
	rec, err := r.Read()
	if err != nil {
		t.Fatalf("Read error: %v", err)
	}

	if rec.SeqID != "ctg123" {
		t.Errorf("SeqID = %q", rec.SeqID)
	}
	if rec.Type != "gene" {
		t.Errorf("Type = %q", rec.Type)
	}
	if rec.Start != 1000 {
		t.Errorf("Start = %d", rec.Start)
	}
}

func TestReaderEmptyFile(t *testing.T) {
	r := NewReader(strings.NewReader(""))
	_, err := r.Read()
	if err != io.EOF {
		t.Errorf("expected EOF, got %v", err)
	}
}

func TestReaderOnlyDirectives(t *testing.T) {
	input := "##gff-version 3.1.26\n##sequence-region ctg123 1 1000\n###\n"
	r := NewReader(strings.NewReader(input))
	_, err := r.Read()
	if err != io.EOF {
		t.Errorf("expected EOF, got %v", err)
	}

	dirs := r.Directives()
	if len(dirs) != 3 {
		t.Errorf("expected 3 directives, got %d", len(dirs))
	}
}

func TestReaderBlankLines(t *testing.T) {
	input := "##gff-version 3.1.26\n\n\nctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001\n"
	r := NewReader(strings.NewReader(input))
	rec, err := r.Read()
	if err != nil {
		t.Fatalf("Read error: %v", err)
	}
	if rec.Type != "gene" {
		t.Errorf("Type = %q", rec.Type)
	}
}

func TestReaderComments(t *testing.T) {
	input := "##gff-version 3.1.26\n# this is a comment\nctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001\n"
	r := NewReader(strings.NewReader(input))
	rec, err := r.Read()
	if err != nil {
		t.Fatalf("Read error: %v", err)
	}
	if rec.Type != "gene" {
		t.Errorf("Type = %q", rec.Type)
	}
}

func TestReaderFASTASection(t *testing.T) {
	input := "##gff-version 3.1.26\nctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001\n##FASTA\n>ctg123\nACGT\n"
	r := NewReader(strings.NewReader(input))

	rec, err := r.Read()
	if err != nil {
		t.Fatalf("Read error: %v", err)
	}
	if rec.Type != "gene" {
		t.Errorf("Type = %q", rec.Type)
	}

	_, err = r.Read()
	if err != io.EOF {
		t.Errorf("expected EOF after FASTA, got %v", err)
	}
}

func TestReaderInvalidLine(t *testing.T) {
	input := "##gff-version 3.1.26\nctg123\t.\tgene\t1000\t9000\t.\tX\t.\tID=gene00001\n"
	r := NewReader(strings.NewReader(input))
	_, err := r.Read()
	if err == nil {
		t.Error("expected error for invalid strand")
	}
}

func TestReaderMissingGFFVersion(t *testing.T) {
	input := "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001\n"
	r := NewReader(strings.NewReader(input))
	_, err := r.Read()
	if err == nil {
		t.Error("expected error for missing ##gff-version")
	}
}

func TestReaderTerminatorStopsReading(t *testing.T) {
	input := "##gff-version 3\nctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001\n###\nctg123\t.\tmRNA\t1050\t9000\t.\t+\t.\tID=mRNA00001\n"
	r := NewReader(strings.NewReader(input))

	rec, err := r.Read()
	if err != nil {
		t.Fatalf("Read error: %v", err)
	}
	if rec.Type != "gene" {
		t.Errorf("expected gene, got %s", rec.Type)
	}

	_, err = r.Read()
	if err != io.EOF {
		t.Errorf("expected EOF after ###, got %v", err)
	}

	if len(r.Directives()) != 2 {
		t.Errorf("expected 2 directives, got %d", len(r.Directives()))
	}
}

func TestReaderTerminatorIncludesDirective(t *testing.T) {
	input := "##gff-version 3\n###\n"
	r := NewReader(strings.NewReader(input))
	_, err := r.Read()
	if err != io.EOF {
		t.Errorf("expected EOF, got %v", err)
	}
	dirs := r.Directives()
	if len(dirs) != 2 {
		t.Fatalf("expected 2 directives, got %d", len(dirs))
	}
	if dirs[1].Kind != DirTerminator {
		t.Errorf("expected DirTerminator, got %v", dirs[1].Kind)
	}
}
