package gff3

import (
	"io"
	"strings"
	"testing"
)

func TestMarshalRecord(t *testing.T) {
	r := &Record{
		SeqID:      "ctg123",
		Source:     ".",
		Type:       "gene",
		Start:      1000,
		End:        9000,
		Score:      ".",
		Strand:     "+",
		Phase:      PhaseUndefined,
		Attributes: Attributes{"ID": {"gene00001"}, "Name": {"EDEN"}},
	}

	line, err := r.Marshal()
	if err != nil {
		t.Fatalf("Marshal error: %v", err)
	}

	expected := "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN"
	if line != expected {
		t.Errorf("Marshal = %q, want %q", line, expected)
	}
}

func TestMarshalCDSRecord(t *testing.T) {
	r := &Record{
		SeqID:      "ctg123",
		Source:     ".",
		Type:       "CDS",
		Start:      1201,
		End:        1500,
		Score:      ".",
		Strand:     "+",
		Phase:      0,
		Attributes: Attributes{"ID": {"cds00001"}, "Parent": {"mRNA00001"}},
	}

	line, err := r.Marshal()
	if err != nil {
		t.Fatalf("Marshal error: %v", err)
	}

	expected := "ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001"
	if line != expected {
		t.Errorf("Marshal = %q, want %q", line, expected)
	}
}

func TestMarshalRecordMissingFields(t *testing.T) {
	r := &Record{}
	_, err := r.Marshal()
	if err == nil {
		t.Error("expected error for invalid record")
	}
}

func TestWriterCanonicalGene(t *testing.T) {
	var buf strings.Builder
	w := NewWriter(&buf)
	r := NewReader(strings.NewReader(canonicalGene))

	var dirs []Directive
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read error: %v", err)
		}
		if err := w.WriteRecord(rec); err != nil {
			t.Fatalf("WriteRecord error: %v", err)
		}
	}
	_ = dirs
}

func TestRoundTripDeep(t *testing.T) {
	tests := []struct {
		name  string
		input string
	}{
		{
			name: "canonical gene (full)",
			input: canonicalGene,
		},
		{
			name: "simplified gene",
			input: simplifiedGene,
		},
		{
			name: "EST alignment pair",
			input: `##gff-version 3.1.26
##sequence-region ctg123 1 1497228
ctg123	.	EST_match	1200	3200	2.2e-30	+	.	ID=match00002;Target=mjm1123.5 5 506;Gap=M301 D1499 M201
ctg123	.	EST_match	7000	9000	7.4e-32	-	.	ID=match00003;Target=mjm1123.3 1 502;Gap=M101 D1499 M401
`,
		},
		{
			name: "cDNA match with gap",
			input: `##gff-version 3.1.26
##sequence-region ctg123 1 1497228
ctg123	.	cDNA_match	1050	9000	6.2e-45	+	.	ID=match00001;Target=cdna0123 12 2964;Gap=M451 D3499 M501 D1499 M2001
`,
		},
		{
			name: "circular genome",
			input: `##gff-version 3.1.26
J02448	GenBank	region	1	6407	.	+	.	ID=J02448;Name=J02448;Is_circular=true
J02448	GenBank	CDS	6006	7238	.	+	0	ID=geneII;Name=II;Note=protein II
`,
		},
		{
			name: "mixed types and phases",
			input: `##gff-version 3.1.26
##sequence-region chr1 1 10000
chr1	.	gene	1000	5000	.	+	.	ID=g1;Name=gene1
chr1	.	mRNA	1000	5000	.	+	.	ID=t1;Parent=g1
chr1	.	exon	1000	1200	.	+	.	ID=e1;Parent=t1
chr1	.	CDS	1000	1200	.	+	0	ID=cds1;Parent=t1
chr1	.	CDS	2000	2500	.	+	2	ID=cds1;Parent=t1
chr1	.	mRNA	1000	5000	.	-	.	ID=t2;Parent=g1
chr1	.	exon	4000	5000	.	-	.	Parent=t2
chr1	.	CDS	4000	5000	.	-	1	ID=cds2;Parent=t2
`,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			firstPass := NewReader(strings.NewReader(tt.input))
			var records []*Record
			for {
				rec, err := firstPass.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatalf("first pass Read error: %v", err)
				}
				records = append(records, rec)
			}
			directives := firstPass.Directives()

			var buf strings.Builder
			w := NewWriter(&buf)
			for _, d := range directives {
				if err := w.WriteDirective(d); err != nil {
					t.Fatalf("WriteDirective error: %v", err)
				}
			}
			for _, rec := range records {
				if err := w.WriteRecord(rec); err != nil {
					t.Fatalf("WriteRecord error: %v", err)
				}
			}

			secondPass := NewReader(strings.NewReader(buf.String()))
			var records2 []*Record
			for {
				rec, err := secondPass.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatalf("second pass Read error: %v", err)
				}
				records2 = append(records2, rec)
			}
			directives2 := secondPass.Directives()

			if len(records) != len(records2) {
				t.Fatalf("record count mismatch: %d vs %d", len(records), len(records2))
			}
			if len(directives) != len(directives2) {
				t.Fatalf("directive count mismatch: %d vs %d", len(directives), len(directives2))
			}

			for i := range directives {
				if directives[i].Kind != directives2[i].Kind {
					t.Errorf("directive %d Kind: %v vs %v", i, directives[i].Kind, directives2[i].Kind)
				}
				if !stringSlicesEqual(directives[i].Args, directives2[i].Args) {
					t.Errorf("directive %d Args: %v vs %v", i, directives[i].Args, directives2[i].Args)
				}
			}

			for i := range records {
				err := recordsDeepEqual(records[i], records2[i])
				if err != "" {
					t.Errorf("record %d: %s", i, err)
				}
			}
		})
	}
}

func recordsDeepEqual(a, b *Record) string {
	if a.SeqID != b.SeqID {
		return "SeqID mismatch: " + a.SeqID + " vs " + b.SeqID
	}
	if a.Source != b.Source {
		return "Source mismatch: " + a.Source + " vs " + b.Source
	}
	if a.Type != b.Type {
		return "Type mismatch: " + a.Type + " vs " + b.Type
	}
	if a.Start != b.Start {
		return "Start mismatch"
	}
	if a.End != b.End {
		return "End mismatch"
	}
	if a.Score != b.Score {
		return "Score mismatch: " + a.Score + " vs " + b.Score
	}
	if a.Strand != b.Strand {
		return "Strand mismatch: " + a.Strand + " vs " + b.Strand
	}
	if a.Phase != b.Phase {
		return "Phase mismatch"
	}
	if len(a.Attributes) != len(b.Attributes) {
		return "Attribute count mismatch"
	}
	for k, va := range a.Attributes {
		vb, ok := b.Attributes[k]
		if !ok {
			return "missing attribute " + k
		}
		if len(va) != len(vb) {
			return "attribute " + k + " value count mismatch"
		}
		for i := range va {
			if va[i] != vb[i] {
				return "attribute " + k + " value[" + string(rune('0'+i)) + "] mismatch: " + va[i] + " vs " + vb[i]
			}
		}
	}
	return ""
}

func stringSlicesEqual(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestValidate(t *testing.T) {
	tests := []struct {
		name string
		r    *Record
		ok   bool
	}{
		{
			name: "valid gene",
			r:    &Record{SeqID: "chr1", Source: ".", Type: "gene", Start: 1000, End: 2000, Score: ".", Strand: "+", Phase: PhaseUndefined},
			ok:   true,
		},
		{
			name: "valid CDS",
			r:    &Record{SeqID: "chr1", Source: ".", Type: "CDS", Start: 100, End: 300, Score: ".", Strand: "-", Phase: 0},
			ok:   true,
		},
		{
			name: "start > end",
			r:    &Record{SeqID: "chr1", Source: ".", Type: "gene", Start: 2000, End: 1000, Score: ".", Strand: "+", Phase: PhaseUndefined},
			ok:   false,
		},
		{
			name: "invalid strand",
			r:    &Record{SeqID: "chr1", Source: ".", Type: "gene", Start: 1000, End: 2000, Score: ".", Strand: "X", Phase: PhaseUndefined},
			ok:   false,
		},
		{
			name: "negative coordinate",
			r:    &Record{SeqID: "chr1", Source: ".", Type: "gene", Start: -1, End: 2000, Score: ".", Strand: "+", Phase: PhaseUndefined},
			ok:   false,
		},
		{
			name: "empty seqid",
			r:    &Record{SeqID: "", Source: ".", Type: "gene", Start: 1000, End: 2000, Score: ".", Strand: "+", Phase: PhaseUndefined},
			ok:   false,
		},
	}

	for _, tt := range tests {
		err := tt.r.Validate()
		if tt.ok && err != nil {
			t.Errorf("%s: expected valid, got %v", tt.name, err)
		}
		if !tt.ok && err == nil {
			t.Errorf("%s: expected error, got nil", tt.name)
		}
	}
}

func TestGroupByID(t *testing.T) {
	r1 := &Record{Attributes: Attributes{"ID": {"gene00001"}}}
	r2 := &Record{Attributes: Attributes{"ID": {"cds00001"}}}
	r3 := &Record{Attributes: Attributes{"ID": {"cds00001"}}}
	r4 := &Record{Attributes: Attributes{}}

	groups := GroupByID([]*Record{r1, r2, r3, r4})
	if len(groups) != 2 {
		t.Errorf("expected 2 groups, got %d", len(groups))
	}
	if len(groups["cds00001"]) != 2 {
		t.Errorf("expected 2 records for cds00001, got %d", len(groups["cds00001"]))
	}
}

func TestDetectCycle(t *testing.T) {
	records := []*Record{
		{Attributes: Attributes{"ID": {"gene1"}, "Parent": nil}},
		{Attributes: Attributes{"ID": {"mRNA1"}, "Parent": {"gene1"}}},
	}
	if err := DetectCycle(records); err != nil {
		t.Errorf("unexpected cycle error: %v", err)
	}
}

func TestDetectCycleFound(t *testing.T) {
	records := []*Record{
		{Attributes: Attributes{"ID": {"gene1"}, "Parent": {"mRNA1"}}},
		{Attributes: Attributes{"ID": {"mRNA1"}, "Parent": {"gene1"}}},
	}
	if err := DetectCycle(records); err == nil {
		t.Error("expected cycle error")
	}
}

func TestDetectCycleSelfReferencing(t *testing.T) {
	records := []*Record{
		{Attributes: Attributes{"ID": {"gene1"}, "Parent": {"gene1"}}},
	}
	if err := DetectCycle(records); err == nil {
		t.Error("expected cycle error for self-reference")
	}
}
