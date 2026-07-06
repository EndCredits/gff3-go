package gff3

import (
	"io"
	"reflect"
	"strings"
	"testing"
)

// Canonical Gene example from the GFF3 specification.
const canonicalGene = `##gff-version 3.1.26
##sequence-region ctg123 1 1497228
ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN
ctg123	.	TF_binding_site	1000	1012	.	+	.	ID=tfbs00001;Parent=gene00001
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001;Name=EDEN.1
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00002;Parent=gene00001;Name=EDEN.2
ctg123	.	mRNA	1300	9000	.	+	.	ID=mRNA00003;Parent=gene00001;Name=EDEN.3
ctg123	.	exon	1300	1500	.	+	.	ID=exon00001;Parent=mRNA00003
ctg123	.	exon	1050	1500	.	+	.	ID=exon00002;Parent=mRNA00001,mRNA00002
ctg123	.	exon	3000	3902	.	+	.	ID=exon00003;Parent=mRNA00001,mRNA00003
ctg123	.	exon	5000	5500	.	+	.	ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123	.	exon	7000	9000	.	+	.	ID=exon00005;Parent=mRNA00001,mRNA00002,mRNA00003
ctg123	.	CDS	1201	1500	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123	.	CDS	3000	3902	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123	.	CDS	5000	5500	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123	.	CDS	7000	7600	.	+	0	ID=cds00001;Parent=mRNA00001;Name=edenprotein.1
ctg123	.	CDS	1201	1500	.	+	0	ID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123	.	CDS	5000	5500	.	+	0	ID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123	.	CDS	7000	7600	.	+	0	ID=cds00002;Parent=mRNA00002;Name=edenprotein.2
ctg123	.	CDS	3301	3902	.	+	0	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123	.	CDS	5000	5500	.	+	1	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123	.	CDS	7000	7600	.	+	1	ID=cds00003;Parent=mRNA00003;Name=edenprotein.3
ctg123	.	CDS	3391	3902	.	+	0	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
ctg123	.	CDS	5000	5500	.	+	1	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
ctg123	.	CDS	7000	7600	.	+	1	ID=cds00004;Parent=mRNA00003;Name=edenprotein.4
`

func TestCanonicalGeneUnmarshal(t *testing.T) {
	lines := splitLines(canonicalGene)
	featureLines := 0
	for _, line := range lines {
		if line == "" || line[0] == '#' {
			continue
		}
		var r Record
		err := r.Unmarshal(line)
		if err != nil {
			t.Fatalf("Unmarshal(%q) error: %v", line, err)
		}
		featureLines++
	}
	if featureLines != 23 {
		t.Errorf("expected 23 feature lines, got %d", featureLines)
	}
}

func TestCanonicalGeneGeneRecord(t *testing.T) {
	line := "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN"
	var r Record
	err := r.Unmarshal(line)
	if err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	if r.SeqID != "ctg123" {
		t.Errorf("SeqID = %q, want %q", r.SeqID, "ctg123")
	}
	if r.Source != "." {
		t.Errorf("Source = %q, want %q", r.Source, ".")
	}
	if r.Type != "gene" {
		t.Errorf("Type = %q, want %q", r.Type, "gene")
	}
	if r.Start != 1000 {
		t.Errorf("Start = %d, want 1000", r.Start)
	}
	if r.End != 9000 {
		t.Errorf("End = %d, want 9000", r.End)
	}
	if r.Score != "." {
		t.Errorf("Score = %q, want %q", r.Score, ".")
	}
	if r.Strand != "+" {
		t.Errorf("Strand = %q, want %q", r.Strand, "+")
	}
	if r.Phase != PhaseUndefined {
		t.Errorf("Phase = %d, want %d", r.Phase, PhaseUndefined)
	}
	if got := r.Attributes.Get("ID"); got != "gene00001" {
		t.Errorf("ID = %q, want %q", got, "gene00001")
	}
	if got := r.Attributes.Get("Name"); got != "EDEN" {
		t.Errorf("Name = %q, want %q", got, "EDEN")
	}
}

func TestCanonicalGeneCDSRecord(t *testing.T) {
	line := "ctg123\t.\tCDS\t1201\t1500\t.\t+\t0\tID=cds00001;Parent=mRNA00001;Name=edenprotein.1"
	var r Record
	err := r.Unmarshal(line)
	if err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	if r.Type != "CDS" {
		t.Errorf("Type = %q, want CDS", r.Type)
	}
	if r.Start != 1201 {
		t.Errorf("Start = %d, want 1201", r.Start)
	}
	if r.End != 1500 {
		t.Errorf("End = %d, want 1500", r.End)
	}
	if r.Phase != 0 {
		t.Errorf("Phase = %d, want 0", r.Phase)
	}
	if got := r.Attributes.Get("ID"); got != "cds00001" {
		t.Errorf("ID = %q, want %q", got, "cds00001")
	}
	if got := r.Attributes.Get("Parent"); got != "mRNA00001" {
		t.Errorf("Parent = %q, want %q", got, "mRNA00001")
	}
	if got := r.Attributes.Get("Name"); got != "edenprotein.1" {
		t.Errorf("Name = %q, want %q", got, "edenprotein.1")
	}
}

func TestCanonicalGeneCDSPhase1(t *testing.T) {
	line := "ctg123\t.\tCDS\t5000\t5500\t.\t+\t1\tID=cds00003;Parent=mRNA00003;Name=edenprotein.3"
	var r Record
	err := r.Unmarshal(line)
	if err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	if r.Phase != 1 {
		t.Errorf("Phase = %d, want 1", r.Phase)
	}
}

func TestCanonicalGeneExonMultipleParents(t *testing.T) {
	line := "ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003"
	var r Record
	err := r.Unmarshal(line)
	if err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	parents := r.Attributes["Parent"]
	expected := []string{"mRNA00001", "mRNA00002", "mRNA00003"}
	if !reflect.DeepEqual(parents, expected) {
		t.Errorf("Parents = %v, want %v", parents, expected)
	}
}

func TestUnmarshalInvalidColumnCount(t *testing.T) {
	line := "ctg123\t.\tgene\t1000\t9000"
	var r Record
	err := r.Unmarshal(line)
	if err == nil {
		t.Error("expected error for invalid column count")
	}
}

func TestUnmarshalInvalidStrand(t *testing.T) {
	line := "ctg123\t.\tgene\t1000\t9000\t.\tX\t.\tID=gene00001"
	var r Record
	err := r.Unmarshal(line)
	if err == nil {
		t.Error("expected error for invalid strand")
	}
}

func TestUnmarshalInvalidStart(t *testing.T) {
	line := "ctg123\t.\tgene\tabc\t9000\t.\t+\t.\tID=gene00001"
	var r Record
	err := r.Unmarshal(line)
	if err == nil {
		t.Error("expected error for invalid start")
	}
}

func TestUnmarshalUndefinedStart(t *testing.T) {
	line := "ctg123\t.\tgene\t.\t9000\t.\t+\t.\tID=gene00001"
	var r Record
	err := r.Unmarshal(line)
	if err == nil {
		t.Error("expected error for undefined start coordinate")
	}
}

// Simplified canonical gene example (features without children may omit ID).
const simplifiedGene = `##gff-version 3.1.26
##sequence-region ctg123 1 1497228
ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN
ctg123	.	TF_binding_site	1000	1012	.	+	.	Parent=gene00001
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00001;Parent=gene00001
ctg123	.	mRNA	1050	9000	.	+	.	ID=mRNA00002;Parent=gene00001
ctg123	.	mRNA	1300	9000	.	+	.	ID=mRNA00003;Parent=gene00001
ctg123	.	exon	1300	1500	.	+	.	Parent=mRNA00003
ctg123	.	exon	1050	1500	.	+	.	Parent=mRNA00001,mRNA00002
ctg123	.	exon	3000	3902	.	+	.	Parent=mRNA00001,mRNA00003
ctg123	.	exon	5000	5500	.	+	.	Parent=mRNA00001,mRNA00002,mRNA00003
ctg123	.	exon	7000	9000	.	+	.	Parent=mRNA00001,mRNA00002,mRNA00003
ctg123	.	CDS	1201	1500	.	+	0	ID=cds00001;Parent=mRNA00001
ctg123	.	CDS	3000	3902	.	+	0	ID=cds00001;Parent=mRNA00001
ctg123	.	CDS	5000	5500	.	+	0	ID=cds00001;Parent=mRNA00001
ctg123	.	CDS	7000	7600	.	+	0	ID=cds00001;Parent=mRNA00001
ctg123	.	CDS	1201	1500	.	+	0	ID=cds00002;Parent=mRNA00002
ctg123	.	CDS	5000	5500	.	+	0	ID=cds00002;Parent=mRNA00002
ctg123	.	CDS	7000	7600	.	+	0	ID=cds00002;Parent=mRNA00002
ctg123	.	CDS	3301	3902	.	+	0	ID=cds00003;Parent=mRNA00003
ctg123	.	CDS	5000	5500	.	+	1	ID=cds00003;Parent=mRNA00003
ctg123	.	CDS	7000	7600	.	+	1	ID=cds00003;Parent=mRNA00003
ctg123	.	CDS	3391	3902	.	+	0	ID=cds00004;Parent=mRNA00003
ctg123	.	CDS	5000	5500	.	+	1	ID=cds00004;Parent=mRNA00003
ctg123	.	CDS	7000	7600	.	+	1	ID=cds00004;Parent=mRNA00003
`

func TestSimplifiedGeneNoIDFeature(t *testing.T) {
	r := NewReader(strings.NewReader(simplifiedGene))
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read error: %v", err)
		}
		if rec.Type == "TF_binding_site" {
			if rec.Attributes.Get("ID") != "" {
				t.Error("TF_binding_site should have no ID in simplified gene")
			}
			if rec.Attributes.Get("Parent") != "gene00001" {
				t.Errorf("TF_binding_site Parent = %q", rec.Attributes.Get("Parent"))
			}
		}
	}
}

func TestMultiLineFeatureAssembly(t *testing.T) {
	var records []*Record

	r := NewReader(strings.NewReader(canonicalGene))
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

	groups := GroupByID(records)

	cds1Records := groups["cds00001"]
	if len(cds1Records) != 4 {
		t.Errorf("cds00001 should have 4 segments, got %d", len(cds1Records))
	}

	starts := make([]int, len(cds1Records))
	for i, r := range cds1Records {
		starts[i] = r.Start
	}

	expectedStarts := []int{1201, 3000, 5000, 7000}
	if len(starts) != len(expectedStarts) {
		t.Fatalf("expected %d CDS segments", len(expectedStarts))
	}
	for i := range expectedStarts {
		if starts[i] != expectedStarts[i] {
			t.Errorf("cds00001 segment %d start = %d, want %d", i, starts[i], expectedStarts[i])
		}
	}
}

func TestDetectCycleOnCanonicalGene(t *testing.T) {
	var records []*Record
	r := NewReader(strings.NewReader(canonicalGene))
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

	if err := DetectCycle(records); err != nil {
		t.Errorf("canonical gene should have no cycles: %v", err)
	}
}

func TestRecordClone(t *testing.T) {
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

	c := r.Clone()
	c.Attributes["ID"][0] = "modified"
	if r.Attributes.Get("ID") != "gene00001" {
		t.Error("Clone modified original Attributes")
	}
}

func splitLines(s string) []string {
	var lines []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == '\n' {
			lines = append(lines, s[start:i])
			start = i + 1
		}
	}
	if start < len(s) {
		lines = append(lines, s[start:])
	}
	return lines
}

func TestESTAlignmentTarget(t *testing.T) {
	line := "ctg123\t.\tEST_match\t1200\t3200\t2.2e-30\t+\t.\tID=match00002;Target=mjm1123.5 5 506;Gap=M301 D1499 M201"
	var r Record
	if err := r.Unmarshal(line); err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	target, err := ParseTarget(r.Attributes.Get("Target"))
	if err != nil {
		t.Fatalf("ParseTarget error: %v", err)
	}
	if target.ID != "mjm1123.5" {
		t.Errorf("Target ID = %q", target.ID)
	}
	if target.Start != 5 {
		t.Errorf("Target Start = %d", target.Start)
	}
	if target.End != 506 {
		t.Errorf("Target End = %d", target.End)
	}

	gap, err := ParseGap(r.Attributes.Get("Gap"))
	if err != nil {
		t.Fatalf("ParseGap error: %v", err)
	}
	if len(gap) != 3 {
		t.Errorf("expected 3 gap ops, got %d", len(gap))
	}
}

func TestESTAlignmentMinusStrand(t *testing.T) {
	line := "ctg123\t.\tEST_match\t7000\t9000\t7.4e-32\t-\t.\tID=match00003;Target=mjm1123.3 1 502;Gap=M101 D1499 M401"
	var r Record
	if err := r.Unmarshal(line); err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	if r.Strand != "-" {
		t.Errorf("Strand = %q, want -", r.Strand)
	}

	target, _ := ParseTarget(r.Attributes.Get("Target"))
	if target.End != 502 {
		t.Errorf("Target End = %d, want 502", target.End)
	}
}

func TestCDNAMatchAlignment(t *testing.T) {
	line := "ctg123\t.\tcDNA_match\t1050\t9000\t6.2e-45\t+\t.\tID=match00001;Target=cdna0123 12 2964;Gap=M451 D3499 M501 D1499 M2001"
	var r Record
	if err := r.Unmarshal(line); err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	target, err := ParseTarget(r.Attributes.Get("Target"))
	if err != nil {
		t.Fatalf("ParseTarget error: %v", err)
	}
	if target.ID != "cdna0123" {
		t.Errorf("Target ID = %q", target.ID)
	}

	gap, err := ParseGap(r.Attributes.Get("Gap"))
	if err != nil {
		t.Fatalf("ParseGap error: %v", err)
	}
	if len(gap) != 5 {
		t.Errorf("expected 5 gap ops, got %d", len(gap))
	}
}

func TestCircularGenomeRecord(t *testing.T) {
	line := "J02448\tGenBank\tCDS\t6006\t7238\t.\t+\t0\tID=geneII;Name=II;Note=protein II"
	var r Record
	if err := r.Unmarshal(line); err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	if r.SeqID != "J02448" {
		t.Errorf("SeqID = %q", r.SeqID)
	}
	if r.Source != "GenBank" {
		t.Errorf("Source = %q", r.Source)
	}
	if r.End != 7238 {
		t.Errorf("End = %d, want 7238", r.End)
	}
}

func TestCircularGenomeIsCircular(t *testing.T) {
	line := "J02448\tGenBank\tregion\t1\t6407\t.\t+\t.\tID=J02448;Name=J02448;Is_circular=true"
	var r Record
	if err := r.Unmarshal(line); err != nil {
		t.Fatalf("Unmarshal error: %v", err)
	}

	if r.Attributes.Get("Is_circular") != "true" {
		t.Errorf("Is_circular = %q, want true", r.Attributes.Get("Is_circular"))
	}
}

