package gff3idx

import (
	"os"
	"strings"
	"testing"

	"github.com/EndCredits/gff3-go"
)

func parseTestRecords(data string) []*gff3.Record {
	r := gff3.NewReader(strings.NewReader(data))
	var records []*gff3.Record
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		if rec == nil {
			break
		}
		records = append(records, rec)
	}
	return records
}

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

func TestBuildAndQuery(t *testing.T) {
	records := parseTestRecords(canonicalGene)
	if len(records) == 0 {
		t.Fatal("no records parsed")
	}

	tmp, err := os.CreateTemp("", "gff3idx_test_*.idx")
	if err != nil {
		t.Fatal(err)
	}
	tmpPath := tmp.Name()
	tmp.Close()
	defer os.Remove(tmpPath)

	if err := Build(records, tmpPath); err != nil {
		t.Fatalf("Build error: %v", err)
	}

	r, err := Open(tmpPath)
	if err != nil {
		t.Fatalf("Open error: %v", err)
	}
	defer r.Close()

	if r.EntryCount() == 0 {
		t.Error("empty entry count")
	}
	t.Logf("entries: %d, genes: %d, chrs: %d", r.EntryCount(), r.GeneCount(), r.ChrCount())

	feat, ok := r.ByID("gene00001")
	if !ok {
		t.Fatal("gene00001 not found")
	}
	if feat.Type != "gene" {
		t.Errorf("expected gene, got %s", feat.Type)
	}
	if feat.Start != 1000 || feat.End != 9000 {
		t.Errorf("gene coords: %d-%d, want 1000-9000", feat.Start, feat.End)
	}
	if feat.Strand != "+" {
		t.Errorf("strand: %s, want +", feat.Strand)
	}

	feat, ok = r.ByID("cds00003")
	if !ok {
		t.Fatal("cds00003 not found")
	}
	if feat.Type != "CDS" {
		t.Errorf("expected CDS, got %s", feat.Type)
	}

	_, ok = r.ByID("nonexistent")
	if ok {
		t.Error("nonexistent ID should not be found")
	}
}

func TestChildrenOf(t *testing.T) {
	records := parseTestRecords(canonicalGene)
	tmp, err := os.CreateTemp("", "gff3idx_test_*.idx")
	if err != nil {
		t.Fatal(err)
	}
	tmpPath := tmp.Name()
	tmp.Close()
	defer os.Remove(tmpPath)

	if err := Build(records, tmpPath); err != nil {
		t.Fatalf("Build error: %v", err)
	}

	r, err := Open(tmpPath)
	if err != nil {
		t.Fatalf("Open error: %v", err)
	}
	defer r.Close()

	gc, ok := r.ChildrenOf("gene00001")
	if !ok {
		t.Fatal("gene00001 not found for children")
	}

	if len(gc.Transcripts) != 3 {
		t.Errorf("expected 3 transcripts, got %d: %v", len(gc.Transcripts), gc.Transcripts)
	}
	if len(gc.CDSs) != 4 {
		t.Errorf("expected 4 CDSs, got %d", len(gc.CDSs))
	}
	if len(gc.Exons) != 5 {
		t.Errorf("expected 5 exons, got %d", len(gc.Exons))
	}

	expectedTx := []string{"mRNA00001", "mRNA00002", "mRNA00003"}
	for i, tx := range expectedTx {
		if i >= len(gc.Transcripts) || gc.Transcripts[i] != tx {
			t.Errorf("transcript[%d] = %q, want %q", i, gc.Transcripts[i], tx)
		}
	}
}

func TestSpatialQuery(t *testing.T) {
	records := parseTestRecords(canonicalGene)
	tmp, err := os.CreateTemp("", "gff3idx_test_*.idx")
	if err != nil {
		t.Fatal(err)
	}
	tmpPath := tmp.Name()
	tmp.Close()
	defer os.Remove(tmpPath)

	if err := Build(records, tmpPath); err != nil {
		t.Fatalf("Build error: %v", err)
	}

	r, err := Open(tmpPath)
	if err != nil {
		t.Fatalf("Open error: %v", err)
	}
	defer r.Close()

	feats := r.InRange("ctg123", 1000, 1500)
	if len(feats) == 0 {
		t.Fatal("no features in range")
	}

	for _, f := range feats {
		if f.End < 1000 || f.Start > 1500 {
			t.Errorf("feature %s (%d-%d) outside query range 1000-1500", f.ID, f.Start, f.End)
		}
	}

	t.Logf("features in [1000,1500]:")
	for _, f := range feats {
		t.Logf("  %s %s %d-%d", f.ID, f.Type, f.Start, f.End)
	}
}

func TestFullBuild(t *testing.T) {
	path := os.Getenv("GFF3_TEST_FILE")
	if path == "" {
		t.Skip("GFF3_TEST_FILE not set; export GFF3_TEST_FILE=/path/to/annotations.gff3")
	}

	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	r := gff3.NewReader(f)
	var records []*gff3.Record
	for {
		rec, err := r.Read()
		if err != nil { break }
		if rec == nil { break }
		records = append(records, rec)
	}
	t.Logf("parsed %d records", len(records))

	tmp, err := os.CreateTemp("", "gff3idx_full_*.idx")
	if err != nil {
		t.Fatal(err)
	}
	tmpPath := tmp.Name()
	tmp.Close()
	defer os.Remove(tmpPath)

	if err := Build(records, tmpPath); err != nil {
		t.Fatalf("Build error: %v", err)
	}

	fi, _ := os.Stat(tmpPath)
	t.Logf("index size: %.1f MB", float64(fi.Size())/1e6)

	idx, err := Open(tmpPath)
	if err != nil {
		t.Fatalf("Open error: %v", err)
	}
	defer idx.Close()

	t.Logf("entries: %d, genes: %d, chrs: %d", idx.EntryCount(), idx.GeneCount(), idx.ChrCount())

	if idx.EntryCount() != 983853 {
		t.Errorf("expected 983853 entries (records with ID), got %d", idx.EntryCount())
	}

	feat, ok := idx.ByID("arahy.Tifrunner.gnm2.ann2.Ah01g000200")
	if !ok {
		t.Fatal("Ah01g000200 not found")
	}
	if feat.Type != "gene" {
		t.Errorf("expected gene, got %s", feat.Type)
	}
	if feat.Start != 19126 || feat.End != 25719 {
		t.Errorf("coords: %d-%d, want 19126-25719", feat.Start, feat.End)
	}
	t.Logf("lookup: %s %s %d-%d %s phase=%d", feat.SeqID, feat.Type, feat.Start, feat.End, feat.Strand, feat.Phase)

	gc, ok := idx.ChildrenOf("arahy.Tifrunner.gnm2.ann2.Ah01g000200")
	if !ok {
		t.Fatal("children not found")
	}
	t.Logf("children: %d tx, %d cds, %d exons", len(gc.Transcripts), len(gc.CDSs), len(gc.Exons))
	if len(gc.Transcripts) != 1 {
		t.Errorf("expected 1 transcript, got %d", len(gc.Transcripts))
	}

	feats := idx.InRange("arahy.Tifrunner.gnm2.chr01", 50_000, 60_000)
	t.Logf("chr01 50k-60k: %d features", len(feats))
	if len(feats) == 0 {
		t.Error("expected features in this range")
	}
	for _, f := range feats {
		t.Logf("  %s %s %d-%d", f.ID, f.Type, f.Start, f.End)
	}
}

func TestMemQuerierCanonicalGene(t *testing.T) {
	records := parseTestRecords(canonicalGene)
	m := Wrap(records)

	feat, ok := m.ByID("gene00001")
	if !ok {
		t.Fatal("gene00001 not found")
	}
	if feat.Type != "gene" || feat.Start != 1000 || feat.End != 9000 {
		t.Errorf("gene: type=%s %d-%d", feat.Type, feat.Start, feat.End)
	}

	gc, ok := m.ChildrenOf("gene00001")
	if !ok {
		t.Fatal("children not found")
	}
	if len(gc.Transcripts) != 3 {
		t.Errorf("transcripts: %d, want 3", len(gc.Transcripts))
	}
	if len(gc.CDSs) != 4 {
		t.Errorf("CDSs: %d, want 4", len(gc.CDSs))
	}
	if len(gc.Exons) != 5 {
		t.Errorf("exons: %d, want 5", len(gc.Exons))
	}

	feats := m.InRange("ctg123", 1000, 1500)
	if len(feats) == 0 {
		t.Fatal("no features in range")
	}
	for _, f := range feats {
		if f.End < 1000 || f.Start > 1500 {
			t.Errorf("%s (%d-%d) has no overlap with [1000,1500]", f.ID, f.Start, f.End)
		}
	}
}

func TestMemQuerierVsBinaryIndex(t *testing.T) {
	records := parseTestRecords(canonicalGene)

	tmp, _ := os.CreateTemp("", "gff3idx_test_*.idx")
	tmpPath := tmp.Name()
	tmp.Close()
	defer os.Remove(tmpPath)

	if err := Build(records, tmpPath); err != nil {
		t.Fatal(err)
	}

	idx, _ := Open(tmpPath)
	defer idx.Close()

	mem := Wrap(records)

	ids := []string{"gene00001", "mRNA00001", "cds00001", "exon00004", "tfbs00001", "nonexistent"}
	for _, id := range ids {
		mf, mOk := mem.ByID(id)
		bf, bOk := idx.ByID(id)
		if mOk != bOk {
			t.Errorf("ByID(%s) mem=%v bin=%v", id, mOk, bOk)
			continue
		}
		if mOk && (mf.Type != bf.Type || mf.Start != bf.Start || mf.End != bf.End) {
			t.Errorf("ByID(%s): mem={%s %d-%d} bin={%s %d-%d}", id, mf.Type, mf.Start, mf.End, bf.Type, bf.Start, bf.End)
		}
	}

	gcMem, _ := mem.ChildrenOf("gene00001")
	gcBin, _ := idx.ChildrenOf("gene00001")
	if !strSlicesEq(gcMem.Transcripts, gcBin.Transcripts) ||
		!strSlicesEq(gcMem.CDSs, gcBin.CDSs) ||
		!strSlicesEq(gcMem.Exons, gcBin.Exons) {
		t.Error("ChildrenOf mismatch between mem and binary index")
	}

	for _, q := range [][2]int{{1000, 1500}, {3000, 4000}, {5000, 5500}} {
		memFeats := mem.InRange("ctg123", q[0], q[1])
		binFeats := idx.InRange("ctg123", q[0], q[1])
		if len(memFeats) != len(binFeats) {
			t.Errorf("InRange(%d,%d): mem=%d bin=%d", q[0], q[1], len(memFeats), len(binFeats))
		}
	}
}

func strSlicesEq(a, b []string) bool {
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