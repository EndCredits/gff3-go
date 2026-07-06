package gff3

import (
	"encoding/json"
	"flag"
	"io"
	"os"
	"os/exec"
	"strings"
	"testing"
)

var gff3File = flag.String("gff3", "", "path to GFF3 file for integration testing")

type crossStats struct {
	TotalRecords int            `json:"total_records"`
	TypeCounts   map[string]int `json:"type_counts"`
	SourceCounts map[string]int `json:"source_counts"`
	StrandCounts map[string]int `json:"strand_counts"`
	UniqueSeqIDs int            `json:"unique_seqids"`
	Errors       int            `json:"errors"`
}

func TestCrossValidateWithPython(t *testing.T) {
	if *gff3File == "" {
		t.Skip("no -gff3 flag provided; use: go test -run CrossValidate -args -gff3 <file>")
	}

	goStats := runGoParser(t, *gff3File)

	t.Logf("Go parser: %d records, %d types, %d sources, %d seqIDs, %d errors",
		goStats.TotalRecords, len(goStats.TypeCounts), len(goStats.SourceCounts),
		goStats.UniqueSeqIDs, goStats.Errors)
	for typ, count := range goStats.TypeCounts {
		t.Logf("  Go  type %s: %d", typ, count)
	}

	pyStats := runPythonLineSplit(t, *gff3File)
	t.Logf("Python line-split: %d records, %d types, %d sources, %d seqIDs, %d errors",
		pyStats.TotalRecords, len(pyStats.TypeCounts), len(pyStats.SourceCounts),
		pyStats.UniqueSeqIDs, pyStats.Errors)
	for typ, count := range pyStats.TypeCounts {
		t.Logf("  Py  type %s: %d", typ, count)
	}

	if goStats.Errors > 0 {
		t.Errorf("Go parser had %d errors", goStats.Errors)
	}
	if pyStats.Errors > 0 {
		t.Errorf("Python parser had %d errors", pyStats.Errors)
	}
	if goStats.TotalRecords != pyStats.TotalRecords {
		t.Errorf("record count mismatch: Go=%d Python=%d", goStats.TotalRecords, pyStats.TotalRecords)
	}
	if goStats.UniqueSeqIDs != pyStats.UniqueSeqIDs {
		t.Errorf("seqID count mismatch: Go=%d Python=%d", goStats.UniqueSeqIDs, pyStats.UniqueSeqIDs)
	}

	for typ, goCount := range goStats.TypeCounts {
		pyCount := pyStats.TypeCounts[typ]
		if goCount != pyCount {
			t.Errorf("type %q count mismatch: Go=%d Python=%d", typ, goCount, pyCount)
		}
	}
	for typ, pyCount := range pyStats.TypeCounts {
		if _, ok := goStats.TypeCounts[typ]; !ok {
			t.Errorf("type %q found in Python (%d) but not in Go", typ, pyCount)
		}
	}

	for src, goCount := range goStats.SourceCounts {
		pyCount := pyStats.SourceCounts[src]
		if goCount != pyCount {
			t.Errorf("source %q count mismatch: Go=%d Python=%d", src, goCount, pyCount)
		}
	}

	for s, goCount := range goStats.StrandCounts {
		pyCount := pyStats.StrandCounts[s]
		if goCount != pyCount {
			t.Errorf("strand %q count mismatch: Go=%d Python=%d", s, goCount, pyCount)
		}
	}

	t.Logf("cross-validation PASSED")
}

func TestCrossValidateBcbio(t *testing.T) {
	if *gff3File == "" {
		t.Skip("no -gff3 flag provided")
	}

	goStats := runGoParser(t, *gff3File)
	bcbioStats := runPythonBcbio(t, *gff3File)

	t.Logf("Go parser:    %d records, %d types", goStats.TotalRecords, len(goStats.TypeCounts))
	t.Logf("BCBio parser: %d records, %d types", bcbioStats.TotalRecords, len(bcbioStats.TypeCounts))

	for typ, bcCount := range bcbioStats.TypeCounts {
		goCount := goStats.TypeCounts[typ]
		t.Logf("  %s: Go=%d BCBio=%d", typ, goCount, bcCount)
	}

	for typ, goCount := range goStats.TypeCounts {
		bcCount := bcbioStats.TypeCounts[typ]
		if goCount != bcCount {
			t.Errorf("type %q count mismatch: Go=%d BCBio=%d", typ, goCount, bcCount)
		}
	}
	for typ, bcCount := range bcbioStats.TypeCounts {
		if _, ok := goStats.TypeCounts[typ]; !ok {
			t.Errorf("type %q found in BCBio (%d) but not in Go", typ, bcCount)
		}
	}

	t.Logf("BCBio cross-validation PASSED")
}

func TestFullParseCycleAndGroup(t *testing.T) {
	if *gff3File == "" {
		t.Skip("no -gff3 flag provided")
	}

	f, err := os.Open(*gff3File)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	r := NewReader(f)
	var records []*Record
	for {
		rec, err := r.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Errorf("parse error: %v", err)
			continue
		}
		if rec == nil {
			break
		}
		records = append(records, rec)
	}

	groups := GroupByID(records)
	multiLineFeatures := 0
	for _, recs := range groups {
		if len(recs) > 1 {
			multiLineFeatures++
		}
	}

	if err := DetectCycle(records); err != nil {
		t.Errorf("cycle detected: %v", err)
	}

	t.Logf("total records: %d", len(records))
	t.Logf("features with ID: %d", len(groups))
	t.Logf("multi-line features: %d", multiLineFeatures)
	t.Logf("directives: %d", len(r.Directives()))
}

func TestRoundTripDeepFile(t *testing.T) {
	if *gff3File == "" {
		t.Skip("no -gff3 flag provided")
	}

	f, err := os.Open(*gff3File)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	r := NewReader(f)
	var records []*Record
	for {
		rec, err := r.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Fatalf("parse error: %v", err)
		}
		if rec == nil {
			break
		}
		records = append(records, rec)
		if len(records) >= 5000 {
			break
		}
	}
	directives := r.Directives()

	var buf strings.Builder
	for _, d := range directives {
		buf.WriteString(directiveToLine(d))
		buf.WriteByte('\n')
	}
	failures := 0
	for i, rec := range records {
		line, err := rec.Marshal()
		if err != nil {
			t.Fatalf("record %d Marshal error: %v", i, err)
		}
		buf.WriteString(line)
		buf.WriteByte('\n')
	}
	output := buf.String()

	r2 := NewReader(strings.NewReader(output))
	var records2 []*Record
	for {
		rec, err := r2.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Fatalf("second pass Read error: %v", err)
		}
		if rec == nil {
			break
		}
		records2 = append(records2, rec)
	}

	if len(records) != len(records2) {
		t.Fatalf("count mismatch: %d vs %d", len(records), len(records2))
	}

	for i := range records {
		if errStr := recordsDeepEqual(records[i], records2[i]); errStr != "" {
			failures++
			if failures <= 10 {
				t.Errorf("record %d: %s", i, errStr)
			}
		}
	}
	if failures > 0 {
		t.Errorf("round-trip failures: %d / %d records", failures, len(records))
	} else {
		t.Logf("round-trip PASS: %d records verified", len(records))
	}
}

func runGoParser(t *testing.T, path string) crossStats {
	f, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()

	r := NewReader(f)
	s := crossStats{
		TypeCounts:   make(map[string]int),
		SourceCounts: make(map[string]int),
		StrandCounts: make(map[string]int),
	}
	seqIDs := make(map[string]bool)

	for {
		rec, err := r.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			s.Errors++
			continue
		}
		if rec == nil {
			break
		}
		s.TotalRecords++
		s.TypeCounts[rec.Type]++
		s.SourceCounts[rec.Source]++
		s.StrandCounts[rec.Strand]++
		seqIDs[rec.SeqID] = true
	}
	s.UniqueSeqIDs = len(seqIDs)
	return s
}

func runPythonLineSplit(t *testing.T, path string) crossStats {
	cmd := exec.Command("python3", "scripts/validate_gff3.py", path)
	cmd.Dir = "."
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("python line-split failed: %v\n%s", err, string(out))
	}
	var s crossStats
	if err := json.Unmarshal(out, &s); err != nil {
		t.Fatalf("json decode: %v\n%s", err, string(out))
	}
	return s
}

func runPythonBcbio(t *testing.T, path string) crossStats {
	cmd := exec.Command("python3", "scripts/validate_gff3.py", path, "--bcbio")
	cmd.Dir = "."
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("python bcbio failed: %v\n%s", err, string(out))
	}
	var s crossStats
	if err := json.Unmarshal(out, &s); err != nil {
		t.Fatalf("json decode: %v\n%s", err, string(out))
	}
	return s
}

func TestMain(m *testing.M) {
	flag.Parse()
	os.Exit(m.Run())
}

func directiveToLine(d Directive) string {
	switch d.Kind {
	case DirGFFVersion:
		if len(d.Args) > 0 {
			return "##gff-version " + d.Args[0]
		}
		return "##gff-version"
	case DirSequenceRegion:
		if len(d.Args) >= 3 {
			return "##sequence-region " + d.Args[0] + " " + d.Args[1] + " " + d.Args[2]
		}
		return "##sequence-region"
	default:
		return "##" + d.Kind.String()
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
				return "attribute " + k + " value mismatch: " + va[i] + " vs " + vb[i]
			}
		}
	}
	return ""
}
