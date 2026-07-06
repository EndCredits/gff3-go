package gff3

import (
	"io"
	"strings"
	"testing"
)

const testFASTA = `##gff-version 3
ctg123	.	gene	1000	9000	.	+	.	ID=gene00001;Name=EDEN
##FASTA
>ctg123
cttctgggcgtacccgattctcggagaacttgccgcaccattccgccttg
tgttcattgctgcctgcatgttcattgtctacctcggctacgtgtggcta
>cdna0123
ttcaagtgctcagtcaatgtgattcacagtatgtcaccaaatattttggc
agctttctcaagggatcaaaattatggatcattatggaatacctcggtgg
`

func TestReadAllFASTA(t *testing.T) {
	records, err := ReadAllFASTA(strings.NewReader(testFASTA))
	if err != nil {
		t.Fatalf("ReadAllFASTA error: %v", err)
	}

	if len(records) != 2 {
		t.Fatalf("expected 2 FASTA records, got %d", len(records))
	}

	if records[0].ID != "ctg123" {
		t.Errorf("record 0 ID = %q, want ctg123", records[0].ID)
	}
	if records[1].ID != "cdna0123" {
		t.Errorf("record 1 ID = %q, want cdna0123", records[1].ID)
	}
	if len(records[0].Sequence) != 100 {
		t.Errorf("record 0 len = %d, want 100", len(records[0].Sequence))
	}
	if len(records[1].Sequence) != 100 {
		t.Errorf("record 1 len = %d, want 100", len(records[1].Sequence))
	}

	if records[0].Sequence[0] != 'C' {
		t.Errorf("first base = %c, want C", records[0].Sequence[0])
	}
	if records[0].SeqString()[:4] != "CTTC" {
		t.Errorf("first 4 bases = %q", records[0].SeqString()[:4])
	}
}

func TestParseFastaHeader(t *testing.T) {
	tests := []struct {
		input string
		id    string
		desc  string
	}{
		{"ctg123", "ctg123", ""},
		{"ctg123 description here", "ctg123", "description here"},
		{"ctg123  extra  spaces", "ctg123", "extra  spaces"},
		{" ctg123 ", "ctg123", ""},
	}

	for _, tt := range tests {
		id, desc := parseFastaHeader(tt.input)
		if id != tt.id {
			t.Errorf("parseFastaHeader(%q) id = %q, want %q", tt.input, id, tt.id)
		}
		if desc != tt.desc {
			t.Errorf("parseFastaHeader(%q) desc = %q, want %q", tt.input, desc, tt.desc)
		}
	}
}

func TestCleanSequence(t *testing.T) {
	tests := []struct {
		input string
		want  string
	}{
		{"ACGT", "ACGT"},
		{"acgt", "ACGT"},
		{"AcGt", "ACGT"},
		{"ACG T", "ACGT"},
		{"ACG\tT", "ACGT"},
		{"ACG\r\nT", "ACGT"},
		{"ACG123T", "ACG123T"},
	}

	for _, tt := range tests {
		got := string(cleanSequence(tt.input))
		if got != tt.want {
			t.Errorf("cleanSequence(%q) = %q, want %q", tt.input, got, tt.want)
		}
	}
}

func TestReaderWithFASTA(t *testing.T) {
	input := "##gff-version 3\nctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=g1\n##FASTA\n>ctg123\nACGT\n"
	r := NewReader(strings.NewReader(input))

	rec, err := r.Read()
	if err != nil {
		t.Fatalf("Read error: %v", err)
	}
	if rec.Type != "gene" {
		t.Errorf("Type = %q, want gene", rec.Type)
	}

	_, err = r.Read()
	if err != io.EOF {
		t.Errorf("expected EOF after FASTA, got %v", err)
	}

	fasta, err := r.ReadFASTA()
	if err != nil {
		t.Fatalf("ReadFASTA error: %v", err)
	}
	if fasta.ID != "ctg123" {
		t.Errorf("ID = %q, want ctg123", fasta.ID)
	}
	if string(fasta.Sequence) != "ACGT" {
		t.Errorf("Sequence = %q, want ACGT", string(fasta.Sequence))
	}
}

func TestReaderFASTAImplied(t *testing.T) {
	input := "##gff-version 3\nctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=g1\n>ctg123\nACGT\n"
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
		t.Errorf("expected EOF after >, got %v", err)
	}
}

func TestFastaRecordSeqString(t *testing.T) {
	rec := &FastaRecord{ID: "test", Sequence: []byte("ACGT")}
	if rec.SeqString() != "ACGT" {
		t.Errorf("SeqString = %q", rec.SeqString())
	}
}

func TestReadAllFASTAEmpty(t *testing.T) {
	records, err := ReadAllFASTA(strings.NewReader(""))
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if len(records) != 0 {
		t.Errorf("expected 0 records, got %d", len(records))
	}
}
