package gff3

import (
	"strings"
	"testing"
)

func BenchmarkRecordUnmarshal(b *testing.B) {
	line := "ctg123\t.\tgene\t1000\t9000\t.\t+\t.\tID=gene00001;Name=EDEN"
	var r Record
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r.Unmarshal(line)
	}
}

func BenchmarkParseAttributes(b *testing.B) {
	s := "ID=cds00001;Parent=mRNA00001,mRNA00002;Name=edenprotein.1"
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		ParseAttributes(s)
	}
}

func BenchmarkReader(b *testing.B) {
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		r := NewReader(strings.NewReader(canonicalGene))
		for {
			_, err := r.Read()
			if err != nil {
				break
			}
		}
	}
}

func BenchmarkMarshal(b *testing.B) {
	rec := &Record{
		SeqID:      "ctg123",
		Source:     ".",
		Type:       "CDS",
		Start:      1201,
		End:        1500,
		Score:      ".",
		Strand:     "+",
		Phase:      0,
		Attributes: Attributes{"ID": {"cds00001"}, "Parent": {"mRNA00001", "mRNA00002"}},
	}
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		rec.Marshal()
	}
}
