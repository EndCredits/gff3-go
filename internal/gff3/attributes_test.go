package gff3

import (
	"reflect"
	"testing"
)

func TestParseAttributes(t *testing.T) {
	tests := []struct {
		name  string
		input string
		want  Attributes
	}{
		{
			name:  "empty dot",
			input: ".",
			want:  Attributes{},
		},
		{
			name:  "empty string",
			input: "",
			want:  Attributes{},
		},
		{
			name:  "single tag",
			input: "ID=gene00001",
			want:  Attributes{"ID": {"gene00001"}},
		},
		{
			name:  "two tags",
			input: "ID=gene00001;Name=EDEN",
			want:  Attributes{"ID": {"gene00001"}, "Name": {"EDEN"}},
		},
		{
			name:  "single parent",
			input: "ID=exon00001;Parent=mRNA00003",
			want:  Attributes{"ID": {"exon00001"}, "Parent": {"mRNA00003"}},
		},
		{
			name:  "multiple parents",
			input: "ID=exon00002;Parent=mRNA00001,mRNA00002",
			want:  Attributes{"ID": {"exon00002"}, "Parent": {"mRNA00001", "mRNA00002"}},
		},
		{
			name:  "three parents",
			input: "ID=exon00004;Parent=mRNA00001,mRNA00002,mRNA00003",
			want:  Attributes{"ID": {"exon00004"}, "Parent": {"mRNA00001", "mRNA00002", "mRNA00003"}},
		},
		{
			name:  "with name",
			input: "ID=gene00001;Name=EDEN",
			want:  Attributes{"ID": {"gene00001"}, "Name": {"EDEN"}},
		},
		{
			name:  "CDS with name",
			input: "ID=cds00001;Parent=mRNA00001;Name=edenprotein.1",
			want:  Attributes{"ID": {"cds00001"}, "Parent": {"mRNA00001"}, "Name": {"edenprotein.1"}},
		},
		{
			name:  "TF binding site",
			input: "ID=tfbs00001;Parent=gene00001",
			want:  Attributes{"ID": {"tfbs00001"}, "Parent": {"gene00001"}},
		},
		{
			name:  "no ID single exon",
			input: "Parent=mRNA00003",
			want:  Attributes{"Parent": {"mRNA00003"}},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseAttributes(tt.input)
			if err != nil {
				t.Fatalf("ParseAttributes(%q) error: %v", tt.input, err)
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("ParseAttributes(%q) = %v, want %v", tt.input, got, tt.want)
			}
		})
	}
}

func TestAttributesGet(t *testing.T) {
	attrs := Attributes{
		"ID":     {"gene00001"},
		"Parent": {"mRNA00001", "mRNA00002"},
	}

	if got := attrs.Get("ID"); got != "gene00001" {
		t.Errorf("Get(ID) = %q, want %q", got, "gene00001")
	}
	if got := attrs.Get("Name"); got != "" {
		t.Errorf("Get(Name) = %q, want empty", got)
	}
}

func TestAttributesClone(t *testing.T) {
	orig := Attributes{
		"ID":     {"gene00001"},
		"Parent": {"mRNA00001", "mRNA00002"},
	}
	clone := orig.Clone()
	clone["ID"][0] = "modified"
	if orig.Get("ID") != "gene00001" {
		t.Error("Clone modified original")
	}
}

func TestParseAttributesWithEscapedChars(t *testing.T) {
	attrs, err := ParseAttributes("Note=contains%3B semicolon%2C comma%3D equals")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	got := attrs.Get("Note")
	want := "contains; semicolon, comma= equals"
	if got != want {
		t.Errorf("got %q, want %q", got, want)
	}
}
