package gff3

import (
	"reflect"
	"testing"
)

func TestParseTarget(t *testing.T) {
	tests := []struct {
		input string
		want  Target
	}{
		{
			input: "cdna0123 12 2964",
			want:  Target{ID: "cdna0123", Start: 12, End: 2964, Strand: "+"},
		},
		{
			input: "cdna0123 12 2964 +",
			want:  Target{ID: "cdna0123", Start: 12, End: 2964, Strand: "+"},
		},
		{
			input: "EST_B 1 500 -",
			want:  Target{ID: "EST_B", Start: 1, End: 500, Strand: "-"},
		},
		{
			input: "EST_C 1 500 +",
			want:  Target{ID: "EST_C", Start: 1, End: 500, Strand: "+"},
		},
		{
			input: "EST_D 1 500 -",
			want:  Target{ID: "EST_D", Start: 1, End: 500, Strand: "-"},
		},
	}

	for _, tt := range tests {
		got, err := ParseTarget(tt.input)
		if err != nil {
			t.Errorf("ParseTarget(%q) error: %v", tt.input, err)
			continue
		}
		if !reflect.DeepEqual(got, tt.want) {
			t.Errorf("ParseTarget(%q) = %+v, want %+v", tt.input, got, tt.want)
		}
	}
}

func TestParseTargetMinimal(t *testing.T) {
	got, err := ParseTarget("target1 1 100")
	if err != nil {
		t.Fatalf("unexpected error: %v", err)
	}
	if got.Strand != "+" {
		t.Errorf("default strand = %q, want +", got.Strand)
	}
}

func TestParseTargetInvalid(t *testing.T) {
	_, err := ParseTarget("target1 1")
	if err == nil {
		t.Error("expected error for insufficient fields")
	}
}

func TestParseTargetInvalidStrand(t *testing.T) {
	_, err := ParseTarget("target1 1 100 X")
	if err == nil {
		t.Error("expected error for invalid strand")
	}
}

func TestParseGap(t *testing.T) {
	tests := []struct {
		input string
		want  []GapOp
	}{
		{
			input: "M8 D3 M6 I1 M6",
			want: []GapOp{
				{Code: 'M', Length: 8},
				{Code: 'D', Length: 3},
				{Code: 'M', Length: 6},
				{Code: 'I', Length: 1},
				{Code: 'M', Length: 6},
			},
		},
		{
			input: "M451 D3499 M501 D1499 M2001",
			want: []GapOp{
				{Code: 'M', Length: 451},
				{Code: 'D', Length: 3499},
				{Code: 'M', Length: 501},
				{Code: 'D', Length: 1499},
				{Code: 'M', Length: 2001},
			},
		},
		{
			input: "M3 I1 M2 D1 M4",
			want: []GapOp{
				{Code: 'M', Length: 3},
				{Code: 'I', Length: 1},
				{Code: 'M', Length: 2},
				{Code: 'D', Length: 1},
				{Code: 'M', Length: 4},
			},
		},
		{
			input: "M3 I1 M2 F1 M4",
			want: []GapOp{
				{Code: 'M', Length: 3},
				{Code: 'I', Length: 1},
				{Code: 'M', Length: 2},
				{Code: 'F', Length: 1},
				{Code: 'M', Length: 4},
			},
		},
		{
			input: "M3 I1 M2 R1 M4",
			want: []GapOp{
				{Code: 'M', Length: 3},
				{Code: 'I', Length: 1},
				{Code: 'M', Length: 2},
				{Code: 'R', Length: 1},
				{Code: 'M', Length: 4},
			},
		},
		{
			input: "",
			want:  nil,
		},
	}

	for _, tt := range tests {
		got, err := ParseGap(tt.input)
		if err != nil {
			t.Errorf("ParseGap(%q) error: %v", tt.input, err)
			continue
		}
		if !reflect.DeepEqual(got, tt.want) {
			t.Errorf("ParseGap(%q) = %+v, want %+v", tt.input, got, tt.want)
		}
	}
}

func TestParseGapInvalidCode(t *testing.T) {
	_, err := ParseGap("X8")
	if err == nil {
		t.Error("expected error for invalid gap code")
	}
}

func TestParseGapInvalidLength(t *testing.T) {
	_, err := ParseGap("M0")
	if err == nil {
		t.Error("expected error for zero length")
	}
}
