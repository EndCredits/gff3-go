package gff3

import "testing"

func TestEscape(t *testing.T) {
	tests := []struct {
		input    string
		expected string
	}{
		{"hello", "hello"},
		{"hello\tworld", "hello%09world"},
		{"line\nbreak", "line%0Abreak"},
		{"100% sure", "100%25 sure"},
		{"control\x00char", "control%00char"},
	}

	for _, tt := range tests {
		got := Escape(tt.input)
		if got != tt.expected {
			t.Errorf("Escape(%q) = %q, want %q", tt.input, got, tt.expected)
		}
	}
}

func TestUnescape(t *testing.T) {
	tests := []struct {
		input    string
		expected string
	}{
		{"hello", "hello"},
		{"hello%09world", "hello\tworld"},
		{"line%0Abreak", "line\nbreak"},
		{"line%0abreak", "line\nbreak"},
		{"100%25sure", "100%sure"},
		{"colon%3Bsemicolon", "colon;semicolon"},
		{"equal%3Dsign", "equal=sign"},
		{"amp%26ersand", "amp&ersand"},
		{"comma%2Cseparated", "comma,separated"},
		{"keep+plus+literal", "keep+plus+literal"},
		{"coord%2Bshift", "coord+shift"},
		{".", "."},
		{"", ""},
		{"incomplete%", "incomplete%"},
		{"invalid%XX", "invalid%XX"},
	}

	for _, tt := range tests {
		got := Unescape(tt.input)
		if got != tt.expected {
			t.Errorf("Unescape(%q) = %q, want %q", tt.input, got, tt.expected)
		}
	}
}

func TestEscapeRoundTrip(t *testing.T) {
	inputs := []string{
		"simple text",
		"text with spaces",
		"gene00001",
		"EDEN.1",
		"Name=edenprotein.1",
		"100% sure",
	}

	for _, input := range inputs {
		encoded := Escape(input)
		decoded := Unescape(encoded)
		if decoded != input {
			t.Errorf("roundtrip failed: %q -> %q -> %q", input, encoded, decoded)
		}
	}
}
