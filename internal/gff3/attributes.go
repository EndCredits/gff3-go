package gff3

import (
	"fmt"
	"strings"
)

// Attributes holds the parsed tag=value pairs from column 9 of a GFF3 line.
type Attributes map[string][]string

// Get returns the first value for the given tag, or empty string if not present.
func (a Attributes) Get(tag string) string {
	vs := a[tag]
	if len(vs) == 0 {
		return ""
	}
	return vs[0]
}

// Clone returns a deep copy of the Attributes.
func (a Attributes) Clone() Attributes {
	c := make(Attributes, len(a))
	for k, v := range a {
		vv := make([]string, len(v))
		copy(vv, v)
		c[k] = vv
	}
	return c
}

// ParseAttributes parses a column 9 attributes string into a tag→values map.
//
// The format is tag=value pairs separated by semicolons. Multiple values for
// the same tag are separated by commas (only for Parent, Alias, Note, Dbxref,
// and Ontology_term; other tags treat commas as literal value content).
// Percent-encoded characters are decoded after splitting on reserved delimiters.
func ParseAttributes(s string) (Attributes, error) {
	if s == "." || s == "" {
		return Attributes{}, nil
	}

	attrs := make(Attributes, 8)
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] != ';' {
			continue
		}
		if err := parsePair(s[start:i], attrs); err != nil {
			return nil, err
		}
		start = i + 1
	}
	if start < len(s) {
		if err := parsePair(s[start:], attrs); err != nil {
			return nil, err
		}
	}
	return attrs, nil
}

func parsePair(pair string, attrs Attributes) error {
	if pair == "" {
		return nil
	}
	eq := strings.IndexByte(pair, '=')
	if eq < 0 {
		return fmt.Errorf("gff3: attribute %q missing '=' separator", pair)
	}
	tag := Unescape(pair[:eq])
	rawValue := pair[eq+1:]

	if multiValueTag(tag) {
		vals := splitMultiValues(rawValue)
		for _, v := range vals {
			attrs[tag] = append(attrs[tag], Unescape(v))
		}
	} else {
		attrs[tag] = append(attrs[tag], Unescape(rawValue))
	}
	return nil
}

func splitMultiValues(s string) []string {
	n := 1
	for i := 0; i < len(s); i++ {
		if s[i] == ',' {
			n++
		}
	}
	result := make([]string, 0, n)
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == ',' {
			result = append(result, s[start:i])
			start = i + 1
		}
	}
	if start < len(s) {
		result = append(result, s[start:])
	}
	return result
}

func multiValueTag(tag string) bool {
	switch tag {
	case "Parent", "Alias", "Note", "Dbxref", "Ontology_term":
		return true
	default:
		return false
	}
}
