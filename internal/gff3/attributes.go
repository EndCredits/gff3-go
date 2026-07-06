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
	attrs := make(Attributes)
	if s == "." || s == "" {
		return attrs, nil
	}

	pairs := splitAttrPairs(s)
	for _, pair := range pairs {
		if pair == "" {
			continue
		}
		eq := strings.IndexByte(pair, '=')
		if eq < 0 {
			return nil, fmt.Errorf("gff3: attribute %q missing '=' separator", pair)
		}
		tag := Unescape(pair[:eq])
		rawValue := pair[eq+1:]

		if multiValueTag(tag) {
			for _, v := range splitMultiValues(rawValue) {
				attrs[tag] = append(attrs[tag], Unescape(v))
			}
		} else {
			attrs[tag] = append(attrs[tag], Unescape(rawValue))
		}
	}
	return attrs, nil
}

func splitAttrPairs(s string) []string {
	var result []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == ';' {
			result = append(result, s[start:i])
			start = i + 1
		}
	}
	if start < len(s) {
		result = append(result, s[start:])
	}
	return result
}

func splitMultiValues(s string) []string {
	var result []string
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
