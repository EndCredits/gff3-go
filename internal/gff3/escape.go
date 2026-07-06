package gff3

import (
	"net/url"
	"strings"
)

// Unescape decodes GFF3 Percent-Encoding according to RFC 3986.
func Unescape(s string) string {
	if !strings.Contains(s, "%") {
		return s
	}
	decoded, err := url.QueryUnescape(s)
	if err != nil {
		return s
	}
	return decoded
}

// Escape encodes file-level reserved characters using Percent-Encoding.
// Escapes: tab, newline, carriage return, %, and control characters.
func Escape(s string) string {
	return escapeWith(s, fileMustEscape)
}

// EscapeAttr encodes both file-level and column-9 reserved characters.
// In addition to file-level escaping, also escapes ; = & , which have
// reserved meanings in GFF3 column 9.
func EscapeAttr(s string) string {
	return escapeWith(s, attrMustEscape)
}

func escapeWith(s string, must func(byte) bool) string {
	for i := 0; i < len(s); i++ {
		if must(s[i]) {
			return doEscape(s, must)
		}
	}
	return s
}

func doEscape(s string, must func(byte) bool) string {
	result := make([]byte, 0, len(s))
	for i := 0; i < len(s); i++ {
		c := s[i]
		if must(c) {
			result = append(result, '%')
			upper := c >> 4
			lower := c & 0x0F
			result = append(result, hexDigit(upper))
			result = append(result, hexDigit(lower))
		} else {
			result = append(result, c)
		}
	}
	return string(result)
}

func fileMustEscape(c byte) bool {
	switch c {
	case '\t', '\n', '\r', '%':
		return true
	default:
		return c < 0x20 || c == 0x7F
	}
}

func attrMustEscape(c byte) bool {
	if fileMustEscape(c) {
		return true
	}
	switch c {
	case ';', '=', '&', ',':
		return true
	default:
		return false
	}
}

func hexDigit(v byte) byte {
	if v < 10 {
		return '0' + v
	}
	return 'A' + (v - 10)
}
