package gff3

import (
	"fmt"
	"strconv"
	"strings"
)

type Target struct {
	ID     string
	Start  int
	End    int
	Strand string
}

func ParseTarget(s string) (Target, error) {
	var t Target
	parts := strings.Fields(s)
	if len(parts) < 3 {
		return t, fmt.Errorf("gff3: target requires at least target_id, start, end")
	}

	t.ID = Unescape(parts[0])

	start, err := strconv.Atoi(parts[1])
	if err != nil {
		return t, fmt.Errorf("gff3: invalid target start: %w", err)
	}
	t.Start = start

	end, err := strconv.Atoi(parts[2])
	if err != nil {
		return t, fmt.Errorf("gff3: invalid target end: %w", err)
	}
	t.End = end

	if len(parts) >= 4 {
		strand := parts[3]
		if strand != "+" && strand != "-" {
			return t, fmt.Errorf("gff3: invalid target strand %q", strand)
		}
		t.Strand = strand
	} else {
		t.Strand = "+"
	}

	return t, nil
}

type GapOp struct {
	Code   byte
	Length int
}

func ParseGap(s string) ([]GapOp, error) {
	if s == "" {
		return nil, nil
	}
	parts := strings.Fields(s)
	ops := make([]GapOp, 0, len(parts))
	for _, part := range parts {
		if len(part) < 2 {
			return nil, fmt.Errorf("gff3: invalid gap segment %q", part)
		}
		code := part[0]
		if !validGapCode(code) {
			return nil, fmt.Errorf("gff3: invalid gap operation %c", code)
		}
		length, err := strconv.Atoi(part[1:])
		if err != nil {
			return nil, fmt.Errorf("gff3: invalid gap length in %q: %w", part, err)
		}
		if length < 1 {
			return nil, fmt.Errorf("gff3: gap length must be positive, got %d", length)
		}
		ops = append(ops, GapOp{Code: code, Length: length})
	}
	return ops, nil
}

func validGapCode(c byte) bool {
	return c == 'M' || c == 'I' || c == 'D' || c == 'F' || c == 'R'
}
