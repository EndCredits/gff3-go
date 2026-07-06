package gff3

import (
	"bufio"
	"fmt"
	"io"
	"strings"
)

type FastaRecord struct {
	ID          string
	Description string
	Sequence    []byte
}

func (f *FastaRecord) SeqString() string {
	return string(f.Sequence)
}

func (r *Reader) ReadFASTA() (*FastaRecord, error) {
	for r.sc.Scan() {
		r.lineNum++
		line := r.sc.Text()

		if line == "" {
			continue
		}

		if line[0] == '>' {
			id, desc := parseFastaHeader(line[1:])
			var seq []byte
			for r.sc.Scan() {
				r.lineNum++
				next := r.sc.Text()
				if next == "" {
					continue
				}
				if next[0] == '>' {
					rec := &FastaRecord{ID: id, Description: desc, Sequence: seq}
					id, desc = parseFastaHeader(next[1:])
					seq = nil
					return rec, nil
				}
				seq = append(seq, cleanSequence(next)...)
			}
			if id != "" && len(seq) > 0 {
				return &FastaRecord{ID: id, Description: desc, Sequence: seq}, nil
			}
			return nil, io.EOF
		}

		seq := cleanSequence(line)
		if len(seq) > 0 {
			return &FastaRecord{ID: "", Sequence: seq}, nil
		}
	}

	if err := r.sc.Err(); err != nil {
		return nil, err
	}
	return nil, io.EOF
}

func parseFastaHeader(s string) (id, desc string) {
	s = strings.TrimSpace(s)
	idx := strings.IndexAny(s, " \t")
	if idx < 0 {
		return s, ""
	}
	return s[:idx], strings.TrimSpace(s[idx+1:])
}

func cleanSequence(s string) []byte {
	var out []byte
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c == ' ' || c == '\t' || c == '\r' || c == '\n' {
			continue
		}
		if c >= 'a' && c <= 'z' {
			c -= 32
		}
		out = append(out, c)
	}
	return out
}

func ReadAllFASTA(r io.Reader) ([]*FastaRecord, error) {
	sc := bufio.NewScanner(r)
	var records []*FastaRecord
	var current *FastaRecord

	for sc.Scan() {
		line := sc.Text()
		if line == "" {
			continue
		}
		if line[0] == '>' {
			if current != nil && len(current.Sequence) > 0 {
				records = append(records, current)
			}
			id, desc := parseFastaHeader(line[1:])
			current = &FastaRecord{ID: id, Description: desc}
			continue
		}
		if current == nil {
			continue
		}
		current.Sequence = append(current.Sequence, cleanSequence(line)...)
	}

	if current != nil && len(current.Sequence) > 0 {
		records = append(records, current)
	}

	if err := sc.Err(); err != nil {
		return nil, fmt.Errorf("gff3: fasta read error: %w", err)
	}
	return records, nil
}
