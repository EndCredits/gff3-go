package main

import (
	"encoding/json"
	"flag"
	"fmt"
	"io"
	"os"

	"gff3-go"
)

type DirectiveSummary struct {
	Kind string   `json:"kind"`
	Args []string `json:"args"`
}

type Stats struct {
	File         string             `json:"file"`
	TotalRecords int                `json:"total_records"`
	TypeCounts   map[string]int     `json:"type_counts"`
	SourceCounts map[string]int     `json:"source_counts"`
	StrandCounts map[string]int     `json:"strand_counts"`
	Directives   []DirectiveSummary `json:"directives"`
	SeqIDs       int                `json:"unique_seqids"`
	Errors       int                `json:"errors"`
	ErrorMsgs    []string           `json:"error_messages,omitempty"`
}

func main() {
	flag.Parse()
	path := flag.Arg(0)
	if path == "" {
		fmt.Fprintln(os.Stderr, "usage: gff3stat <file.gff3>")
		os.Exit(1)
	}

	f, err := os.Open(path)
	if err != nil {
		fmt.Fprintf(os.Stderr, "open: %v\n", err)
		os.Exit(1)
	}
	defer f.Close()

	stats := Stats{
		File:         path,
		TypeCounts:   make(map[string]int),
		SourceCounts: make(map[string]int),
		StrandCounts: make(map[string]int),
	}

	r := gff3.NewReader(f)
	seqIDs := make(map[string]bool)

	for {
		rec, err := r.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			stats.Errors++
			stats.ErrorMsgs = append(stats.ErrorMsgs, err.Error())
			continue
		}
		if rec == nil {
			break
		}
		stats.TotalRecords++
		stats.TypeCounts[rec.Type]++
		stats.SourceCounts[rec.Source]++
		stats.StrandCounts[rec.Strand]++
		seqIDs[rec.SeqID] = true
	}

	stats.SeqIDs = len(seqIDs)

	for _, d := range r.Directives() {
		stats.Directives = append(stats.Directives, DirectiveSummary{
			Kind: d.Kind.String(),
			Args: d.Args,
		})
	}

	enc := json.NewEncoder(os.Stdout)
	enc.SetIndent("", "  ")
	enc.Encode(stats)
}

