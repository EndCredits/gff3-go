package main

import (
	"fmt"
	"os"

	"github.com/EndCredits/gff3-go"
	"github.com/EndCredits/gff3-go/gff3idx"
)

func main() {
	if len(os.Args) != 3 {
		fmt.Fprintf(os.Stderr, "Usage: gff3index <input.gff3> <output.gff3idx>\n")
		os.Exit(1)
	}

	f, err := os.Open(os.Args[1])
	if err != nil {
		fmt.Fprintf(os.Stderr, "open: %v\n", err)
		os.Exit(1)
	}
	defer f.Close()

	r := gff3.NewReader(f)
	var records []*gff3.Record
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		if rec == nil {
			break
		}
		records = append(records, rec)
	}

	if err := gff3idx.Build(records, os.Args[2]); err != nil {
		fmt.Fprintf(os.Stderr, "build: %v\n", err)
		os.Exit(1)
	}

	fmt.Printf("Indexed %d records → %s\n", len(records), os.Args[2])
}
