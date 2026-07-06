# gff3-go

Go library for parsing and writing [GFF3](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) (Generic Feature Format Version 3) files.

**Zero external dependencies.** Standard library only.

## Install

```bash
go get gff3-go
```

## Quick start

```go
package main

import (
    "fmt"
    "io"
    "os"

    "gff3-go"
)

func main() {
    f, _ := os.Open("annotations.gff3")
    defer f.Close()

    r := gff3.NewReader(f)
    for {
        rec, err := r.Read()
        if err == io.EOF {
            break
        }
        fmt.Printf("%s\t%s\t%d\t%d\n", rec.SeqID, rec.Type, rec.Start, rec.End)
    }

    // Directives collected during parsing
    for _, d := range r.Directives() {
        fmt.Printf("##%s %v\n", d.Kind, d.Args)
    }
}
```

## Features

- **Parse** GFF3 files with a streaming `Reader`
- **Write** GFF3 files with a `Writer` (round-trip safe)
- **Percent-Encoding** for file-level and column-9 reserved characters
- **Attribute parsing**: tag=value pairs, multi-value splitting
- **Sub-parsers**: `Target`, `Gap` (CIGAR-style)
- **Validation**: `Record.Validate()`, `DetectCycle()`
- **FASTA section** parsing
- **CLI tool**: `cmd/gff3stat` for quick statistics

## Performance

Benchmarks on Apple M1 (single core):

| File | Size | Records | Time | Throughput |
|------|------|---------|------|-----------|
| *A. hypogaea* genome | 215 MB | 983,853 | 0.95s | 226 MB/s |

Reproduce with:

```bash
go run ./cmd/gff3stat/ annotations.gff3 | jq .total_records
```

Micro-benchmarks:

```bash
go test -bench=. -benchmem ./internal/gff3/
```

## Validating your GFF3 files

### Quick statistics (CLI)

```bash
go run ./cmd/gff3stat/ your_annotations.gff3
```

Outputs JSON with record counts by type, source, strand, unique seqIDs, and any parse errors:

```json
{
  "total_records": 983853,
  "type_counts": {"gene": 83107, "mRNA": 83107, "exon": 417771, "CDS": 399868},
  "source_counts": {"maker": 237762, "AUGUSTUS": 157012, ...},
  "strand_counts": {"+": 490545, "-": 493308},
  "unique_seqids": 140,
  "errors": 0
}
```

### Validate programmatically

```go
f, _ := os.Open("annotations.gff3")
defer f.Close()

r := gff3.NewReader(f)
var records []*gff3.Record
for {
    rec, err := r.Read()
    if err == io.EOF { break }
    if err != nil { log.Fatal(err) }
    records = append(records, rec)
}

// Check every record
for _, rec := range records {
    if err := rec.Validate(); err != nil {
        log.Printf("invalid record: %v", err)
    }
}

// Check for circular Parent/ID relationships
if err := gff3.DetectCycle(records); err != nil {
    log.Printf("parent cycle: %v", err)
}

// Check discontiguous features
groups := gff3.GroupByID(records)
for id, recs := range groups {
    if len(recs) > 1 && recs[0].Type == "CDS" {
        log.Printf("multi-segment CDS: %s (%d segments)", id, len(recs))
    }
}
```

### Cross-validate with Python

```bash
python3 scripts/validate_gff3.py your_annotations.gff3
```

Compares feature counts, source distribution, and strand balance against our Go parser. Use `--bcbio` for a second independent parser:

```bash
python3 scripts/validate_gff3.py your_annotations.gff3 --bcbio
```

### Round-trip integrity

```bash
go test -run TestRoundTripDeepFile -args -gff3 your_annotations.gff3
```

Parses the first 5000 records, writes them back, re-parses, and verifies all 9 columns plus every attribute value are identical.

### Run unit tests

```bash
go test -cover ./internal/gff3/
```

## Documentation

- [Design overview](docs/design.md)
- [API reference](docs/api.md)
- [Package docs (pkg.go.dev)](https://pkg.go.dev/gff3-go)

## License

MIT
