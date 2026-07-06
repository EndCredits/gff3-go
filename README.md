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
| *A. hypogaea* genome | 215 MB | 983,853 | 1.08s | 199 MB/s |

Reproduce with:

```bash
go run ./cmd/gff3stat/ annotations.gff3 | jq .total_records
```

Micro-benchmarks:

```bash
go test -bench=. -benchmem ./internal/gff3/
```

## Documentation

- [Design overview](docs/design.md)
- [API reference](docs/api.md)
- [Package docs (pkg.go.dev)](https://pkg.go.dev/gff3-go)

## License

MIT
