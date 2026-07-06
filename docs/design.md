# Design

## Architecture

```
gff3.go              Public API (type aliases, single import surface)
├── internal/gff3/   GFF3 parser + writer (module: gff3-go, zero deps)
└── gff3idx/         Binary index — mmap, O(1) lookup (module: github.com/EndCredits/gff3-go/gff3idx)
```

`gff3.go` re-exports all public symbols from `internal/gff3/` via Go type aliases (`type Record = gff3.Record`). This keeps the API surface visible in one file while the implementation is organized across multiple internal files.

## Parsing pipeline

```
io.Reader
  → bufio.Scanner (line-by-line)
    → line type detection
      ├── ## → ParseDirective() → Directive
      ├── #  → skip (comment)
      ├── >  → enter FASTA mode → ReadFASTA()
└── feature line → Record.Unmarshal()
         ├── splitColumns()       (tab → 9 fields)
         ├── Unescape()           (col 1,2,3,6,7: % decoding)
         ├── parseInt()           (col 4,5: start, end)
         ├── parsePhase()         (col 8: 0|1|2|undefined)
         └── ParseAttributes()    (col 9)
                ├── single pass inline (no intermediate pair slice)
                ├── split by '='  → tag, value
                ├── Unescape()    → decoded tag/value
                └── split by ','  → multi-value (Parent, Alias, Note...)
```

## Escape layers

GFF3 has two layers of Percent-Encoding:

| Layer | Characters | Function |
|-------|-----------|----------|
| File-level | tab, newline, CR, `%`, control chars | `Escape()` |
| Column 9 | `;` `=` `&` `,` (in addition to file-level) | `EscapeAttr()` |

`EscapeAttr()` is used for all attribute tag names and values. The `Unescape()` function is a custom RFC 3986 decoder that does **not** treat `+` as space (the spec declares `+` encoding deprecated).

## Reader state machine

```
                    ┌─────────────────────────────┐
                    │                             │
  Start ──→ Scanning ──→ Read() returns *Record   │
                │                                  │
                ├── ## (directive) ──→ collect     │
                ├── #  (comment)   ──→ skip        │
                ├── >  (implied FASTA) ──→ FASTA   │
                ├── ##FASTA ──→ inFASTA=true ──→ EOF
                ├── ###  ──→ terminated=true ──→ EOF
                └── EOF ──→ return io.EOF          │
                                                    │
  First non-blank/non-comment line MUST be         │
  ##gff-version. Returns error otherwise.           │
```

## Writer round-trip

The Writer produces output that, when re-parsed, yields identical records. Key design decisions:

- **Attribute order**: tags are sorted alphabetically for deterministic output
- **Multi-value tags**: values joined with `,`, each value escaped with `EscapeAttr()`
- **Single-value tags**: each value pair emitted separately, escaped with `EscapeAttr()`
- **Phase**: emitted as `0`/`1`/`2` or `.` for `PhaseUndefined`

## Discontiguous features

GFF3 allows a single feature to span multiple lines with the same ID (e.g., a CDS split across exons). The parser returns each segment as a separate `Record`. Downstream consumers must merge them:

```go
groups := gff3.GroupByID(records)
// groups["cds00001"] = [Record{1201..1500}, Record{3000..3902}, ...]
```

The binary index merges extent coordinates (min start, max end) into a single `EntryRecord` per ID. Individual segments are preserved in the spatial index for accurate interval queries.

## Binary index (`gff3idx`)

### File layout

```
┌─ Header (64 bytes)
│    Magic "GFFI", Version 1, counts, offsets to all sections
├─ Entry Hash Table
│    Open-addressing, linear probe, load factor 50%
│    Hash: xxh3 (zeebo/xxh3), non-zero remapped (0→1 to avoid sentinel collision)
│    HashSlot: {Hash uint64, Val uint64}
│    Val: upper 32 = record index, lower 32 = string pool offset of ID
├─ Entry Records (40 bytes each, sorted by ID)
│    {Start int64, End int64, ChrOffset, SourceOffset, TypeOffset,
│     ScoreOffset, StrandOffset uint32, Phase int32}
├─ Gene Hash Table (same design as entries)
├─ Gene Records (24 bytes each, sorted by gene ID)
│    {TranscriptCount, CDSCount, ExonCount, DataOffset}
├─ Gene Data
│    For each gene: arrays of uint32 string-pool offsets for
│    transcript IDs, CDS IDs, exon IDs
├─ Spatial Index
│    Per-chromosome: {ChrOffset, FeatureCount, DataOffset}
│    Per-feature: {Start int64, End int64, IDOffset, TypeOffset}
│    Features sorted by Start within each chromosome
└─ String Pool
     All strings deduplicated, null-terminated, concatenated
     Empty strings → sentinel 0xFFFFFFFF
```

### Query API

Two backends, one interface:

```go
type Querier interface {
    ByID(id string) (*Feature, bool)
    ChildrenOf(geneID string) (*GeneChildren, bool)
    InRange(chr string, minStart, maxEnd int) []SpatialFeat
}
```

**In-memory** (no build step, ~500MB for 983K records):

```go
q := gff3idx.Wrap(records)       // *MemQuerier
feat, _ := q.ByID("Ah01g000200")
```

**Binary index** (persistent, mmap, ~50MB resident):

```go
idx, _ := gff3idx.Open("genes.gff3idx")
defer idx.Close()

feat, _ := idx.ByID("Ah01g000200")        // O(1) hash lookup
children, _ := idx.ChildrenOf("Ah01g000200") // gene → transcripts/CDSs/exons
feats := idx.InRange("chr01", 1_000_000, 2_000_000) // binary search + scan
```

Both implement `Querier`, so downstream code can accept the interface and work with either backend.

### Build

Programmatic:

```go
gff3idx.Build(records, "genes.gff3idx")
```

Or via CLI:

```bash
go run ./gff3idx/cmd/gff3index/ genes.gff3 genes.gff3idx
```

The builder:
1. Collects entries (merging discontiguous feature extents), genes (two-level hierarchy traversal), and spatial features (individual segments)
2. Builds string pool with deduplication
3. Builds open-addressing hash tables (xxh3, 50% load factor, linear probe)
4. Writes all sections in order, finally writing the header

### Design decisions

- **int64 coordinates**: entry and spatial records use 64-bit coordinates for future-proofing against large synthetic genomes
- **mmap**: the entire index is mapped read-only via `unix.Mmap`. No heap allocations for index data. Platform: Unix (Linux, macOS)
- **bounds validation**: `validateBounds()` checks all section offsets against file size on open. Individual record accesses also check before `unsafe.Pointer` casts
- **hash sentinel**: `hashNonZero()` remaps xxh3 hash 0→1 so the zero slot unambiguously means "empty"
- **spatial dedup**: spatial index stores raw segment coordinates per GFF3 line, not merged extents. This ensures accurate interval queries for discontiguous features
- **empty string sentinel**: `0xFFFFFFFF` in string pool offsets means empty string. The string pool never reaches 4 GB in practice
