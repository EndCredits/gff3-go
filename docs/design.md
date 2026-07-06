# Design

## Architecture

The project consists of two Go modules:

```
gff3.go                   Public API (type aliases, single import surface)
├── internal/gff3/        GFF3 parser + writer (module: github.com/EndCredits/gff3-go, zero deps)
└── gff3idx/              Binary index — mmap, O(1) lookup (module: github.com/EndCredits/gff3-go/gff3idx)
```

`gff3.go` re-exports all public symbols from `internal/gff3/` via Go type aliases (`type Record = gff3.Record`) and variable references (`var NewReader = gff3.NewReader`). This keeps the API surface visible in one file while the implementation is organized across multiple internal files.

`gff3idx/` is a **separate Go module** with its own `go.mod` (`github.com/EndCredits/gff3-go/gff3idx`). It depends on the main module via a versioned `require` (`github.com/EndCredits/gff3-go v0.1.0`).

## Parsing pipeline

```
io.Reader
  → bufio.Scanner (line-by-line)
    → line type detection (by first character)
      ├── "##" → ParseDirective() → Directive
      ├── "#"  → skip (comment)
      ├── ">"  → enter FASTA mode → ReadFASTA()
      └── else → Record.Unmarshal()
            ├── splitColumns()      (tab → 9 fields)
            ├── Unescape()          (cols 1,2,3,6,7: % decoding)
            ├── parseInt()          (cols 4,5: start, end)
            ├── parsePhase()        (col 8: 0|1|2|undefined)
            └── ParseAttributes()   (col 9)
                  ├── split by ';' → tag=value pairs
                  ├── split by '=' → tag, value
                  ├── for multi-value: split by ',' → then Unescape each
                  └── for single-value: Unescape value
```

Column index (1-based) and processing:

| Col | Field      | Processing              |
|-----|------------|-------------------------|
| 1   | seqid      | `Unescape()`            |
| 2   | source     | `Unescape()`            |
| 3   | type       | `Unescape()`            |
| 4   | start      | `parseInt()` (integer)  |
| 5   | end        | `parseInt()` (integer)  |
| 6   | score      | `Unescape()`            |
| 7   | strand     | `Unescape()`            |
| 8   | phase      | `parsePhase()`          |
| 9   | attributes | `ParseAttributes()`     |

## Escape layers

GFF3 has two layers of Percent-Encoding:

| Layer    | Characters                               | Function       |
|----------|------------------------------------------|----------------|
| File-level | `\t`, `\n`, `\r`, `%`, control chars (0x00–0x1F, 0x7F) | `Escape()`     |
| Column 9   | `;` `=` `&` `,` (in addition to file-level) | `EscapeAttr()` |

`Escape()` is used for columns 1,2,3,6 (seqid, source, type, score). `EscapeAttr()` is used for all attribute tag names and values in column 9. The `Unescape()` function is a custom RFC 3986 percent-decoder that does **not** treat `+` as space (the GFF3 spec declares `+` encoding deprecated).

Both `Escape()` and `EscapeAttr()` short-circuit: if no character in the string needs escaping, the input is returned as-is without allocation. When escaping is needed, a new byte slice is allocated and each reserved byte is encoded as `%XX` using uppercase hex.

## Writer round-trip

The Writer produces output that, when re-parsed, yields identical records. Key design decisions:

- **SeqID, Source, Type, Score**: escaped with `Escape()` (file-level)
- **Strand**: written directly (always one of `+`, `-`, `.`, `?`)
- **Phase**: written as `0`/`1`/`2` or `.` for `PhaseUndefined`
- **Attribute tags**: sorted alphabetically for deterministic output
- **Multi-value tags** (`Parent`, `Alias`, `Note`, `Dbxref`, `Ontology_term`): values joined with `,`, each value individually escaped with `EscapeAttr()`
- **Single-value tags**: each `tag=value` pair emitted as a separate semicolon-delimited segment, value escaped with `EscapeAttr()`
- **Empty attributes**: written as `.`
- `WriteRecord()` calls `Marshal()`, which calls `Validate()` first — invalid records are rejected before any bytes are written

## Reader state machine

```
                     ┌──────────────────────────────┐
                     │                              │
  Start ──→ Scanning ──→ Read() returns *Record     │
                 │                                  │
                 ├── ##gff-version → sawVersion=true │
                 ├── ## (other)    → collect directive│
                 ├── #  (comment)  → skip             │
                 ├── >  (implied FASTA) → inFASTA=true│
                 ├── ##FASTA → inFASTA=true → EOF     │
                 ├── ###  → terminated=true → EOF     │
                 ├── blank line → skip                │
                 └── EOF → return io.EOF              │
                                                      │
  The first non-blank, non-comment line MUST be a     │
  ##gff-version directive. Returns error otherwise.   │
```

After `inFASTA` or `terminated` is set, `Read()` immediately returns `io.EOF`. Subsequent calls to `ReadFASTA()` parse FASTA sequences (which may contain `>` headers). `ReadFASTA()` handles both explicit `##FASTA` and implicit `>` entries, including multiple FASTA records in sequence.

The `Directives()` method returns all `##` directives collected during parsing, available at any time.

## Discontiguous features

GFF3 allows a single feature to span multiple lines with the same ID (e.g., a CDS split across exons). The parser returns each segment as a separate `Record`. Downstream consumers must merge them:

```go
groups := gff3.GroupByID(records)
// groups["cds00001"] = [Record{1201..1500}, Record{3000..3902}, ...]
```

**Merge strategies differ by context:**

- **`GroupByID`**: groups records by ID, preserving each segment as-is. Callers process segments individually.
- **Binary index `ByID`**: merges extent coordinates (min Start across all segments with same ID, max End) into a single `EntryRecord`. Returns one `Feature` per ID.
- **`InRange` (spatial queries)**: uses individual segment coordinates from `SpatialFeatureRec`, not merged extents. Each GFF3 line with an ID becomes a separate entry in the spatial index. This ensures accurate interval queries — a query overlapping only part of a discontiguous feature returns only the matching segments.

`MemQuerier` applies the same logic: `ByID` returns merged extents, `InRange` returns individual segments.

## Binary index (`gff3idx`, Unix only)

The binary index is a mmap-based persistent file built from parsed GFF3 records. It provides O(1) by-ID lookup, O(1) gene-children lookup, and O(log n + k) spatial queries.

**Dependencies:** `zeebo/xxh3` (hashing), `golang.org/x/sys/unix` (mmap). Does not use cgo.

### File layout

```
┌─ Header (64 bytes) ─────────────────────────────────┐
│  Magic "GFFI" (4 bytes), Version=1, EntryCount,       │
│  SpatialChrs, GeneCount, StringPoolSize               │
│  Offsets: EntriesOffset, SpatialOffset,               │
│           GenesOffset, StringPoolOff                  │
├─ Entry Hash Table ───────────────────────────────────┤
│  Open-addressing, linear probe, load factor 50%       │
│  Hash: xxh3.HashString → hashNonZero (0→1 remap)     │
│  HashSlot: {Hash uint64, Val uint64}                 │
│    Val: upper 32 bit = record index                  │
│         lower 32 bit = string pool offset of ID       │
│  HashSlotSize = 16 bytes                             │
├─ Entry Records (40 bytes each, sorted by ID) ────────┤
│  {Start int64, End int64,                             │
│   ChrOffset, SourceOffset, TypeOffset,               │
│   ScoreOffset, StrandOffset uint32,                   │
│   Phase int32}                                        │
│  All string fields stored as offsets into string pool │
├─ Gene Hash Table (same design as entry table) ───────┤
├─ Gene Records (24 bytes each, sorted by gene ID) ───┤
│  {TranscriptCount, CDSCount, ExonCount uint32,       │
│   _ uint32, DataOffset uint64}                       │
├─ Gene Data ──────────────────────────────────────────┤
│  Per gene: concatenated uint32 arrays of              │
│  string-pool offsets for:                             │
│    transcript IDs, CDS IDs, exon IDs                  │
│  DataOffset in GeneRecord points here                 │
├─ Spatial Index ──────────────────────────────────────┤
│  uint32 chr_count                                    │
│  Per chromosome: SpatialHeader                       │
│    {ChrOffset uint32, FeatureCount uint32,           │
│     DataOffset uint64}                                │
│  Per feature: SpatialFeatureRec                      │
│    {Start int64, End int64,                          │
│     IDOffset uint32, TypeOffset uint32}              │
│  Features sorted by Start within each chromosome     │
└─ String Pool ────────────────────────────────────────┘
     All strings deduplicated, null-terminated,
     concatenated in order of first occurrence.
     Empty strings → sentinel 0xFFFFFFFF.
     StringAt reads from given offset to next null byte.
```

### Hash table design

- **Hash function**: xxh3 (128-bit hash via `zeebo/xxh3.HashString`), truncated to uint64. This is a fast non-cryptographic hash with excellent distribution.
- **hashNonZero**: remaps hash value 0 → 1. Since 0 is the empty-slot sentinel in the hash table, a valid slot must never have Hash==0. This avoids false-negative probes.
- **Open addressing**: linear probing. On collision, walk `(pos + 1) % capacity` until an empty slot is found.
- **Load factor**: capacity is `nextPow2(count * 2)`, yielding ~50% load factor. The power-of-two table size makes modulo fast (bitmask), and linear probing stays efficient at this density.
- **Lookup**: hash the query ID, probe linearly, compare Hash == stored hash, then verify the ID string from the string pool matches. Stored Val encodes both the record index and the ID string's pool offset for O(1) verification.

### String pool

All strings are stored once (first occurrence), null-terminated, concatenated in the string pool at the end of the file. The pool is referenced by 32-bit offsets. An empty string is represented by the sentinel `0xFFFFFFFF` — since no valid string pool offset can be this large (the pool never reaches 4 GB), the reader returns `""` when it encounters this sentinel.

### mmap and bounds validation

The index file is mapped into memory read-only via `unix.Mmap` (`PROT_READ | MAP_SHARED`). No heap allocations for index data — all reads go directly to the mmap'd region using `unsafe.Pointer` casts and `unsafe.Slice` for safe slices.

On `Open()`, `validateBounds()` checks every section's offset and size against the file size before any data is accessed. Individual record accesses also do bounds checks before `unsafe.Pointer` casts. The `.Close()` method calls `unix.Munmap` and closes the file descriptor.

## MemQuerier vs Reader

Two query backends, one interface:

```go
type Querier interface {
    ByID(id string) (*Feature, bool)
    ChildrenOf(geneID string) (*GeneChildren, bool)
    InRange(chr string, minStart, maxEnd int) []SpatialFeat
}
```

| Aspect | `MemQuerier` | `*Reader` |
|--------|-------------|-----------|
| Construction | `gff3idx.Wrap(records)` | `gff3idx.Open(path)` |
| Storage | In-memory maps (`map[string]*gff3.Record`) | mmap'd binary file |
| Lookup | O(1) hashmap | O(1) hash table (linear probe) |
| Spatial | Binary search over sorted slice | Binary search over mmap'd records |
| Build step | None (uses parsed records) | `gff3idx.Build()` |
| Persistence | No | Yes (file on disk) |
| Destroy on close | GC'd | File remains, `Close()` unmaps |

The `Wrap()` function also performs ID extent merging for `ByID()` and builds the gene parent→child hierarchy, mirroring the builder's logic. Both backends produce identical query results for the same input.

## Build

### Programmatic

```go
gff3idx.Build(records, "genes.gff3idx")
```

### CLI

```bash
go run ./gff3idx/cmd/gff3index/ genes.gff3 genes.gff3idx
```

The builder:

1. **Collect** (`builder.collect`): iterates all records, merging discontiguous feature extents (min Start, max End per ID), building parent→child maps, and collecting spatial feature entries per chromosome.
2. **Build string pool**: each unique string is assigned an offset. Empty strings get sentinel `0xFFFFFFFF`.
3. **Build hash tables**: open-addressing tables for entry IDs and gene IDs, using xxh3 hashes with linear probing.
4. **Write sections** in order: header reservation → entry hash table → entry records → gene hash table → gene records → gene data (transcript/CDS/exon ID arrays) → spatial index → string pool.
5. **Write header**: seek back to offset 0 and write the header with final counts and offsets.

## Design decisions

- **int64 coordinates in binary format**: `EntryRecord` and `SpatialFeatureRec` use `int64` for Start/End, future-proofing against large synthetic genomes. The public `Feature` and `Record` types use `int` (converted at the boundary).
- **hashNonZero**: xxh3 can produce hash value 0, which conflicts with the empty-slot sentinel in the hash table. `hashNonZero` remaps 0→1 so that Hash==0 unambiguously means "empty slot".
- **Empty string sentinel**: `0xFFFFFFFF` as string pool offset means empty string. The pool can never reach 4 GB in practice. The reader's `stringAt()` method's bounds check naturally returns `""` for this sentinel.
- **Spatial segments not merged**: the spatial index stores individual GFF3 line coordinates (one entry per line with an ID), not merged extents. This ensures `InRange()` returns accurate results for discontiguous features — a query overlapping one exon of a multi-exon CDS returns only that exon's entry.
- **Zero external dependencies in core module**: `internal/gff3/` (and thus the main module) has zero third-party dependencies. Only `gff3idx` pulls in `zeebo/xxh3` and `golang.org/x/sys`.
- **Unix-only mmap**: `gff3idx` uses `golang.org/x/sys/unix` for mmap and is therefore restricted to Unix platforms (Linux, macOS). `MemQuerier` has no platform restrictions.
