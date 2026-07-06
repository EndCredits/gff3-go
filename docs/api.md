# API Reference

Import the main module:

```go
import "github.com/EndCredits/gff3-go"
```

Import the binary index module (Unix only):

```go
import "github.com/EndCredits/gff3-go/gff3idx"
```

## Module paths

| Module | Import path |
|--------|------------|
| Main (parser/writer) | `github.com/EndCredits/gff3-go` |
| Binary index | `github.com/EndCredits/gff3-go/gff3idx` |

These are two separate Go modules with versioned dependencies.

---

## Types

### Record

A single GFF3 feature line (9 columns).

```go
type Record struct {
    SeqID      string
    Source     string
    Type       string
    Start      int
    End        int
    Score      string     // "." or numeric string
    Strand     string     // "+", "-", ".", "?"
    Phase      int        // 0, 1, 2, or PhaseUndefined (-1)
    Attributes Attributes
}
```

**Methods:**

| Method | Description |
|--------|-------------|
| `Unmarshal(line string) error` | Parse a tab-delimited GFF3 line into the record |
| `Marshal() (string, error)` | Serialize back to a GFF3 line (validates first) |
| `Validate() error` | Check: non-empty seqid (no spaces, no `>` prefix), start ≤ end, coords positive, valid strand, CDS requires defined phase, source and type non-empty |
| `Clone() *Record` | Deep copy (including Attributes) |

### Attributes

Parsed column 9 tag=value pairs.

```go
type Attributes map[string][]string
```

**Methods:**

| Method | Description |
|--------|-------------|
| `Get(tag string) string` | First value for the tag, or `""` if not present |
| `Clone() Attributes` | Deep copy of all tags and values |

Multi-value tags (`Parent`, `Alias`, `Note`, `Dbxref`, `Ontology_term`) store comma-separated values as separate entries in the slice. Other tags treat commas as literal content.

### Reader

Streaming GFF3 parser wrapping an `io.Reader` via `bufio.Scanner`.

```go
r := gff3.NewReader(file)
for {
    rec, err := r.Read()
    if err == io.EOF { break }
    // use rec
}
dirs := r.Directives()  // collected ## lines
```

**Methods:**

| Method | Description |
|--------|-------------|
| `Read() (*Record, error)` | Next feature record. Returns `io.EOF` at end of file, on `##FASTA`, on `###` (terminator), or on implicit `>` |
| `ReadFASTA() (*FastaRecord, error)` | Next FASTA sequence. Call after `Read()` returns `io.EOF` due to FASTA section. Handles multiple records in sequence |
| `Directives() []Directive` | All `##` directives collected during parsing |

Skips blank lines and `#` comments. The first non-blank, non-comment line must be `##gff-version` — returns an error (wrapped with line number) otherwise.

### Writer

Serializes records and directives back to GFF3 text.

```go
w := gff3.NewWriter(file)
w.WriteDirective(gff3.Directive{Kind: gff3.DirGFFVersion, Args: []string{"3"}})
w.WriteRecord(rec)
```

**Methods:**

| Method | Description |
|--------|-------------|
| `WriteDirective(Directive) error` | Write a `##` directive line. Recognizes all known kinds and `DirUnknown` |
| `WriteRecord(*Record) error` | Write a feature line (validates before writing) |

### Directive

A `##`-prefixed meta-data line.

```go
type Directive struct {
    Kind DirectiveKind
    Args []string
}
```

For `DirUnknown` directives, the first element of `Args` is the unknown keyword itself (preserved verbatim).

### DirectiveKind

```go
type DirectiveKind int

// Kinds:
DirGFFVersion       // ##gff-version
DirSequenceRegion   // ##sequence-region
DirFeatureOntology  // ##feature-ontology
DirAttributeOntology // ##attribute-ontology
DirSourceOntology   // ##source-ontology
DirSpecies          // ##species
DirGenomeBuild      // ##genome-build
DirTerminator       // ###
DirFASTA            // ##FASTA
DirUnknown          // any unrecognized ## keyword
```

`DirectiveKind.String()` converts a kind to its canonical keyword string (e.g., `DirGFFVersion` → `"gff-version"`, `DirTerminator` → `"###"`).

### LineType

```go
type LineType int

const (
    LineDirective LineType = iota  // ## line
    LineComment                     // # line
    LineFeature                     // feature data line
    LineBlank                       // empty line
    LineFASTA                       // > sequence header
)
```

### FastaRecord

A single FASTA sequence record.

```go
type FastaRecord struct {
    ID          string
    Description string
    Sequence    []byte   // upper-cased, whitespace stripped
}

func (f *FastaRecord) SeqString() string  // Sequence as string
```

### Target

Parsed `Target` attribute for alignment features. Format: `target_id start end [strand]`.

```go
type Target struct {
    ID     string
    Start  int
    End    int
    Strand string   // "+" or "-" (defaults to "+" if omitted)
}
```

### GapOp

A single CIGAR-style gap operation (e.g. `M8`, `D3`).

```go
type GapOp struct {
    Code   byte    // M, I, D, F, R
    Length int
}
```

- `M` — match
- `I` — insert gap into reference
- `D` — delete from reference (gap in target)
- `F` — frameshift forward
- `R` — frameshift reverse

### SequenceRegion

Parsed `##sequence-region` directive.

```go
type SequenceRegion struct {
    SeqID string
    Start int
    End   int
}
```

Format: `##sequence-region seqid start end`. Requires at least 3 arguments.

### CycleError

Returned by `DetectCycle()` when a circular Parent/ID relationship is found.

```go
type CycleError struct {
    Node string   // the ID where the cycle was detected
}

func (e *CycleError) Error() string
```

---

## Functions

### Parsing

| Function | Signature | Description |
|----------|-----------|-------------|
| `ParseAttributes` | `func(s string) (Attributes, error)` | Parse a column 9 attributes string into a tag→values map. Splits on `;`, then each pair on `=`. Multi-value tags split on `,`. Percent-decoded after splitting. |
| `ParseDirective` | `func(line string) (Directive, error)` | Parse a `##`-prefixed line. Returns `DirTerminator` for `###`, `DirFASTA` for `##FASTA`. Unknown keywords get `DirUnknown` with keyword as first arg. |
| `ParseSequenceRegion` | `func(d Directive) (SequenceRegion, error)` | Parse `##sequence-region` directive. Requires at least 3 args: seqid, start, end. |
| `ParseTarget` | `func(s string) (Target, error)` | Parse a Target attribute value. Requires at least 3 space-separated fields. Strand defaults to `"+"`. |
| `ParseGap` | `func(s string) ([]GapOp, error)` | Parse a Gap attribute value. Each segment is `[code][length]`, space-separated. Valid codes: M, I, D, F, R. |

### Escape

| Function | Signature | Description |
|----------|-----------|-------------|
| `Escape` | `func(s string) string` | Encode file-level reserved chars: `\t`, `\n`, `\r`, `%`, control chars (0x00–0x1F, 0x7F). Short-circuits if no escaping needed. |
| `EscapeAttr` | `func(s string) string` | Encode file-level + column-9 reserved chars: additionally escapes `;`, `=`, `&`, `,`. Short-circuits if no escaping needed. |
| `Unescape` | `func(s string) string` | Decode Percent-Encoding (`%XX`). Does **not** treat `+` as space. Short-circuits if no `%` found. |

### Utilities

| Function | Signature | Description |
|----------|-----------|-------------|
| `ReadAllFASTA` | `func(r io.Reader) ([]*FastaRecord, error)` | Parse all FASTA sequences from a standalone FASTA reader. Whitespace stripped, bases upper-cased. |
| `GroupByID` | `func(records []*Record) map[string][]*Record` | Group records by ID attribute. Records without an ID are excluded. Each group preserves all segments as-is. |
| `DetectCycle` | `func(records []*Record) error` | Check for circular Parent/ID relationships. Returns `*CycleError` on failure, `nil` otherwise. Uses DFS with visited set. |

### Constructors

| Function | Signature | Description |
|----------|-----------|-------------|
| `NewReader` | `func(r io.Reader) *Reader` | Create a GFF3 reader using `bufio.NewScanner` |
| `NewWriter` | `func(w io.Writer) *Writer` | Create a GFF3 writer |

---

## Constants

### Strand

```go
StrandPlus    = "+"
StrandMinus   = "-"
StrandNone    = "."
StrandUnknown = "?"
```

### Phase

```go
PhaseUndefined = -1   // sentinel for features without a phase (non-CDS)
```

### Directive kinds

```go
DirGFFVersion, DirSequenceRegion, DirFeatureOntology,
DirAttributeOntology, DirSourceOntology, DirSpecies,
DirGenomeBuild, DirTerminator, DirFASTA, DirUnknown
```

### Line types

```go
LineDirective, LineComment, LineFeature, LineBlank, LineFASTA
```

---

## Binary index (`gff3idx`, Unix only)

The binary index is a persistent, mmap-based lookup structure built from parsed GFF3 records. Import:

```go
import "github.com/EndCredits/gff3-go/gff3idx"
```

**Dependencies:** `zeebo/xxh3` (hashing), `golang.org/x/sys/unix` (mmap). No cgo. Platform: Linux, macOS, and other Unix systems.

### Build

```go
func Build(records []*gff3.Record, outPath string) error
```

Builds a binary index file from parsed GFF3 records. Internally:
- Merges discontiguous feature extents per ID
- Builds gene parent→child hierarchies (gene → mRNA → CDS/exon)
- Collects spatial feature data per chromosome
- Constructs xxh3-based open-addressing hash tables
- Writes all sections to file

### Open

```go
func Open(path string) (*Reader, error)
```

Mmaps the index file read-only. Validates the magic number (`GFFI`), version (1), and all section bounds against the file size. Returns an error if any check fails. Caller must call `Close()` to unmap.

### Wrap (in-memory)

```go
func Wrap(records []*gff3.Record) *MemQuerier
```

Creates an in-memory querier from parsed records. No build step or file I/O. Applies the same extent merging and gene hierarchy logic as `Build`. Suitable for ad-hoc queries without persisting an index.

### Querier interface

Both `*Reader` (binary) and `*MemQuerier` (in-memory) implement this interface:

```go
type Querier interface {
    ByID(id string) (*Feature, bool)
    ChildrenOf(geneID string) (*GeneChildren, bool)
    InRange(chr string, minStart, maxEnd int) []SpatialFeat
}
```

Code that accepts `gff3idx.Querier` works with either backend.

### Feature

Returned by `ByID()`. Represents a single feature with merged extent coordinates.

```go
type Feature struct {
    SeqID  string
    Source string
    Type   string
    Start  int     // min start across all segments with this ID
    End    int     // max end across all segments with this ID
    Score  string
    Strand string
    Phase  int
}
```

### GeneChildren

Returned by `ChildrenOf()`. Contains sorted ID lists of a gene's child features.

```go
type GeneChildren struct {
    Transcripts []string  // mRNA IDs (two-level: gene → mRNA)
    CDSs        []string  // CDS feature IDs
    Exons       []string  // exon feature IDs
}
```

Children are discovered by traversing the Parent hierarchy. mRNA children of a gene are direct; CDS and exon children can be either direct children of the gene or grandchildren via mRNA. All lists are sorted alphabetically by ID.

### SpatialFeat

Returned by `InRange()`. Represents a single segment of a feature within a chromosome range.

```go
type SpatialFeat struct {
    Start int
    End   int
    ID    string
    Type  string
}
```

Results are unordered but bounded by the query range. Features are scanned from the first entry with Start ≥ minStart until Start > maxEnd.

### Reader methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `Close` | `func() error` | Unmap the index and close the file descriptor |
| `EntryCount` | `func() uint32` | Number of indexed feature entries |
| `GeneCount` | `func() uint32` | Number of indexed genes |
| `ChrCount` | `func() uint32` | Number of chromosomes in the spatial index |
| `ByID` | `func(id string) (*Feature, bool)` | O(1) hash lookup by feature ID |
| `ChildrenOf` | `func(geneID string) (*GeneChildren, bool)` | O(1) hash lookup for gene child features |
| `InRange` | `func(chr string, minStart, maxEnd int) []SpatialFeat` | Binary search + linear scan for features in a chromosome range |

### MemQuerier methods

| Method | Signature | Description |
|--------|-----------|-------------|
| `ByID` | `func(id string) (*Feature, bool)` | O(1) map lookup by feature ID |
| `ChildrenOf` | `func(geneID string) (*GeneChildren, bool)` | O(1) map lookup for gene child features |
| `InRange` | `func(chr string, minStart, maxEnd int) []SpatialFeat` | Binary search over sorted slice |

### Format constants (gff3idx/format.go)

```go
Magic   = "GFFI"          // 4-byte magic
Version = 1               // format version
```

These are also exported as functions in `gff3idx`:

```go
func ByteOrder() binary.ByteOrder  // binary.LittleEndian
func NextPow2(n uint32) uint32     // smallest power of 2 ≥ n
```

### Tools

| Command | Path | Description |
|---------|------|-------------|
| `gff3stat` | `cmd/gff3stat/main.go` | Parse a GFF3 file → JSON statistics (record counts by type/source/strand, unique seqids, directives). Import: `github.com/EndCredits/gff3-go` |
| `gff3index` | `gff3idx/cmd/gff3index/main.go` | Build a binary index from a GFF3 file. Usage: `gff3index <input.gff3> <output.gff3idx>` |
| `gff3verify` | `gff3idx/cmd/gff3verify/main.go` | Verify a binary index against an in-memory reference. Builds both, compares all entries, gene children, and spatial queries. Usage: `gff3verify <input.gff3>` |
