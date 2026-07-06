# API Reference

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
    Phase      int        // 0, 1, 2, or PhaseUndefined
    Attributes Attributes
}
```

**Methods:**

| Method | Description |
|--------|-------------|
| `Unmarshal(line string) error` | Parse a tab-delimited GFF3 line |
| `Marshal() (string, error)` | Serialize back to GFF3 line (validates first) |
| `Validate() error` | Check required fields and constraints |
| `Clone() *Record` | Deep copy |

### Attributes

Parsed column 9 tag=value pairs.

```go
type Attributes map[string][]string
```

**Methods:**

| Method | Description |
|--------|-------------|
| `Get(tag string) string` | First value for tag, or `""` |
| `Clone() Attributes` | Deep copy |

Multi-value tags (`Parent`, `Alias`, `Note`, `Dbxref`, `Ontology_term`) store comma-separated values as separate slice entries.

### Reader

Streaming GFF3 parser wrapping an `io.Reader`.

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
| `Read() (*Record, error)` | Next feature record |
| `ReadFASTA() (*FastaRecord, error)` | Next FASTA sequence (after `##FASTA`) |
| `Directives() []Directive` | All collected directives |

Automatically skips blank lines and comments. Stops at `##FASTA` or implicit `>seqid`. ReadFASTA() works correctly when preceded by `##FASTA`; with implicit `>`, the header line has already been consumed.

### Writer

Serializes records and directives back to GFF3 text.

```go
w := gff3.NewWriter(file)
w.WriteDirective(d)
w.WriteRecord(rec)
```

**Methods:**

| Method | Description |
|--------|-------------|
| `WriteDirective(Directive) error` | Write a `##` directive line |
| `WriteRecord(*Record) error` | Write a feature line |

### Directive

A `##`-prefixed meta-data line.

```go
type Directive struct {
    Kind DirectiveKind
    Args []string
}
```

**Kinds:** `DirGFFVersion`, `DirSequenceRegion`, `DirFeatureOntology`, `DirAttributeOntology`, `DirSourceOntology`, `DirSpecies`, `DirGenomeBuild`, `DirTerminator`, `DirFASTA`, `DirUnknown`

### FastaRecord

A single FASTA sequence.

```go
type FastaRecord struct {
    ID          string
    Description string
    Sequence    []byte
}

func (f *FastaRecord) SeqString() string
```

### Target

Parsed `Target` attribute for alignment features. Format: `target_id start end [strand]`.

```go
type Target struct {
    ID     string
    Start  int
    End    int
    Strand string
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

### SequenceRegion

Parsed `##sequence-region` directive.

```go
type SequenceRegion struct {
    SeqID string
    Start int
    End   int
}
```

### CycleError

Returned by `DetectCycle` when a circular Parent/ID relationship is found.

---

## Functions

### Parsing

| Function | Description |
|----------|-------------|
| `ParseAttributes(s string) (Attributes, error)` | Parse column 9 attributes |
| `ParseDirective(line string) (Directive, error)` | Parse a `##` directive line |
| `ParseSequenceRegion(Directive) (SequenceRegion, error)` | Parse `##sequence-region` |
| `ParseTarget(s string) (Target, error)` | Parse Target attribute value |
| `ParseGap(s string) ([]GapOp, error)` | Parse Gap attribute value |

### Escape

| Function | Description |
|----------|-------------|
| `Escape(s string) string` | Encode file-level reserved chars |
| `EscapeAttr(s string) string` | Encode file + column-9 reserved chars |
| `Unescape(s string) string` | Decode Percent-Encoding |

### Utilities

| Function | Description |
|----------|-------------|
| `ReadAllFASTA(r io.Reader) ([]*FastaRecord, error)` | Parse all FASTA sequences |
| `GroupByID(records []*Record) map[string][]*Record` | Group records by ID attribute |
| `DetectCycle(records []*Record) error` | Check for circular Parent/ID relationships (returns `*CycleError` on failure) |

### Constructors

| Function | Description |
|----------|-------------|
| `NewReader(r io.Reader) *Reader` | Create a GFF3 reader |
| `NewWriter(w io.Writer) *Writer` | Create a GFF3 writer |

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
PhaseUndefined = -1  // features without a phase (non-CDS)
```

### Directive kinds

`DirGFFVersion`, `DirSequenceRegion`, `DirFeatureOntology`, `DirAttributeOntology`, `DirSourceOntology`, `DirSpecies`, `DirGenomeBuild`, `DirTerminator`, `DirFASTA`, `DirUnknown`

```go
// DirectiveKind.String() converts a kind to its keyword name.
// e.g. DirGFFVersion → "gff-version", DirTerminator → "###"
func (k DirectiveKind) String() string
```

### Line types

```go
LineDirective  // ## line
LineComment    // # line
LineFeature    // feature data line
LineBlank      // empty line
LineFASTA      // > sequence header
```

## Binary index (`gff3idx`, Unix only)

The binary index is a mmap-based persistent index built from parsed records. Import as `gff3-go/gff3idx`.

> **Dependencies:** `zeebo/xxh3`, `golang.org/x/sys/unix` (optional; only needed if you use this sub-package)

### Build

```go
func Build(records []*gff3.Record, outPath string) error
```

Builds a binary index file from parsed GFF3 records.

### Open

```go
func Open(path string) (*Reader, error)
```

Mmaps the index file and returns a queryable Reader.

### Reader methods

```go
func (r *Reader) Close() error
func (r *Reader) EntryCount() uint32
func (r *Reader) GeneCount() uint32
func (r *Reader) ChrCount() uint32

func (r *Reader) ByID(id string) (*Feature, bool)
func (r *Reader) ChildrenOf(geneID string) (*GeneChildren, bool)
func (r *Reader) InRange(chr string, minStart, maxEnd int) []SpatialFeat
```

### Feature

```go
type Feature struct {
    SeqID  string
    Source string
    Type   string
    Start  int
    End    int
    Score  string
    Strand string
    Phase  int
}
```

### GeneChildren

```go
type GeneChildren struct {
    Transcripts []string  // mRNA IDs
    CDSs        []string  // CDS feature IDs
    Exons       []string  // exon feature IDs
}
```

### SpatialFeat

```go
type SpatialFeat struct {
    Start int
    End   int
    ID    string
    Type  string
}
```

### MemQuerier (in-memory)

Zero build cost. Accepts parsed records directly via `Wrap`.

```go
func Wrap(records []*gff3.Record) *MemQuerier

func (m *MemQuerier) ByID(id string) (*Feature, bool)
func (m *MemQuerier) ChildrenOf(geneID string) (*GeneChildren, bool)
func (m *MemQuerier) InRange(chr string, minStart, maxEnd int) []SpatialFeat
```

### Querier interface

Both `MemQuerier` and `*Reader` satisfy:

```go
type Querier interface {
    ByID(id string) (*Feature, bool)
    ChildrenOf(geneID string) (*GeneChildren, bool)
    InRange(chr string, minStart, maxEnd int) []SpatialFeat
}
```

### Tools

| Command | Description |
|---------|-------------|
| `cmd/gff3stat` | Parse GFF3 → JSON statistics |
| `gff3idx/cmd/gff3index` | Build binary index from GFF3 |
| `gff3idx/cmd/gff3verify` | Verify binary index against in-memory reference |
