# Design

## Architecture

```
gff3.go              Public API (type aliases, single import surface)
└── internal/gff3/   Implementation (not importable by external packages)
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
         ├── Unescape()           (col 1,2,3,5,6: % decoding)
         ├── parseInt()           (col 4,5: start, end)
         ├── parsePhase()         (col 8: 0|1|2|undefined)
         └── ParseAttributes()    (col 9)
                ├── split by ';'  → tag=value pairs
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

`EscapeAttr()` is used when writing attribute values. `Escape()` is sufficient for non-attribute fields and tag names.

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
                └── EOF ──→ return io.EOF          │
                                                    │
  inFASTA=true ──→ ReadFASTA() returns *FastaRecord │
```

## Writer round-trip

The Writer produces output that, when re-parsed, yields identical records. Key design decisions:

- **Attribute order**: tags are sorted alphabetically for deterministic output
- **Multi-value tags**: values joined with `,`, each value escaped with `EscapeAttr()`
- **Single-value tags**: each value pair emitted separately, escaped with `EscapeAttr()`
- **Phase**: emitted as `0`/`1`/`2` or `.` for `PhaseUndefined`
