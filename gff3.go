// Package gff3 provides a parser, writer, and utilities for the
// Generic Feature Format Version 3 (GFF3) file format.
//
// GFF3 is a nine-column, tab-delimited plain text format used in
// bioinformatics to represent genomic features such as genes, exons,
// and coding sequences.
//
// # Basic usage
//
// Read a GFF3 file:
//
//	r := gff3.NewReader(file)
//	for {
//	    rec, err := r.Read()
//	    if err == io.EOF {
//	        break
//	    }
//	    fmt.Printf("%s\t%s\t%d-%d\t%s\n",
//	        rec.SeqID, rec.Type, rec.Start, rec.End, rec.Strand)
//	}
//	for _, d := range r.Directives() {
//	    fmt.Printf("##%s\n", d.Kind)
//	}
//
// Write a GFF3 file:
//
//	w := gff3.NewWriter(file)
//	w.WriteDirective(gff3.Directive{Kind: gff3.DirGFFVersion, Args: []string{"3"}})
//	w.WriteRecord(&gff3.Record{
//	    SeqID:  "chr1",
//	    Source: ".",
//	    Type:   "gene",
//	    Start:  1000,
//	    End:    9000,
//	    Strand: gff3.StrandPlus,
//	    Attributes: gff3.Attributes{"ID": {"gene1"}, "Name": {"EDEN"}},
//	})
//
// # Validation and utilities
//
//	if err := rec.Validate(); err != nil {
//	    log.Fatal(err)
//	}
//	groups := gff3.GroupByID(records)
//	if err := gff3.DetectCycle(records); err != nil {
//	    log.Fatal(err)
//	}
//
// # FASTA section
//
// After the Reader encounters ##FASTA, call ReadFASTA():
//
//	for {
//	    seq, err := r.ReadFASTA()
//	    if err == io.EOF {
//	        break
//	    }
//	    fmt.Printf(">%s\n%s\n", seq.ID, seq.SeqString())
//	}
//
// Or parse standalone FASTA:
//
//	seqs, _ := gff3.ReadAllFASTA(reader)
//
// # Alignment sub-parsers
//
//	target, _ := gff3.ParseTarget(rec.Attributes.Get("Target"))
//	gap, _ := gff3.ParseGap(rec.Attributes.Get("Gap"))
//
// The package has zero external dependencies.
package gff3

import (
	"github.com/EndCredits/gff3-go/internal/gff3"
)

// Record represents a single GFF3 feature line (9 columns).
//
//   Column 1: seqid     — landmark ID (chromosome, scaffold)
//   Column 2: source    — algorithm or database name
//   Column 3: type      — feature type (gene, mRNA, exon, CDS, ...)
//   Column 4: start     — 1-based start coordinate
//   Column 5: end       — 1-based end coordinate
//   Column 6: score     — floating point or "."
//   Column 7: strand    — +, -, ., or ?
//   Column 8: phase     — 0, 1, 2 (CDS), or PhaseUndefined
//   Column 9: attributes — tag=value pairs parsed into Attributes
type Record = gff3.Record

// Attributes holds the parsed tag=value pairs from column 9.
//
// Reserved tags: ID, Name, Alias, Parent, Target, Gap, Derives_from,
// Note, Dbxref, Ontology_term, Is_circular.
//
// Multi-value tags (Parent, Alias, Note, Dbxref, Ontology_term) store
// comma-separated values as individual strings in the slice.
type Attributes = gff3.Attributes

// Directive represents a ##-prefixed pragma line.
//
// The Kind field identifies the directive type. Args contains the
// space-separated arguments following the keyword.
type Directive = gff3.Directive

// DirectiveKind identifies the type of a ## directive line.
type DirectiveKind = gff3.DirectiveKind

// LineType distinguishes between kinds of GFF3 lines.
type LineType = gff3.LineType

// SequenceRegion holds the parsed ##sequence-region directive.
//
// Format: ##sequence-region seqid start end
type SequenceRegion = gff3.SequenceRegion

// Target holds the parsed Target attribute for alignment features.
//
// Format: target_id start end [strand]
type Target = gff3.Target

// GapOp is a single operation in a Gap (CIGAR-style) attribute.
//
//   M — match
//   I — insert gap into reference
//   D — delete from reference (gap in target)
//   F — frameshift forward
//   R — frameshift reverse
type GapOp = gff3.GapOp

// FastaRecord is a single FASTA sequence record.
//
// Sequence bases are stored upper-cased with whitespace removed.
type FastaRecord = gff3.FastaRecord

// CycleError indicates a circular Parent/ID relationship among features.
type CycleError = gff3.CycleError

// Reader reads GFF3 feature records from an io.Reader.
//
// Reader skips blank lines and comments. Directives are collected
// internally and available via Directives(). When a ##FASTA directive
// or an implicit > line is encountered, Read() returns io.EOF and
// subsequent FASTA sequences can be read via ReadFASTA().
//
// The first non-blank, non-comment line must be a ##gff-version
// directive. A ### (terminator) directive stops reading immediately.
type Reader = gff3.Reader

// Writer writes GFF3 formatted records and directives to an io.Writer.
//
// Output is round-trip safe: records written and re-parsed produce
// identical values.
type Writer = gff3.Writer

// Directive kind constants.
const (
	DirGFFVersion       = gff3.DirGFFVersion
	DirSequenceRegion   = gff3.DirSequenceRegion
	DirFeatureOntology  = gff3.DirFeatureOntology
	DirAttributeOntology = gff3.DirAttributeOntology
	DirSourceOntology   = gff3.DirSourceOntology
	DirSpecies          = gff3.DirSpecies
	DirGenomeBuild      = gff3.DirGenomeBuild
	DirTerminator       = gff3.DirTerminator
	DirFASTA            = gff3.DirFASTA
	DirUnknown          = gff3.DirUnknown
)

// Line type constants.
const (
	LineDirective = gff3.LineDirective
	LineComment   = gff3.LineComment
	LineFeature   = gff3.LineFeature
	LineBlank     = gff3.LineBlank
	LineFASTA     = gff3.LineFASTA
)

// Strand constants.
const (
	StrandPlus    = gff3.StrandPlus
	StrandMinus   = gff3.StrandMinus
	StrandNone    = gff3.StrandNone
	StrandUnknown = gff3.StrandUnknown
)

// PhaseUndefined is the sentinel value for features without a phase (non-CDS).
const PhaseUndefined = gff3.PhaseUndefined

// ParseAttributes parses a column 9 attributes string into a tag→values map.
//
// The format is tag=value pairs separated by semicolons. Multiple values for
// the same tag are separated by commas (only for Parent, Alias, Note, Dbxref,
// and Ontology_term). Percent-encoded characters are decoded after splitting
// on reserved delimiters.
var ParseAttributes = gff3.ParseAttributes

// ParseDirective parses a ##-prefixed directive line.
var ParseDirective = gff3.ParseDirective

// ParseSequenceRegion parses a ##sequence-region directive into a SequenceRegion.
var ParseSequenceRegion = gff3.ParseSequenceRegion

// ParseTarget parses a Target attribute value.
//
// Format: target_id start end [strand]. Strand defaults to "+" if omitted.
var ParseTarget = gff3.ParseTarget

// ParseGap parses a Gap attribute value into a slice of GapOp.
//
// Format: space-separated (code,length) pairs, e.g. "M8 D3 M6".
var ParseGap = gff3.ParseGap

// Escape encodes file-level reserved characters using Percent-Encoding.
//
// Escapes: tab, newline, carriage return, %, and control characters.
var Escape = gff3.Escape

// EscapeAttr encodes both file-level and column-9 reserved characters.
//
// In addition to file-level escaping, also escapes ; = & , which have
// reserved meanings in GFF3 column 9. Use this when writing attribute values.
var EscapeAttr = gff3.EscapeAttr

// Unescape decodes GFF3 Percent-Encoding according to RFC 3986.
var Unescape = gff3.Unescape

// ReadAllFASTA reads all FASTA sequences from an io.Reader.
//
// The reader is assumed to contain only FASTA-formatted data
// (lines starting with > followed by sequence lines).
var ReadAllFASTA = gff3.ReadAllFASTA

// GroupByID groups records by their ID attribute.
//
// Records without an ID are excluded. Returns a map from ID to
// all records sharing that ID (discontiguous features).
var GroupByID = gff3.GroupByID

// DetectCycle checks for circular Parent/ID relationships.
//
// Returns CycleError if a cycle is found, nil otherwise.
var DetectCycle = gff3.DetectCycle

// NewReader creates a GFF3 Reader from an io.Reader.
var NewReader = gff3.NewReader

// NewWriter creates a GFF3 Writer writing to an io.Writer.
var NewWriter = gff3.NewWriter
