// Package gff3 provides a parser, writer, and utilities for the
// Generic Feature Format Version 3 (GFF3) file format.
//
// GFF3 files are nine-column, tab-delimited plain text files used
// in bioinformatics to represent genomic features.
//
// Basic usage:
//
//	r := gff3.NewReader(file)
//	for {
//	    rec, err := r.Read()
//	    if err == io.EOF {
//	        break
//	    }
//	    fmt.Println(rec.SeqID, rec.Type, rec.Start, rec.End)
//	}
package gff3

import (
	"gff3-go/internal/gff3"
)

// Record represents a single GFF3 feature line (9 columns).
type Record = gff3.Record

// Attributes holds the parsed tag=value pairs from column 9.
type Attributes = gff3.Attributes

// Directive represents a ##-prefixed meta-data line.
type Directive = gff3.Directive

// DirectiveKind identifies the type of a ## directive line.
type DirectiveKind = gff3.DirectiveKind

// LineType distinguishes between kinds of GFF3 lines.
type LineType = gff3.LineType

// SequenceRegion holds the parsed ##sequence-region directive.
type SequenceRegion = gff3.SequenceRegion

// Target holds the parsed Target attribute for alignment features.
type Target = gff3.Target

// GapOp is a single operation in a Gap attribute (e.g. M8, D3).
type GapOp = gff3.GapOp

// FastaRecord is a single FASTA sequence record.
type FastaRecord = gff3.FastaRecord

// CycleError indicates a circular Parent/ID relationship.
type CycleError = gff3.CycleError

// Reader reads GFF3 feature records from an io.Reader.
type Reader = gff3.Reader

// Writer writes GFF3 formatted records and directives to an io.Writer.
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

// PhaseUndefined is the sentinel value for features without a phase.
const PhaseUndefined = gff3.PhaseUndefined

// Parser functions.
var (
	ParseAttributes    = gff3.ParseAttributes
	ParseDirective     = gff3.ParseDirective
	ParseSequenceRegion = gff3.ParseSequenceRegion
	ParseTarget        = gff3.ParseTarget
	ParseGap           = gff3.ParseGap
)

// Escape functions.
var (
	Escape     = gff3.Escape
	EscapeAttr = gff3.EscapeAttr
	Unescape   = gff3.Unescape
)

// Utility functions.
var (
	ReadAllFASTA = gff3.ReadAllFASTA
	GroupByID    = gff3.GroupByID
	DetectCycle  = gff3.DetectCycle
)

// Constructor functions.
var (
	NewReader = gff3.NewReader
	NewWriter = gff3.NewWriter
)
