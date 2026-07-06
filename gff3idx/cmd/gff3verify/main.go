package main

import (
	"fmt"
	"os"
	"sort"

	"github.com/EndCredits/gff3-go"
	"github.com/EndCredits/gff3-go/gff3idx"
)

func main() {
	if len(os.Args) != 2 {
		fmt.Fprintln(os.Stderr, "Usage: gff3verify <input.gff3>")
		os.Exit(1)
	}
	if err := run(os.Args[1]); err != nil {
		fmt.Fprintf(os.Stderr, "FAILED: %v\n", err)
		os.Exit(1)
	}
}

func run(gff3Path string) error {
	fmt.Fprintf(os.Stderr, "Parsing %s...\n", gff3Path)
	f, err := os.Open(gff3Path)
	if err != nil {
		return fmt.Errorf("open: %w", err)
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
	fmt.Fprintf(os.Stderr, "Parsed %d records\n", len(records))

	fmt.Fprintf(os.Stderr, "Building reference index...\n")
	seen := make(map[string]*idRecord)
	parentToChildren := make(map[string][]string)
	for _, rec := range records {
		id := rec.Attributes.Get("ID")
		if id == "" {
			continue
		}
		start := int64(rec.Start)
		end := int64(rec.End)
		if prev, ok := seen[id]; ok {
			if start < prev.minStart {
				prev.minStart = start
			}
			if end > prev.maxEnd {
				prev.maxEnd = end
			}
		} else {
			seen[id] = &idRecord{r: rec, minStart: start, maxEnd: end}
		}
		for _, pid := range rec.Attributes["Parent"] {
			parentToChildren[pid] = append(parentToChildren[pid], id)
		}
	}

	recordByID := make(map[string]*gff3.Record)
	for id, ir := range seen {
		r2 := ir.r.Clone()
		r2.Start = int(ir.minStart)
		r2.End = int(ir.maxEnd)
		recordByID[id] = r2
	}

	type refGeneChildren struct {
		transcripts []string
		cdss        []string
		exons       []string
	}
	refGenes := make(map[string]*refGeneChildren)
	for geneID, geneRec := range recordByID {
		if geneRec.Type != "gene" {
			continue
		}
		txSet := make(map[string]bool)
		cdsSet := make(map[string]bool)
		exonSet := make(map[string]bool)

		for _, childID := range parentToChildren[geneID] {
			child, ok := recordByID[childID]
			if !ok {
				continue
			}
			switch child.Type {
			case "mRNA":
				txSet[childID] = true
				for _, grandchildID := range parentToChildren[childID] {
					gc, ok2 := recordByID[grandchildID]
					if !ok2 {
						continue
					}
					switch gc.Type {
					case "CDS":
						cdsSet[grandchildID] = true
					case "exon":
						exonSet[grandchildID] = true
					}
				}
			case "CDS":
				cdsSet[childID] = true
			case "exon":
				exonSet[childID] = true
			}
		}

		refGenes[geneID] = &refGeneChildren{
			transcripts: sortedKeys(txSet),
			cdss:        sortedKeys(cdsSet),
			exons:       sortedKeys(exonSet),
		}
	}

	chrMap := make(map[string][]refSpatialFeat)
	for _, rec := range records {
		id := rec.Attributes.Get("ID")
		if id == "" {
			continue
		}
		chrMap[rec.SeqID] = append(chrMap[rec.SeqID], refSpatialFeat{
			Start: rec.Start,
			End:   rec.End,
			ID:    id,
			Type:  rec.Type,
		})
	}
	chrKeys := make([]string, 0, len(chrMap))
	for k := range chrMap {
		chrKeys = append(chrKeys, k)
		sort.Slice(chrMap[k], func(i, j int) bool { return chrMap[k][i].Start < chrMap[k][j].Start })
	}
	sort.Strings(chrKeys)

	binPath := gff3Path + ".verify.idx"
	fmt.Fprintf(os.Stderr, "Building binary index at %s...\n", binPath)
	if err := gff3idx.Build(records, binPath); err != nil {
		os.Remove(binPath)
		return fmt.Errorf("build: %w", err)
	}
	defer os.Remove(binPath)

	fmt.Fprintf(os.Stderr, "Opening binary index...\n")
	idx, err := gff3idx.Open(binPath)
	if err != nil {
		return fmt.Errorf("open: %w", err)
	}
	defer idx.Close()

	failures := 0

	fmt.Fprintf(os.Stderr, "Comparing entries (%d in reference)...\n", len(recordByID))
	ids := make([]string, 0, len(recordByID))
	for id := range recordByID {
		ids = append(ids, id)
	}
	sort.Strings(ids)

	entryChecked := 0
	for _, id := range ids {
		ref := recordByID[id]
		feat, ok := idx.ByID(id)
		if !ok {
			fmt.Printf("MISS  entry  %s (binary missing)\n", id)
			failures++
			continue
		}
		if ref.SeqID != feat.SeqID || ref.Source != feat.Source || ref.Type != feat.Type ||
			ref.Start != feat.Start || ref.End != feat.End ||
			ref.Score != feat.Score || ref.Strand != feat.Strand || ref.Phase != feat.Phase {
			fmt.Printf("DIFF  entry  %s: ref={%s %s %s %d-%d %s %s %d} bin={%s %s %s %d-%d %s %s %d}\n",
				id,
				ref.SeqID, ref.Source, ref.Type, ref.Start, ref.End, ref.Score, ref.Strand, ref.Phase,
				feat.SeqID, feat.Source, feat.Type, feat.Start, feat.End, feat.Score, feat.Strand, feat.Phase)
			failures++
		}
		entryChecked++
	}
	fmt.Fprintf(os.Stderr, "  Compared %d entries\n", entryChecked)

	fmt.Fprintf(os.Stderr, "Comparing gene children (%d genes in reference)...\n", len(refGenes))
	geneChecked := 0
	for geneID, ref := range refGenes {
		gc, ok := idx.ChildrenOf(geneID)
		if !ok {
			fmt.Printf("MISS  gene %s (binary missing)\n", geneID)
			failures++
			continue
		}
		if !strSliceEq(ref.transcripts, gc.Transcripts) {
			fmt.Printf("DIFF  gene %s transcripts: ref=%v bin=%v\n", geneID, ref.transcripts, gc.Transcripts)
			failures++
		}
		if !strSliceEq(ref.cdss, gc.CDSs) {
			fmt.Printf("DIFF  gene %s CDSs: ref=%v bin=%v\n", geneID, ref.cdss, gc.CDSs)
			failures++
		}
		if !strSliceEq(ref.exons, gc.Exons) {
			fmt.Printf("DIFF  gene %s exons: ref=%v bin=%v\n", geneID, ref.exons, gc.Exons)
			failures++
		}
		geneChecked++
	}
	fmt.Fprintf(os.Stderr, "  Compared %d genes\n", geneChecked)

	fmt.Fprintf(os.Stderr, "Comparing spatial queries...\n")
	sampleChrs := chrKeys
	if len(sampleChrs) > 10 {
		sampleChrs = sampleChrs[:10]
	}

	spatialChecked := 0
	for _, chr := range sampleChrs {
		refFeats := chrMap[chr]
		if len(refFeats) == 0 {
			continue
		}

		minStart := refFeats[0].Start
		maxEnd := refFeats[len(refFeats)-1].End
		span := maxEnd - minStart
		if span <= 0 {
			continue
		}

		step := span / 4
		if step <= 0 {
			step = 1
		}

		for rng := 0; rng < 3; rng++ {
			qStart := minStart + rng*step
			qEnd := qStart + step - 1
			if qEnd > maxEnd {
				qEnd = maxEnd
			}

			lo := sortSearchSpatialRef(refFeats, qStart)
			var refResult []string
			for j := lo; j < len(refFeats); j++ {
				if refFeats[j].Start > qEnd {
					break
				}
				refResult = append(refResult, refFeats[j].ID)
			}
			sort.Strings(refResult)

			binResult := idx.InRange(chr, qStart, qEnd)
			binIDs := make([]string, len(binResult))
			for i, f := range binResult {
				binIDs[i] = f.ID
			}
			sort.Strings(binIDs)

			if !strSliceEq(refResult, binIDs) {
				fmt.Printf("DIFF  spatial %s [%d-%d]: ref=%d features bin=%d features\n", chr, qStart, qEnd, len(refResult), len(binIDs))
				refSet := make(map[string]bool)
				for _, id := range refResult {
					refSet[id] = true
				}
				for _, id := range binIDs {
					if !refSet[id] {
						fmt.Printf("  extra in bin: %s\n", id)
					}
				}
				binSet := make(map[string]bool)
				for _, id := range binIDs {
					binSet[id] = true
				}
				for _, id := range refResult {
					if !binSet[id] {
						fmt.Printf("  missing from bin: %s\n", id)
					}
				}
				failures++
			}
			spatialChecked++
		}
	}
	fmt.Fprintf(os.Stderr, "  Compared %d spatial queries\n", spatialChecked)

	if failures == 0 {
		fmt.Println("VERIFIED: GFF3 and binary index produce identical results.")
		return nil
	}
	return fmt.Errorf("%d mismatches found", failures)
}

type refSpatialFeat struct {
	Start int
	End   int
	ID    string
	Type  string
}

func sortSearchSpatialRef(recs []refSpatialFeat, target int) int {
	lo, hi := 0, len(recs)-1
	for lo <= hi {
		mid := (lo + hi) / 2
		if recs[mid].Start < target {
			lo = mid + 1
		} else {
			hi = mid - 1
		}
	}
	return lo
}

func sortedKeys(m map[string]bool) []string {
	result := make([]string, 0, len(m))
	for k := range m {
		result = append(result, k)
	}
	sort.Strings(result)
	return result
}

func strSliceEq(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

type idRecord struct {
	r        *gff3.Record
	minStart int64
	maxEnd   int64
}
