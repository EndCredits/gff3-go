package gff3idx

import (
	"sort"

	"github.com/EndCredits/gff3-go"
)

type MemQuerier struct {
	byID      map[string]*gff3.Record
	geneToTx  map[string][]string
	geneToCDS map[string][]string
	geneToEx  map[string][]string
	spatial   map[string][]spatialEntry
}

type spatialEntry struct {
	start int
	end   int
	id    string
	typ   string
}

type idExtent struct {
	r        *gff3.Record
	minStart int
	maxEnd   int
}

func Wrap(records []*gff3.Record) *MemQuerier {
	m := &MemQuerier{
		byID:      make(map[string]*gff3.Record),
		geneToTx:  make(map[string][]string),
		geneToCDS: make(map[string][]string),
		geneToEx:  make(map[string][]string),
		spatial:   make(map[string][]spatialEntry),
	}

	parentToChildren := make(map[string][]string)
	extents := make(map[string]*idExtent)

	for _, r := range records {
		id := r.Attributes.Get("ID")
		if id == "" {
			continue
		}
		if prev, ok := extents[id]; ok {
			if r.Start < prev.minStart {
				prev.minStart = r.Start
			}
			if r.End > prev.maxEnd {
				prev.maxEnd = r.End
			}
		} else {
			extents[id] = &idExtent{r: r, minStart: r.Start, maxEnd: r.End}
		}
		m.spatial[r.SeqID] = append(m.spatial[r.SeqID], spatialEntry{
			start: r.Start,
			end:   r.End,
			id:    id,
			typ:   r.Type,
		})
		for _, pid := range r.Attributes["Parent"] {
			parentToChildren[pid] = append(parentToChildren[pid], id)
		}
	}

	for id, ex := range extents {
		r := ex.r.Clone()
		r.Start = ex.minStart
		r.End = ex.maxEnd
		m.byID[id] = r
	}

	for geneID, rec := range m.byID {
		if rec.Type != "gene" {
			continue
		}
		txSet := make(map[string]bool)
		cdsSet := make(map[string]bool)
		exonSet := make(map[string]bool)

		for _, childID := range parentToChildren[geneID] {
			child := m.byID[childID]
			if child == nil {
				continue
			}
			switch child.Type {
			case "mRNA":
				txSet[childID] = true
				for _, grandchildID := range parentToChildren[childID] {
					gc := m.byID[grandchildID]
					if gc == nil {
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

		m.geneToTx[geneID] = sortedSet(txSet)
		m.geneToCDS[geneID] = sortedSet(cdsSet)
		m.geneToEx[geneID] = sortedSet(exonSet)
	}

	for chr := range m.spatial {
		sort.Slice(m.spatial[chr], func(i, j int) bool {
			return m.spatial[chr][i].start < m.spatial[chr][j].start
		})
	}

	return m
}

func (m *MemQuerier) ByID(id string) (*Feature, bool) {
	r, ok := m.byID[id]
	if !ok {
		return nil, false
	}
	return &Feature{
		SeqID:  r.SeqID,
		Source: r.Source,
		Type:   r.Type,
		Start:  r.Start,
		End:    r.End,
		Score:  r.Score,
		Strand: r.Strand,
		Phase:  r.Phase,
	}, true
}

func (m *MemQuerier) ChildrenOf(geneID string) (*GeneChildren, bool) {
	tx, txOk := m.geneToTx[geneID]
	cds, cdsOk := m.geneToCDS[geneID]
	ex, exOk := m.geneToEx[geneID]
	if !txOk && !cdsOk && !exOk {
		return nil, false
	}
	if tx == nil {
		tx = []string{}
	}
	if cds == nil {
		cds = []string{}
	}
	if ex == nil {
		ex = []string{}
	}
	return &GeneChildren{Transcripts: tx, CDSs: cds, Exons: ex}, true
}

func (m *MemQuerier) InRange(chr string, minStart, maxEnd int) []SpatialFeat {
	entries := m.spatial[chr]
	if len(entries) == 0 {
		return nil
	}

	lo := sort.Search(len(entries), func(i int) bool {
		return entries[i].start >= minStart
	})

	var result []SpatialFeat
	for i := lo; i < len(entries); i++ {
		e := entries[i]
		if e.start > maxEnd {
			break
		}
		result = append(result, SpatialFeat{
			Start: e.start,
			End:   e.end,
			ID:    e.id,
			Type:  e.typ,
		})
	}
	return result
}

func sortedSet(m map[string]bool) []string {
	result := make([]string, 0, len(m))
	for k := range m {
		result = append(result, k)
	}
	sort.Strings(result)
	return result
}
