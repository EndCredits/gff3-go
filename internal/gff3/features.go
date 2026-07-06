package gff3

func GroupByID(records []*Record) map[string][]*Record {
	groups := make(map[string][]*Record)
	for _, r := range records {
		id := r.Attributes.Get("ID")
		if id == "" {
			continue
		}
		groups[id] = append(groups[id], r)
	}
	return groups
}

func DetectCycle(records []*Record) error {
	children := make(map[string][]string)
	for _, r := range records {
		id := r.Attributes.Get("ID")
		if id == "" {
			continue
		}
		for _, parent := range r.Attributes["Parent"] {
			children[parent] = append(children[parent], id)
		}
	}

	for id := range children {
		visited := make(map[string]bool)
		if hasCycle(children, id, visited) {
			return &CycleError{Node: id}
		}
	}
	return nil
}

type CycleError struct {
	Node string
}

func (e *CycleError) Error() string {
	return "gff3: parent cycle detected involving " + e.Node
}

func hasCycle(graph map[string][]string, node string, visiting map[string]bool) bool {
	if visiting[node] {
		return true
	}
	visiting[node] = true
	for _, child := range graph[node] {
		if hasCycle(graph, child, visiting) {
			return true
		}
	}
	visiting[node] = false
	return false
}
