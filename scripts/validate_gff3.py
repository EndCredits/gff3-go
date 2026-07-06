#!/usr/bin/env python3
"""
Cross-validate a GFF3 file using BCBio-GFF parser.

Usage: python3 validate_gff3.py <file.gff3>

Outputs JSON with:
  - total_records
  - type_counts: {feature_type: count}
  - source_counts: {source: count}
  - strand_counts: {strand: count}
  - unique_seqids: int
  - errors: int
"""
import json
import sys
from collections import Counter

def validate(path):
    type_counts = Counter()
    source_counts = Counter()
    strand_counts = Counter()
    seq_ids = set()
    errors = 0
    total = 0
    directives = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("##"):
                directives.append(line)
                continue
            if line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 9:
                errors += 1
                continue
            total += 1
            seq_ids.add(parts[0])
            source_counts[parts[1]] += 1
            type_counts[parts[2]] += 1
            # parts[3] start, parts[4] end
            strand_counts[parts[6]] += 1

    result = {
        "total_records": total,
        "type_counts": dict(type_counts),
        "source_counts": dict(source_counts),
        "strand_counts": dict(strand_counts),
        "unique_seqids": len(seq_ids),
        "errors": errors,
        "directives": directives,
    }
    return result


def validate_bcbio(path):
    """Validate using BCBio-GFF parser."""
    from BCBio import GFF

    type_counts = Counter()
    source_counts = Counter()
    strand_counts = Counter()
    seq_ids = set()
    errors = 0
    total = 0

    with open(path) as f:
        for rec in GFF.parse(f):
            for feature in rec.features:
                _count_feature(feature, rec.id, type_counts, source_counts, strand_counts, seq_ids)
                total += 1
                total += _count_subfeatures(feature.sub_features, rec.id, type_counts, source_counts, strand_counts, seq_ids)

    return {
        "total_records": total,
        "type_counts": dict(type_counts),
        "source_counts": dict(source_counts),
        "strand_counts": dict(strand_counts),
        "unique_seqids": len(seq_ids),
        "errors": errors,
    }


def _count_feature(f, seqid, tc, sc, strc, sids):
    sids.add(seqid)
    tc[f.type] += 1

    src = f.qualifiers.get("source", ["."])
    sc[src[0] if isinstance(src, list) else src] += 1

    s = "."
    try:
        st = f.location.strand
        if st == 1:
            s = "+"
        elif st == -1:
            s = "-"
        else:
            s = "."
    except Exception:
        s = "."
    strc[s] += 1


def _count_subfeatures(subs, seqid, tc, sc, strc, sids):
    total = 0
    for sf in subs:
        _count_feature(sf, seqid, tc, sc, strc, sids)
        total += 1
        total += _count_subfeatures(sf.sub_features, seqid, tc, sc, strc, sids)
    return total


def main():
    if len(sys.argv) < 2:
        print("usage: validate_gff3.py <file.gff3> [--bcbio]", file=sys.stderr)
        sys.exit(1)

    path = sys.argv[1]
    use_bcbio = "--bcbio" in sys.argv

    if use_bcbio:
        result = validate_bcbio(path)
        result["parser"] = "bcbio-gff"
    else:
        result = validate(path)
        result["parser"] = "line-split"

    json.dump(result, sys.stdout, indent=2)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
