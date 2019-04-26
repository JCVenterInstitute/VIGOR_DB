#!/usr/bin/env python
"""
Compare genes, CDS's and features from two GFF3 files
"""
from __future__ import print_function
import sys
import argparse
import itertools

# overlap to consider two elements to be the same element
MINIMUM_COVERAGE = .5

# GFF3 field constants
SEQID = 0
SOURCE = 1
TYPE = 2
START = 3
END = 4
SCORE = 5
STRAND = 6
PHASE=7
ATTRIBUTES = 8

# GFF3 attributes
A_PARENT = "Parent"
A_ID = "ID"

# parsed entity constants
E_SEQNAME=0
E_ID = 1
E_TYPE = 2
E_RANGE=3
E_PHASE=4
E_STRAND=5
E_ATTRIBUTES = 6
E_LINE = 7
E_SOURCE = 8
E_LINENUMBER = 9

# tree constants
T_SOURCE = "source"
T_ID = "by_id"
T_ENTITIES = "entities"
T_PARENT = "by_parent"
T_CHILD = "by_child"
T_TYPE = "by_type"
T_SEQNAME = "by_seqname"

TYPE_ORDER = dict(gene=1, CDS=2)

class Difference:
    ERROR  = "ERROR"
    WARNING="WARNING"

    def __init__(self, message, level, *entries):
        self.message = message
        self.level = level
        self.entries = entries

    @staticmethod
    def error(message, *entries):
        return Difference(message, Difference.ERROR, *entries)

    @staticmethod
    def warning(message, *entries):
        return Difference(message, Difference.WARNING, *entries)

    def __str__(self):
        msg = ["{}: {}".format(self.level, self.message)]
        for e in self.entries:
            msg.append("{}:{}".format(e[E_SOURCE],e[E_LINENUMBER]))
            msg.append(e[E_LINE])

        return "\n".join(msg)


def phase(phase):
    if phase in ("", "."):
        return 0
    return int(phase)

def parse_attributes(data):
    attribs = dict()
    for kv in data.split(";"):
        t = kv.split("=")
        attribs[t[0]] = t[1] if len(t) > 1 else ""
    return attribs

def get_entry(line, linenumber=-1, source="notset"):
    data = line.split("\t",8)
    try:
        attributes = parse_attributes(data[ATTRIBUTES])
        return (data[SEQID],
                attributes.get(A_ID),
                data[TYPE],
                (int(data[START]), int(data[END])),
                phase(data[PHASE]),
                data[STRAND],
                attributes,
                line,
                source,
                linenumber)
    except (KeyError,IndexError) as e:
        print("ERROR: {}:{} formatting error in line {}".format( source, linenumber, line))
        raise e

def entries(infile, source="notset"):
    structure = {}
    for (linenumber, line) in enumerate(infile):
        if line.startswith("#"):
            continue
        line = line.strip()
        if len(line) == 0:
            continue
        entry = get_entry(line,linenumber, source)
        yield entry

def build_tree(infile, source="notset"):
    tree = {
        T_SOURCE:source,
        T_ENTITIES:[],
        T_ID:dict(),
        T_TYPE:{},
        T_PARENT:{},
        T_CHILD:{},
        T_SEQNAME:{}
    }

    for entry in entries(infile, source=source):
        tree[T_ENTITIES].append(entry)
        tree[T_SEQNAME].setdefault(entry[E_SEQNAME],[]).append(entry)
        e_id = entry[E_ID]
        tree[T_ID][e_id] = entry
        tree[T_TYPE].setdefault(entry[E_SEQNAME],set()).add( (entry[E_TYPE],entry[E_RANGE]) )
        parent =  entry[E_ATTRIBUTES].get(A_PARENT)
        if parent:
            tree[T_PARENT][e_id] = tree[T_ID][parent]
            tree[T_CHILD].setdefault(parent,[]).append(entry)
    return tree

def check_phase( e1,e2):
    if e1[E_PHASE] != e2[E_PHASE]:
        return '{} {} phase {} != phase {}'.format( e1[E_SEQNAME], e1[E_TYPE], e1[E_PHASE], e2[E_PHASE])

def check_strand(e1,e2):
    if e1[E_STRAND] != e2[E_STRAND]:
        return '{} {} strand {} != strand {}'.format( e1[E_SEQNAME], e1[E_TYPE], e1[E_STRAND], e2[E_STRAND])

def check_range(e1,e2):
    if e1[E_RANGE] != e2[E_RANGE]:
        intersection = range_intersection(e1[E_RANGE], e2[E_RANGE])
        return "{} {} range difference ({}) != ({}) overlap ({})".format(e1[E_SEQNAME],
                                                                         e1[E_TYPE],
                                                                         format_range(e1[E_RANGE]),
                                                                         format_range(e2[E_RANGE]),
                                                                         format_range(intersection))
    
def parse_args():
    parser = argparse.ArgumentParser(description="Compare features from two gff3 files")
    parser.add_argument("gff1",metavar="FIRST", help="first GFF3 file to compare")
    parser.add_argument("gff2",metavar="SECOND", help="second GFF3 file to compare")
    parser.add_argument("-c","--coverage", metavar="COVERAGE",
                      help="percent range (0-1) overlap of an element for another element of the same type  to be considered the same",
                      default=.5,
                      dest="coverage",
                      type=float
    )
    return parser.parse_args()

def range_coverage(erange, intersection):
    if intersection:
        return float (intersection[1] - intersection[0]) / float(erange[1] - erange[0])
    return 0

def format_range(erange):
    return "{}-{}/{}".format(erange[0], erange[1],erange[1]-erange[0])

def find_corresponding_entry(entry1, entries):
    entry_range = entry1[E_RANGE]

    entry = (0,None)
    for (coverage,entry2) in ( (range_coverage(entry_range, range_intersection(entry_range, g[E_RANGE])), g) for g in entries):
        if MINIMUM_COVERAGE < coverage and entry[0] < coverage:
            entry = (coverage, entry2)

    return entry[1]

def compare_entries(e1, e2, comparisons):
    messages = []
    for comparison in comparisons:
        msg = comparison(e1,e2)
        if msg:
            messages.append(msg)
    if messages:
        yield "\n".join(messages)

def compare_trees(tree1, tree2):
    """
    Compare two GFF trees. Yields Difference objects
    """
    source1 = tree1[T_SOURCE]
    source2 = tree2[T_SOURCE]

    for (seqname,entries1) in ( (key, ( e for e in value if e[E_TYPE] == 'gene')) for (key,value) in tree1[T_SEQNAME].iteritems()):
        
        entries2 = [e for e in tree2[T_SEQNAME].get(seqname,[]) if e[E_TYPE] == "gene"]
        if not entries2:
            yield Difference.error("{} no genes in {}".format(seqname, source2))
            continue

        for entry1 in (e for e in entries1 if e[E_TYPE] == "gene"):
            e_type = entry1[E_TYPE]
            e_range = entry1[E_RANGE]
            e_phase = entry1[E_PHASE]
            e_strand = entry1[E_STRAND]

            entry2 = find_corresponding_entry(entry1, entries2)
            if not entry2:
                yield Difference.error("{} no {} with range ({})".format(seqname,  e_type,format_range(e_range)), entry1)
                continue

            for message in compare_entries(entry1, entry2, [check_range, check_strand]):
                yield Difference.error(message, entry1, entry2)

            # compare CDS
            # Note: the current reference GFF3 contain duplicate IDs so also filter on sequence name
            cdss = [ c for c in tree2[T_CHILD].get(entry2[E_ID],[]) if c[E_TYPE] == 'CDS' and c[E_SEQNAME] == seqname]

            # Note: the current reference GFF3 contain duplicate IDs so also filter on sequence name
            for cds1 in (c for c in tree1[T_CHILD].get(entry1[E_ID],[]) if c[E_TYPE] == 'CDS' and c[E_SEQNAME] == seqname):
                cds2 = find_corresponding_entry(cds1, cdss)
                if not cds2:
                    yield Difference.error("{} no corresponding CDS found {}".format(cds1[E_SEQNAME], format_range(cds1[E_RANGE])), cds1)
                    continue

                for message in compare_entries(cds1,cds2, [check_range, check_strand, check_phase]):
                    yield Difference.error(message, cds1, cds2)

def range_intersection(range1, range2):
    """
    return overlap, if any of two tuples of (start, end)
    returns overlap as (start, end) or empty tuple of no overlap
    """

    if range1[0] <= range2[1]:
        if range1[1] >= range2[0]:
            overlap =  (
                max(range1[0], range2[0]),
                min(range1[1], range2[1])
            )
            if overlap[0] == overlap[1]:
                return tuple()
            return overlap
        else:
            return tuple()
    else:
        return range_intersection(range2, range1)

def difference_signature(difference):
    eids = [ str(e) for e in sorted([id(e) for e in difference.entries]) ]
    return ";".join(eids)

if __name__ == '__main__':
    args = parse_args()
    tree1 = build_tree(open(args.gff1), args.gff1)
    tree2 = build_tree(open(args.gff2), args.gff2)
    MINIMUM_COVERAGE = args.coverage
    ok = True
    difference_filter = set()
    for difference in itertools.chain(
            compare_trees(tree1,tree2),
            compare_trees(tree2,tree1)):
        d_signature = difference_signature(difference)
        if d_signature in difference_filter:
            continue
        difference_filter.add(d_signature)
        print(difference)
        print()
        ok &= difference.level != Difference.ERROR
    sys.exit(0 if ok else 1)
