#!/usr/bin/env python3

import sys
import pysam

MITO_CONTIG = "chrM"
DEL_SIZE = 500
GAP_MERGE_WINDOW = 10


def get_merged_deletions(read, gap_merge_window=GAP_MERGE_WINDOW):
    if read.cigartuples is None:
        return []

    ref_pos = read.reference_start
    query_pos = 0
    events = []
    del_start = None
    del_total = 0
    del_query_pos = None

    for op, length in read.cigartuples:
        if op == 2:  # D
            if del_start is None:
                del_start = ref_pos
                del_query_pos = query_pos
            del_total += length
            ref_pos += length
        elif op in (0, 7, 8):  # M, =, X
            if del_total and length > gap_merge_window:
                events.append((del_start, del_total, del_query_pos))
                del_start, del_total, del_query_pos = None, 0, None
            ref_pos += length
            query_pos += length
        elif op == 3:  # N
            ref_pos += length
        elif op == 1:  # I
            query_pos += length
        elif op == 4:  # S
            query_pos += length

    if del_total:
        events.append((del_start, del_total, del_query_pos))
    return events


def extract_large_deletions_to_tsv(input_bam, output_tsv):
    bam_in = pysam.AlignmentFile(input_bam, "rb")

    with open(output_tsv, "w") as out:
        out.write("read_name\tleft_bp\tright_bp\tsize\tquery_pos\n")

        for read in bam_in.fetch(MITO_CONTIG):
            if read.is_unmapped or read.is_secondary or read.cigartuples is None:
                continue

            for left, length, query_pos in get_merged_deletions(read, GAP_MERGE_WINDOW):
                if length >= DEL_SIZE:
                    right = left + length
                    out.write(f"{read.query_name}\t{left}\t{right}\t{length}\t{query_pos}\n")

    bam_in.close()
    return output_tsv


def main():
    input_bam = sys.argv[1]
    deletions_tsv = sys.argv[2]

    extract_large_deletions_to_tsv(input_bam, deletions_tsv)


if __name__ == "__main__":
    main()
