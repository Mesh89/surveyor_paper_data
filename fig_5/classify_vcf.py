#!/usr/bin/env python3
import argparse
import pysam
import os
from collections import defaultdict

MIN_COV = 2000

def has_large_indel_in_region(read, region_start, region_end, indel_threshold):
    """Check if there is an indel greater than `indel_threshold` in either of the specified regions."""

    added_bp = 0
    ref_pos = read.reference_start + 1  # Convert 0-based to 1-based
    for op in read.cigartuples:
        op_char, op_len = op

        if op_char in {pysam.CINS, pysam.CDEL}:  # 'I' = Insertion, 'D' = Deletion
            indel_start, indel_end = ref_pos, ref_pos+1
            if op_char == pysam.CDEL:
                indel_end += op_len - 1
            if min(indel_end, region_end) - max(indel_start, region_start) > 0:
                if op_char == pysam.CINS:
                    added_bp += op_len
                elif op_char == pysam.CDEL:
                    added_bp -= op_len
        if op_char in {pysam.CMATCH, pysam.CDEL, pysam.CEQUAL, pysam.CDIFF}:  # Match/Mismatch, Deletion, etc., contribute to ref_pos
            ref_pos += op_len

    return abs(added_bp) >= indel_threshold

def overlaps_breakpoint(aln, inv_start, inv_end):
    """Check if the read overlaps with the inversion breakpoints."""
    read_inv_overlap = max(0, min(aln.reference_end, inv_end) - max(aln.reference_start, inv_start))
    lf_overlap = max(0, inv_start - aln.reference_start)
    rf_overlap = max(0, aln.reference_end - inv_end)
    flank_overlap = lf_overlap + rf_overlap
    return flank_overlap >= MIN_COV and read_inv_overlap >= min(inv_end - inv_start, MIN_COV)

def calc_min_diff(read, contig_len, min_diff_fraction):
    """Calculate the minimum difference in AS score as a fraction of the inversion size."""
    ref_split = read.reference_name.split('_')
    lf_size, rf_size = int(ref_split[-2]), int(ref_split[-1])
    inv_start, inv_end = lf_size, contig_len - rf_size
    read_inv_overlap = max(0, min(read.reference_end, inv_end) - max(read.reference_start, inv_start))
    lf_overlap = max(0, inv_start - read.reference_start)
    rf_overlap = max(0, read.reference_end - inv_end)
    flank_overlap = lf_overlap + rf_overlap
    return min_diff_fraction * min(read_inv_overlap, flank_overlap)

def main():
    parser = argparse.ArgumentParser(
        description="Classify inversions."
    )
    parser.add_argument("alt_bam", help="Path to ALT BAM file")
    parser.add_argument("ref_bam", help="Path to REF BAM file")
    parser.add_argument("full_ref_bam", help="Path to FULL-REF BAM file")
    parser.add_argument("inv_vcf", help="Path to VCF file with inversion calls")
    parser.add_argument("output_vcf_t1", type=str, help="VCF.gz file to output supported T1 inversions with at least 3 reads")
    parser.add_argument("output_vcf_t2", type=str, help="VCF.gz file to output supported T2 inversions with at least 3 reads")
    parser.add_argument("output_vcf_fp", type=str, help="VCF.gz file to output unsupported inversions")
    parser.add_argument("--min_diff_fraction", type=float, default=0.5, help="Minimum difference in AS score as a fraction of the inversion size")
    parser.add_argument("--max-indel-dist", type=int, default=100, help="Maximum distance from breakpoint to consider for indels for T2 reads.")
    parser.add_argument("--indel_threshold", type=int, default=50, help="Minimum indel size to classify a read as T2 (default: 10 bp)")
    args = parser.parse_args()
    
    inv_ids = []
    with pysam.VariantFile(args.inv_vcf, 'r') as vcf:
        for record in vcf:
            inv_ids.append(record.id)

    inversion_support_reads = {
        inv_id: {'T1_LB': list(), 'T1_RB': list(), 'T2_LB': list(), 'T2_RB': list()}
        for inv_id in inv_ids
    }
    inversion_deny_reads = {
        inv_id: {'LB': list(), 'RB': list()}
        for inv_id in inv_ids
    }

    # Surprisingly, some reads have higher AS in REF than in FULL-REF, despite REF being a subset of FULL-REF. 
    # Problem with minimap2 heuristic or am I missing something?
    # For this reason, we check for the max between the two files
    best_as_ref, best_contig_ref = defaultdict(int), dict()
    with pysam.AlignmentFile(args.full_ref_bam, "rb") as full_ref_bam:
        for aln in full_ref_bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue

            ref_bam_as = aln.get_tag("AS")
            best_as_ref[aln.qname] = max(best_as_ref[aln.qname], ref_bam_as)

    with pysam.AlignmentFile(args.ref_bam, "rb") as ref_bam:
        for aln in ref_bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue

            ref_bam_as = aln.get_tag("AS")
            if ref_bam_as >= best_as_ref[aln.qname]:
                best_as_ref[aln.qname] = ref_bam_as
                best_contig_ref[aln.qname] = aln.reference_name

    best_as_alt = defaultdict(int)
    alt_accept_reads, ref_accept_reads = list(), list()
    alt_bam = pysam.AlignmentFile(args.alt_bam, "rb")
    for aln in alt_bam.fetch(until_eof=True):
        if aln.is_unmapped:
            continue

        as_score = aln.get_tag("AS")
        best_as_alt[aln.qname] = max(best_as_alt[aln.qname], as_score)

        contig_len = alt_bam.get_reference_length(aln.reference_name)
        ref_split = aln.reference_name.split('_')
        lf_size, rf_size = int(ref_split[-2]), int(ref_split[-1])
        inv_start, inv_end = lf_size, contig_len - rf_size
        
        if overlaps_breakpoint(aln, inv_start, inv_end):
            min_diff = calc_min_diff(aln, contig_len, args.min_diff_fraction)
            if as_score - best_as_ref[aln.qname] >= min_diff:
                alt_accept_reads.append(aln)

    with pysam.AlignmentFile(args.ref_bam, "rb") as ref_bam:
        for aln in ref_bam.fetch(until_eof=True):
            if aln.is_unmapped:
                continue

            as_score = aln.get_tag("AS")
            contig_len = ref_bam.get_reference_length(aln.reference_name)
            ref_split = aln.reference_name.split('_')
            lf_size, rf_size = int(ref_split[-2]), int(ref_split[-1])
            inv_start, inv_end = lf_size, contig_len - rf_size

            if overlaps_breakpoint(aln, inv_start, inv_end):
                min_diff = calc_min_diff(aln, contig_len, args.min_diff_fraction)
                if as_score >= best_as_ref[aln.qname] and as_score - best_as_alt[aln.qname] >= min_diff:
                    ref_accept_reads.append(aln)

    for aln in alt_accept_reads:
        contig_len = alt_bam.get_reference_length(aln.reference_name)
        ref_split = aln.reference_name.split('_')
        lf_size, rf_size = int(ref_split[-2]), int(ref_split[-1])
        inv_start, inv_end = lf_size, contig_len - rf_size
        inv_name = '_'.join(ref_split[:-2])
        tier = "T1"
        lb, rb = False, False
        if aln.reference_start < inv_start-MIN_COV and aln.reference_end > inv_start+MIN_COV:
            lb = True
            if has_large_indel_in_region(aln, inv_start-args.max_indel_dist, inv_start+args.max_indel_dist, args.indel_threshold):
                tier = "T2"
        if aln.reference_start < inv_end-MIN_COV and aln.reference_end > inv_end+MIN_COV:
            rb = True
            if has_large_indel_in_region(aln, inv_end-args.max_indel_dist, inv_end+args.max_indel_dist, args.indel_threshold):
                tier = "T2"

        if lb:
            inversion_support_reads[inv_name][f"{tier}_LB"].append(aln)
        if rb:
            inversion_support_reads[inv_name][f"{tier}_RB"].append(aln)
    
    for aln in ref_accept_reads:
        contig_len = alt_bam.get_reference_length(aln.reference_name)
        ref_split = aln.reference_name.split('_')
        lf_size, rf_size = int(ref_split[-2]), int(ref_split[-1])
        inv_start, inv_end = lf_size, contig_len - rf_size
        
        if aln.reference_start < inv_start-MIN_COV and aln.reference_end > inv_start+MIN_COV:
            lb = True
        if aln.reference_start < inv_end-MIN_COV and aln.reference_end > inv_end+MIN_COV:
            rb = True

        inv_name = '_'.join(ref_split[:-2])
        if lb:
            inversion_deny_reads[inv_name]["LB"].append(aln)
        if rb:
            inversion_deny_reads[inv_name]["RB"].append(aln)

    with pysam.VariantFile(args.inv_vcf, 'r') as vcf_in:
        header = vcf_in.header
        with pysam.VariantFile(args.output_vcf_t1, 'w', header=header) as vcf_out_t1, \
            pysam.VariantFile(args.output_vcf_t2, 'w', header=header) as vcf_out_t2, \
            pysam.VariantFile(args.output_vcf_fp, 'w', header=header) as vcf_out_fp:
            for record in vcf_in:
                inv_id = record.id
                alt_counts = inversion_support_reads[inv_id]
                ref_counts = inversion_deny_reads[inv_id]
                lb_t1, rb_t1 = len(alt_counts['T1_LB']), len(alt_counts['T1_RB'])
                lb_t2, rb_t2 = len(alt_counts['T2_LB']), len(alt_counts['T2_RB'])
                lb_neg, rb_neg = len(ref_counts['LB']), len(ref_counts['RB'])
                print(f"{inv_id}: {lb_t1}/{lb_t2}/{lb_neg} LB, {rb_t1}/{rb_t2}/{rb_neg} RB")
                if lb_t1 >= 3 and lb_t1/(lb_t1+lb_neg) >= 0.25 and rb_t1 >= 3 and rb_t1/(rb_t1+rb_neg) >= 0.25:
                    vcf_out_t1.write(record)
                    outdir = os.path.dirname(args.output_vcf_t1)
                    reads = sorted(alt_counts['T1_LB'] + alt_counts['T1_RB'], key=lambda x: (x.reference_id, x.reference_start))
                    with pysam.AlignmentFile(f"{outdir}/{inv_id}.T1.bam", "wb", header=alt_bam.header) as out_bam:
                        for read in reads:
                            out_bam.write(read)
                    reads = sorted(alt_counts['T2_LB'] + alt_counts['T2_RB'], key=lambda x: (x.reference_id, x.reference_start))
                    with pysam.AlignmentFile(f"{outdir}/{inv_id}.T2.bam", "wb", header=alt_bam.header) as out_bam:
                        for read in reads:
                            out_bam.write(read)
                    reads = sorted(ref_counts['LB'] + ref_counts['RB'], key=lambda x: (x.reference_id, x.reference_start))
                    with pysam.AlignmentFile(f"{outdir}/{inv_id}.ref.bam", "wb", header=alt_bam.header) as out_bam:
                        for read in reads:
                            out_bam.write(read)
                elif lb_t1+lb_t2 >= 3 and (lb_t1+lb_t2)/(lb_t1+lb_t2+lb_neg) >= 0.25 and rb_t1+rb_t2 >= 3 and (rb_t1+rb_t2)/(rb_t1+rb_t2+rb_neg) >= 0.25:
                    vcf_out_t2.write(record)
                else:
                    vcf_out_fp.write(record)
    alt_bam.close()

if __name__ == "__main__":
    main()
