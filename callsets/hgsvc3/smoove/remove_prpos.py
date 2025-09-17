#!/usr/bin/env python3
import sys
import pysam

if len(sys.argv) != 3:
    sys.stderr.write(f"Usage: {sys.argv[0]} <in.vcf.gz> <out.vcf.gz>\n")
    sys.exit(1)

in_vcf = sys.argv[1]
out_vcf = sys.argv[2]

with pysam.VariantFile(in_vcf, "r") as invcf:
    with pysam.VariantFile(out_vcf, "wz", header=invcf.header) as outvcf:
        for rec in invcf:
            # Remove PRPOS and PREND if present
            for tag in ("PRPOS", "PREND"):
                if tag in rec.info:
                    del rec.info[tag]
            outvcf.write(rec)
