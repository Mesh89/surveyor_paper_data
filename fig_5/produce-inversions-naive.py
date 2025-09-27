#!/usr/bin/env python3
import argparse
import pysam
import pyfaidx

def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(complement)[::-1]

def main():
    parser = argparse.ArgumentParser(
        description="For each inversion in the VCF, output a FASTA contig with 100kbp upstream + inverted sequence + 100kbp downstream, "
                    "as well as a FASTA contig with the reference contig and the 100kbp flanking regions."
    )
    parser.add_argument("vcf", help="Input VCF file with inversions")
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument("output_alt", help="Output FASTA file")
    parser.add_argument("output_ref", help="Output FASTA file")
    parser.add_argument("--flank", type=int, default=100000, help="Flanking size (default: 100000 bp)")
    args = parser.parse_args()

    # Open FASTA using pyfaidx (1-based indexing)
    fasta = pyfaidx.Fasta(args.fasta)

    # Open VCF using pysam.VariantFile
    vcf = pysam.VariantFile(args.vcf)

    with open(args.output_alt, "w") as out_alt, open(args.output_ref, "w") as out_ref:
        for record in vcf:
            chrom = record.chrom
            start = record.pos          # VCF positions are 1-based
            end = record.stop

            upstream_start = max(1, start - args.flank)
            upstream_end = start

            contig_length = len(fasta[chrom])
            downstream_start = end
            downstream_end = min(contig_length, end + args.flank)

            upstream_seq = fasta[chrom][upstream_start:upstream_end].seq
            ref_inversion_seq = fasta[chrom][start:end].seq
            if 'SVINSSEQ' in record.info:
                alt_inversion_seq = record.info['SVINSSEQ']
            else:
                inv_start, inv_end = start, end
                if 'INVPOS' in record.info:
                    inv_start, inv_end = record.info['INVPOS']
                    inv_start -= 1  # Adjust to 0-based for pysam
                alt_inversion_seq = reverse_complement(fasta[chrom][inv_start:inv_end].seq)
            downstream_seq = fasta[chrom][downstream_start:downstream_end].seq

            contig_id = f"{record.id}_{len(upstream_seq)}_{len(downstream_seq)}"

            # Write ALT
            alt_seq = upstream_seq + alt_inversion_seq + downstream_seq
            out_alt.write(f">{contig_id}\n")
            out_alt.write(alt_seq + "\n")

            # Write REF
            ref_seq = upstream_seq + ref_inversion_seq + downstream_seq
            out_ref.write(f">{contig_id}\n")
            out_ref.write(ref_seq + "\n")

if __name__ == "__main__":
    main()
