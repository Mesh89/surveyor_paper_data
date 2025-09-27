import pysam, ssw
import argparse

aligner = ssw.Aligner(matrix=ssw.NucleotideScoreMatrix(match=1, mismatch=-4), gap_open=6, gap_extend=1)

def reverse_complement(seq):
    complement_table = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return seq.translate(complement_table)[::-1]

def process_inversions(vcf_file, fasta_file, output_ref_file, output_alt_file, init_flank_size=2000, starting_inv_seq_len=2000):
    # Open the VCF file using pysam
    vcf = pysam.VariantFile(vcf_file)
    
    # Open the reference fasta file using pysam
    fasta = pysam.FastaFile(fasta_file)
    
    # Open the output FASTA file
    with open(output_ref_file, 'w') as output_ref_fasta, open(output_alt_file, 'w') as output_alt_fasta:
        for record in vcf:
            chrom = record.chrom
            start = record.pos - 1  # VCF is 1-based, adjust to 0-based for pysam
            end = record.stop  # This should be provided in the VCF for inversions
            
            # Fetch upstream, inverted sequence, and downstream
            flank_size = init_flank_size

            while True:
                upstream = fasta.fetch(chrom, max(0, start-flank_size), start)
                downstream = fasta.fetch(chrom, end, end + flank_size)

                # align upstream and the rc of downstream
                alignment = aligner.align(upstream, reverse_complement(downstream))
                aln_size = max(alignment.reference_end-alignment.reference_begin, alignment.query_end-alignment.query_begin)
                if start-flank_size < 0 or end+flank_size > fasta.get_reference_length(chrom) or aln_size < flank_size-1000:
                    break
                else:
                    flank_size += 1000

            inv_seq_len = starting_inv_seq_len
            ref_inverted_sequence = fasta.fetch(chrom, start, end)
            if 'SVINSSEQ' in record.info:
                alt_inverted_sequence = record.info['SVINSSEQ']
            else:
                inv_start, inv_end = start, end
                if 'INVPOS' in record.info:
                    inv_start, inv_end = record.info['INVPOS']
                    inv_start -= 1  # Adjust to 0-based for pysam
                alt_inverted_sequence = reverse_complement(fasta.fetch(chrom, inv_start, inv_end))
                while len(ref_inverted_sequence) > 2*inv_seq_len:
                    inv_prefix, inv_suffix_rc = ref_inverted_sequence[:inv_seq_len], reverse_complement(ref_inverted_sequence[-inv_seq_len:])
                    alignment = aligner.align(inv_prefix, inv_suffix_rc)
                    aln_size = max(alignment.reference_end-alignment.reference_begin, alignment.query_end-alignment.query_begin)
                    if aln_size < inv_seq_len-1000:
                        break
                    inv_seq_len += 2000

            if len(ref_inverted_sequence) > 2*inv_seq_len:
                output_ref_fasta.write(f">{record.id}_REF_{flank_size}_LB\n")
                output_ref_fasta.write(f"{upstream + ref_inverted_sequence[:inv_seq_len]}\n")

                output_ref_fasta.write(f">{record.id}_REF_{flank_size}_RB\n")
                output_ref_fasta.write(f"{ref_inverted_sequence[-inv_seq_len:] + downstream}\n")
            else:
                output_ref_fasta.write(f">{record.id}_REF_{flank_size}\n")
                output_ref_fasta.write(f"{upstream + ref_inverted_sequence + downstream}\n")

            if len(alt_inverted_sequence) > 2*inv_seq_len:
                output_alt_fasta.write(f">{record.id}_ALT_{flank_size}_LB\n")
                output_alt_fasta.write(f"{upstream + alt_inverted_sequence[:inv_seq_len]}\n")
                
                output_alt_fasta.write(f">{record.id}_ALT_{flank_size}_RB\n")
                output_alt_fasta.write(f"{alt_inverted_sequence[-inv_seq_len:] + downstream}\n")
            else:
                output_alt_fasta.write(f">{record.id}_ALT_{flank_size}\n")
                output_alt_fasta.write(f"{upstream + alt_inverted_sequence + downstream}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF inversions and extract sequences to a FASTA file.")
    parser.add_argument("vcf_file", help="Path to the VCF file describing inversions.")
    parser.add_argument("fasta_file", help="Path to the reference FASTA file.")
    parser.add_argument("output_ref_file", help="Path to the output FASTA file with REF alleles.")
    parser.add_argument("output_alt_file", help="Path to the output FASTA file with ALT alleles.")
    parser.add_argument("--flank_size", type=int, default=2000, help="Starting size of the flanking regions to extract (default: 2000 bp).")

    args = parser.parse_args()

    process_inversions(args.vcf_file, args.fasta_file, args.output_ref_file, args.output_alt_file, args.flank_size)
