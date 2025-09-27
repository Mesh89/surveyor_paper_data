import pysam
import argparse

def process_reads(bam_file, output_file, tolerance=50):
    # Open the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Get contig lengths from the BAM header
    contig_lengths = {contig['SN']: contig['LN'] for contig in bam.header['SQ']}

    # Open output BAM file, using the same header as the input BAM
    with pysam.AlignmentFile(output_file, "wb", header=bam.header) as out_bam:
        for read in bam:
            contig_len = contig_lengths.get(read.reference_name, None)
            if contig_len is None:
                continue
            
            # Get read start and end positions
            start = read.reference_start + 1  # 0-based to 1-based
            end = read.reference_end  # already 1-based in pysam

            # Read spans the entire contig within tolerance
            if (1 <= start <= tolerance + 1 and contig_len - tolerance <= end <= contig_len):
                out_bam.write(read)
                continue

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract reads based on complex criteria including tolerance, clipping, and mapped length.")
    parser.add_argument("bam_file", help="Input BAM file")
    parser.add_argument("output_file", help="Output BAM file")
    parser.add_argument("--tolerance", type=int, default=50, help="Tolerance in bp for the start and end positions of the reads (default: 50 bp)")

    args = parser.parse_args()
    process_reads(args.bam_file, args.output_file, args.tolerance)
