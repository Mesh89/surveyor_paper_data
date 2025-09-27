import pysam
import sys
import gzip

def extract_unique_reads(bam_files, output_fq_gz):
    seen = set()
    # Open the output in text mode with gzip
    with gzip.open(output_fq_gz, "wt") as fq_out:
        for bam_file in bam_files:
            with pysam.AlignmentFile(bam_file, "rb") as bam:
                for read in bam:
                    # Optionally skip unmapped, secondary, or supplementary reads
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue

                    if read.qname in seen:
                        continue
                    seen.add(read.qname)

                    # Convert base quality scores to ASCII by adding 33
                    qual_str = "".join(chr(q + 33) for q in read.query_qualities)
                    fq_out.write(f"@{read.qname}\n{read.query_sequence}\n+\n{qual_str}\n")

if __name__ == "__main__":
    # Usage: python extract_unique_reads.py bam1 bam2 output.fq.gz
    if len(sys.argv) < 4:
        print("Usage: python extract_unique_reads.py <bam1> <bam2> <output.fq.gz>")
        sys.exit(1)
    
    bam_files = sys.argv[1:-1]
    output_fq_gz = sys.argv[-1]
    extract_unique_reads(bam_files, output_fq_gz)
