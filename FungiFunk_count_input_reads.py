#!/usr/bin/env python3

# Lets us define and read command-line arguments
import argparse

# Opens gzip-compressed FASTQ files
import gzip

# Provides access to system-specific functions, such as exiting the script
import sys

# Gives us a cleaner and safer way to handle file paths
from pathlib import Path


# Parse command-line arguments
def parse_args():
    """
    Parse command-line arguments.
    """

    # Create the argument parser
    parser = argparse.ArgumentParser(
        description="Count input reads in paired-end FASTQ files for FungiFunk."
    )

    # R1 input FASTQ file
    parser.add_argument(
        "--reads1",
        required=True,
        help="Path to paired-end reads 1 FASTQ or FASTQ.GZ file.",
    )

    # R2 input FASTQ file
    parser.add_argument(
        "--reads2",
        required=True,
        help="Path to paired-end reads 2 FASTQ or FASTQ.GZ file.",
    )

    # Sample prefix used by FungiFunk
    parser.add_argument(
        "--prefix",
        required=True,
        help="Sample prefix used for the output file name.",
    )

    # Output directory
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory where the input read count table will be written.",
    )

    # Return parsed arguments
    return parser.parse_args()


# Print warning messages to stderr so they are separated from regular output
def warn_msg(message):
    """
    Print a warning message to stderr.
    """

    print(f"WARNING: {message}", file=sys.stderr, flush=True)


# Print an error message to stderr and stop the script with a non-zero exit code
def error_exit(message):
    """
    Print an error message and stop the program.
    """

    print(f"ERROR: {message}", file=sys.stderr, flush=True)
    sys.exit(1)


# Open plain-text or gzip-compressed text files
def open_text_maybe_gzip(path):
    """
    Open a plain-text or gzip-compressed file in text mode.
    """

    # Use gzip.open for files ending in .gz
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")

    # Use regular open for uncompressed files
    return open(path, "rt", encoding="utf-8", errors="replace")


# Count reads in one FASTQ file
def count_fastq_reads(path):
    """
    Count FASTQ records by counting lines and dividing by four.
    """

    # Convert path to a pathlib object
    path = Path(path)

    # Stop if the input FASTQ file does not exist
    if not path.exists():
        error_exit(f"Input FASTQ file does not exist: {path}")

    # Store the number of lines found in the FASTQ file
    line_count = 0

    # Count lines without loading the full FASTQ file into memory
    try:
        with open_text_maybe_gzip(path) as handle:
            for line_count, _ in enumerate(handle, 1):
                pass

    # Stop if the FASTQ file cannot be read
    except Exception as exc:
        error_exit(f"Could not read FASTQ file {path}: {exc}")

    # Warn if the FASTQ line count does not look valid
    if line_count % 4 != 0:
        warn_msg(
            f"FASTQ line count is not divisible by 4 for {path}. "
            "The reported read count uses floor(line_count / 4)."
        )

    # FASTQ files have four lines per read record
    return line_count // 4


# Write input read counts for one FungiFunk sample
def write_input_read_counts(reads1, reads2, prefix, outdir):
    """
    Count R1 and R2 input reads and write a small TSV summary.
    """

    # Convert the output directory to a pathlib object
    outdir = Path(outdir)

    # Create the output directory if it does not exist
    outdir.mkdir(parents=True, exist_ok=True)

    # Count reads in R1 and R2
    R1_input_reads = count_fastq_reads(reads1)
    R2_input_reads = count_fastq_reads(reads2)

    # Calculate the total number of input reads used by FungiFunk
    total_input_reads = R1_input_reads + R2_input_reads

    # Check whether the paired-end files contain the same number of reads
    if R1_input_reads == R2_input_reads:
        paired_read_count_status = "matched"
    else:
        paired_read_count_status = "mismatched"

    # Define the output table path
    output_path = outdir / f"{prefix}_input_read_counts.tsv"

    # Write the input read count table
    with open(output_path, "w", encoding="utf-8") as handle:
        handle.write(
            "sample\tR1_input_reads\tR2_input_reads\ttotal_input_reads\tpaired_read_count_status\n"
        )
        handle.write(
            f"{prefix}\t{R1_input_reads}\t{R2_input_reads}\t{total_input_reads}\t{paired_read_count_status}\n"
        )

    # Return the output path so the wrapper can report it
    return output_path


# Main program
def main():
    """
    Run the input-read counting workflow.
    """

    # Read command-line arguments
    args = parse_args()

    # Count input reads and write the output table
    write_input_read_counts(
        reads1=args.reads1,
        reads2=args.reads2,
        prefix=args.prefix,
        outdir=args.outdir,
    )



# Run the program
if __name__ == "__main__":
    main()