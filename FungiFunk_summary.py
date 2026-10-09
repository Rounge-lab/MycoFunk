#!/usr/bin/env python3

"""
FungiFunk_summary.py

Create comparative summary tables and an HTML report from multiple FungiFunk runs.

This module reads FungiFunk annotation/count outputs from several samples and creates:
- raw count matrices
- optional CPM-normalised matrices
- protein metadata
- annotation summaries
- detected feature summaries
- an HTML report by default
"""

# Lets us define and read command-line arguments
import argparse

# Provides regular expressions for matching and extracting filename patterns
import re

# Provides access to system-specific functions, such as exiting the script
import sys

# Provides dictionaries with default values, useful for grouping files by sample or table type
from collections import defaultdict

# Import datetime to keep track of the time for the timestamps
from datetime import datetime

# Gives us a cleaner and safer way to handle file paths
from pathlib import Path

# Imports pandas, which we use to read, merge, and save tables
import pandas as pd

# Escape HTML special characters before inserting values into the HTML report
from html import escape

# Build Plotly figures with manual traces
import plotly.graph_objects as go


# Define a function that reads and returns the command-line arguments
def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Summarise and compare outputs from multiple FungiFunk runs."
    )

    # Add the required argument --input_dir
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing FungiFunk outputs. The directory is searched recursively.",
    )

    # Add the required argument --outdir
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for FungiFunk summary results.",
    )

    # Add the optional argument --threads
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help=(
            "Number of worker processes used by FungiFunk_summary. "
            "Mainly used when summarising annotation tables. "
            "Higher values may reduce runtime but increase RAM usage. "
            "Default: 1."
        ),
    )

    # Add the optional argument --normalise_cpm
    parser.add_argument(
        "--normalise_cpm",
        action="store_true",
        help=(
            "Write additional CPM-normalised matrices. CPM is calculated within each "
            "count table as count / total sample count * 1,000,000."
        ),
    )

    # Add the optional argument --store_true
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Require all expected FungiFunk output tables for every detected sample. "
            "Without this option, missing files are reported as warnings and the summary continues."
        ),
    )

    # Add the optional argument --skip_html_report
    parser.add_argument(
        "--skip_html_report",
        action="store_true",
        help="Do not generate the HTML report. By default, the HTML report is generated.",
    )

    # Add the optional argument --html_top_n
    parser.add_argument(
        "--html_top_n",
        type=int,
        default=50,
        help=(
            "Number of top functional features to display in the HTML report for each "
            "annotation type. Features are ranked by total count across samples. "
            "This only affects the HTML report; complete TSV matrices are always written. "
            "Default: 50."
        ),
    )

    # Add the optional argument --html_no_top
    parser.add_argument(
        "--html_no_top",
        action="store_true",
        help=(
            "Do not limit the number of functional features shown in the HTML report. "
            "This may create a very large HTML file and may make the browser slow."
        ),
    )

    parser.add_argument(
        "--calculate_pfam_eggnog_context",
        action="store_true",
        help=(
            "Calculate dataset-derived context for PFAM and eggNOG OG features. "
            "This adds top-associated protein name and description columns to PFAM "
            "and eggNOG OG metadata and HTML tables. "
            "This step scans annotation tables and may increase runtime."
        ),
    )

    # Add the optional argument --timestamps
    parser.add_argument(
        "--timestamps",
        action="store_true",
        help="Add timestamps to FungiFunk_summary progress messages.",
    )

    # Run metadata log written by the FungiFunk_summary wrapper
    parser.add_argument(
        "--run_info_log",
        default=None,
        help="Path to the FungiFunk_summary run information log.",
    )

    # Parse the arguments provided by the user and return them
    return parser.parse_args()


# Values that should be treated as missing annotations
MISSING_VALUES = {
    "",
    "-",
    "--",
    "NA",
    "N/A",
    "na",
    "NaN",
    "nan",
    "None",
    "NONE",
    "null",
    "NULL",
}

# Count table definitions used by FungiFunk
COUNT_TABLES = {
    "protein": {
        "suffix": "_protein_read_counts.tsv",
        "id_col": "protein_id",
        "raw_filename": "protein_read_count_matrix.tsv",
        "cpm_filename": "protein_read_count_cpm_matrix.tsv",
        "title": "Protein features",
    },
    "ko": {
        "suffix": "_ko_read_counts.tsv",
        "id_col": "KEGG_ko",
        "raw_filename": "ko_read_count_matrix.tsv",
        "cpm_filename": "ko_read_count_cpm_matrix.tsv",
        "title": "KEGG orthologs",
    },
    "eggnog_ogs": {
        "suffix": "_eggnog_ogs_read_counts.tsv",
        "id_col": "eggNOG_OGs",
        "raw_filename": "eggnog_ogs_read_count_matrix.tsv",
        "cpm_filename": "eggnog_ogs_read_count_cpm_matrix.tsv",
        "title": "eggNOG orthologous groups",
    },
    "cazy": {
        "suffix": "_cazy_read_counts.tsv",
        "id_col": "CAZy",
        "raw_filename": "cazy_read_count_matrix.tsv",
        "cpm_filename": "cazy_read_count_cpm_matrix.tsv",
        "title": "CAZy families",
    },
    "pfam": {
        "suffix": "_pfam_read_counts.tsv",
        "id_col": "PFAMs",
        "raw_filename": "pfam_read_count_matrix.tsv",
        "cpm_filename": "pfam_read_count_cpm_matrix.tsv",
        "title": "PFAM features",
    },
}

# Additional FungiFunk output tables
EXTRA_TABLES = {
    "input_read_counts": {
        "suffix": "_input_read_counts.tsv",
    },
    "annotation": {
        "suffix": "_annotation.tsv",
    },
    "unannotated": {
        "suffix": "_unannotated_hits.tsv",
    },
}

# Columns used to create protein metadata
PROTEIN_METADATA_COLUMNS = [
    "protein_id",
    "Description",
    "Preferred_name",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
]

# Describe broad CAZy classes using their standard family prefixes
CAZY_CLASS_DESCRIPTIONS = {
    "GH": "Glycoside Hydrolase",
    "GT": "Glycosyltransferase",
    "PL": "Polysaccharide Lyase",
    "CE": "Carbohydrate Esterase",
    "AA": "Auxiliary Activity",
    "CBM": "Carbohydrate-Binding Module",
}

# Limit how many associated values are written into one metadata cell
MAX_ASSOCIATED_VALUES = 10


# Print progress messages during the summary workflow
def log_msg(message, timestamps=False):
    """
    Print a progress message.

    If timestamps is True, the current time is added before the message.
    """
    # If timestamps is enabled, each message is printed with the current time
    if timestamps:
        print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)
    else:
        print(message, flush=True)


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


# Read the column names from a TSV file without loading the full table into memory
def read_header(path):
    """
    Read only the header line of a TSV file.
    """
    try:
        # Read only the header to avoid loading large TSV files into memory
        return list(pd.read_csv(path, sep="\t", nrows=0).columns)

    except Exception as exc:
        # Add a clearer error message while preserving the original error
        raise RuntimeError(f"Could not read header from {path}: {exc}") from exc


# Check whether a value should be treated as missing in FungiFunk tables
def is_missing_value(value):
    """
    Return True if a value should be treated as missing.
    """
    # Treat pandas missing values, such as NaN, as missing
    if pd.isna(value):
        return True

    # Convert the value to text and remove leading/trailing spaces
    value = str(value).strip()

    # Check whether the cleaned value is one of the predefined missing-value labels
    return value in MISSING_VALUES


# Split annotation fields with several comma-separated values into a clean list
def split_annotation_values(value):
    """
    Split a comma-separated annotation field into clean individual values.
    """
    # If the full field is missing, return an empty list
    if is_missing_value(value):
        return []

    # Create an empty list to store the cleaned annotation values
    values = []

    # Split the field by commas because some annotations contain multiple values
    for item in str(value).split(","):
        # Remove extra spaces around each annotation value
        item = item.strip()

        # Keep only real annotation values and ignore missing-value labels
        if not is_missing_value(item):
            values.append(item)

    # Return the list of cleaned annotation values
    return values


# Extract all unique annotation values from a pandas column
def extract_unique_values(series):
    """
    Extract unique annotation values from a pandas Series.
    """
    # Use a set so repeated annotation values are stored only once
    unique_values = set()

    # Drop pandas missing values before processing the column
    for value in series.dropna():
        # Split comma-separated annotation fields and add all values to the set
        unique_values.update(split_annotation_values(value))

    # Return the complete set of unique annotation values
    return unique_values


# Count the number of data rows in a TSV file without loading it with pandas
def count_data_lines(path):
    """
    Count data lines in a TSV file, excluding the header.
    """
    # Open the file safely as text and replace problematic characters if needed
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        # Count all lines in the file, including the header
        line_count = sum(1 for _ in handle)

    # Subtract the header line and avoid returning negative values for empty files
    return max(line_count - 1, 0)


# Find the first real value in a list-like object
def first_valid_value(values):
    """
    Return the first non-missing value from a list-like object.
    """
    # Check each value in order
    for value in values:
        # Return the first value that is not considered missing
        if not is_missing_value(value):
            return str(value)

    # If all values are missing, return the standard missing-value label
    return "-"


# Extract a sample name from a FungiFunk output filename
def extract_sample_from_filename(filename, suffix):
    """
    Extract the sample name by removing a known suffix from a filename.
    """
    # If the filename does not match the expected suffix, it is not the right file type
    if not filename.endswith(suffix):
        return None

    # Remove the known suffix and keep everything before it as the sample name
    return filename[: -len(suffix)]


# Find FungiFunk output tables in the input directory and group them by sample
def discover_fungifunk_outputs(input_dir, strict=False):
    """
    Discover FungiFunk output files inside an input directory.

    The function searches recursively and extracts sample names from filenames.
    """
    # Store detected files by sample and table type
    samples = defaultdict(dict)
    # Store one record per detected file for the output inventory table
    records = []
    # Store duplicated files so they can be reported later
    duplicate_records = []

    # Combine count table definitions and extra table definitions into one dictionary
    all_definitions = {}
    all_definitions.update(COUNT_TABLES)
    all_definitions.update(EXTRA_TABLES)

    # Search recursively for all TSV files inside the input directory
    for path in sorted(input_dir.rglob("*.tsv")):
        # Get only the filename, without the full directory path
        filename = path.name

        # Check whether the filename matches any known FungiFunk output suffix
        for table_key, info in all_definitions.items():
            # Get the expected suffix for this table type
            suffix = info["suffix"]
            # Extract the sample name if the filename ends with this suffix
            sample = extract_sample_from_filename(filename, suffix)
            # Skip this file if it does not match the current table type
            if sample is None:
                continue

            # Check whether this sample already has a file for this table type
            if table_key in samples[sample]:
                # Store the duplicated file information
                duplicate_records.append((sample, table_key, samples[sample][table_key], path))
                # In strict mode, duplicated files are treated as fatal errors
                if strict:
                    error_exit(
                        f"Duplicate files detected for sample '{sample}' and table '{table_key}': "
                        f"{samples[sample][table_key]} and {path}"
                    )
                # In non-strict mode, keep the first file and warn the user
                warn_msg(
                    f"Duplicate files detected for sample '{sample}' and table '{table_key}'. "
                    f"Keeping first file: {samples[sample][table_key]}"
                )
                # Skip the duplicated file
                continue

            # Store the detected file path for this sample and table type
            samples[sample][table_key] = path

            # Add this detected file to the inventory records
            records.append(
                {
                    "sample": sample,
                    "table_type": table_key,
                    "path": str(path),
                }
            )

    # Convert the inventory records into a pandas table
    inventory_df = pd.DataFrame(records)

    # Return the grouped files, the inventory table, and duplicated file records
    return dict(samples), inventory_df, duplicate_records


# Check that the discovered FungiFunk outputs are complete enough to summarise
def validate_discovered_outputs(samples, strict=False):
    """
    Validate discovered output files.

    In strict mode, every sample must have all expected tables.
    """
    # Build the list of all table types expected for each sample
    expected_tables = list(COUNT_TABLES.keys()) + list(EXTRA_TABLES.keys())

    # Stop immediately if no FungiFunk output files were found
    if not samples:
        error_exit("No FungiFunk output files were detected in --input_dir.")

    # Check the detected files for each sample
    for sample, files in samples.items():
        # Find which expected tables are missing for this sample
        missing = [table for table in expected_tables if table not in files]

        # In strict mode, missing tables are treated as fatal errors
        if missing and strict:
            error_exit(
                f"Sample '{sample}' is missing expected output table(s): {', '.join(missing)}"
            )

        # In non-strict mode, missing tables are reported as warnings
        if missing:
            warn_msg(
                f"Sample '{sample}' is missing output table(s): {', '.join(missing)}"
            )


# Create a natural sorting key for sample names containing numbers
def natural_sort_key(value):
    """
    Create a sorting key that orders text with numbers in a human-friendly way.
    """

    # Convert the value to text
    value = str(value)

    # Split text into numeric and non-numeric parts
    parts = re.split(r"(\d+)", value)

    # Convert numeric parts to integers and text parts to lowercase strings
    return [
        int(part) if part.isdigit() else part.lower()
        for part in parts
    ]


# Sort a table by sample name using natural alphanumeric order
def sort_table_by_sample(table):
    """
    Sort a table by sample name using natural alphanumeric order.
    """

    # Return the table unchanged if it is missing, empty, or has no sample column
    if table is None or table.empty or "sample" not in table.columns:
        return table

    # Sort rows using natural ordering of sample names
    sorted_row_index = sorted(
        table.index,
        key=lambda row_index: natural_sort_key(table.loc[row_index, "sample"]),
    )

    # Return the sorted table with a clean index
    return table.loc[sorted_row_index].reset_index(drop=True)


# Sort matrix columns using natural alphanumeric order
def sort_matrix_columns_natural(matrix, protected_columns=None):
    """
    Sort matrix columns using natural alphanumeric order.

    Protected columns, such as feature or taxon ID columns, are kept at the
    beginning and are not included in the sorting.
    """

    # Return the matrix unchanged if it is missing or empty
    if matrix is None or matrix.empty:
        return matrix

    # Use an empty list if no protected columns were provided
    if protected_columns is None:
        protected_columns = []

    # Keep only protected columns that are actually present in the matrix
    protected_columns = [
        column for column in protected_columns
        if column in matrix.columns
    ]

    # Define the columns that should be sorted
    sortable_columns = [
        column for column in matrix.columns
        if column not in protected_columns
    ]

    # Sort the remaining columns in human-friendly alphanumeric order
    sorted_columns = sorted(
        sortable_columns,
        key=natural_sort_key,
    )

    # Return the matrix with protected columns first and sorted sample columns after
    return matrix.loc[:, protected_columns + sorted_columns]


# Build a per-sample overview of which expected FungiFunk output files were found
def build_sample_summary(samples):
    """
    Build a summary table describing which files were detected for each sample.
    """
    # Store one summary row per sample
    rows = []

    # Build the list of all table types expected for each sample
    expected_tables = list(COUNT_TABLES.keys()) + list(EXTRA_TABLES.keys())

    # Process samples in natural alphabetical order for reproducible output
    for sample in sorted(samples, key=natural_sort_key):
        # Start the row with the sample name
        row = {"sample": sample}

        # Add presence/absence information and file paths for each expected table
        for table_key in expected_tables:
            # Get the detected file path for this table type, or None if it is missing
            path = samples[sample].get(table_key)

            # Record whether this table was found for the sample
            row[f"has_{table_key}"] = path is not None
            # Record the file path if available, otherwise leave the field empty
            row[f"{table_key}_path"] = str(path) if path is not None else ""

        # Add the completed sample row to the summary
        rows.append(row)

    # Convert all summary rows into a pandas table
    return pd.DataFrame(rows)


# Read a FungiFunk count table and return counts indexed by feature ID
def read_count_table(path, id_col):
    """
    Read one FungiFunk count table and return a clean Series of counts.
    """
    # Read only the header first to check that the expected columns exist
    header = read_header(path)

    # Define the columns that every count table must contain
    required_cols = [id_col, "read_count"]

    # Check whether any required columns are missing from the file
    missing_cols = [col for col in required_cols if col not in header]

    # Stop with a clear error if the table does not have the expected format
    if missing_cols:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(missing_cols)}"
        )

    # Read only the required columns to reduce memory usage
    df = pd.read_csv(
        path,
        sep="\t",
        usecols=required_cols,
        dtype={id_col: "string"},
    )

    # Clean feature IDs by converting them to text and removing extra spaces
    df[id_col] = df[id_col].astype("string").str.strip()

    # Remove rows where the feature ID is missing
    df = df[~df[id_col].isna()]
    # Remove rows where the feature ID is one of the predefined missing labels
    df = df[~df[id_col].isin(MISSING_VALUES)]

    # Convert read_count to numeric values and replace invalid counts with zero
    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce").fillna(0)

    # Sum counts by feature ID in case the same feature appears more than once
    counts = df.groupby(id_col)["read_count"].sum()

    # Return a pandas Series where the index is the feature ID and the value is the count
    return counts


# Build a count matrix for one FungiFunk annotation/count table type
def build_count_matrix(samples, table_key, strict=False):
    """
    Build one count matrix for a specific annotation type.
    """
    # Get the table definition, including the feature ID column name
    table_info = COUNT_TABLES[table_key]
    id_col = table_info["id_col"]

    # Store one count Series per sample
    series_list = []

    # Process samples in natural alphabetical order for reproducible output
    for sample in sorted(samples, key=natural_sort_key):
        # Get the file path for this table type in the current sample
        path = samples[sample].get(table_key)

        # Skip samples where this table type was not detected
        if path is None:
            continue

        # Try to read and clean the count table
        try:
            counts = read_count_table(path, id_col)
        # Handle unreadable or malformed count tables
        except Exception as exc:
            # In strict mode, stop immediately if a table cannot be read
            if strict:
                error_exit(f"Could not read {table_key} count table for sample '{sample}': {exc}")

            # In non-strict mode, warn the user and skip this table
            warn_msg(f"Skipping {table_key} count table for sample '{sample}': {exc}")
            continue

        # Name the Series with the sample name so it becomes the matrix column name
        counts.name = sample
        # Add this sample's counts to the list of matrix columns
        series_list.append(counts)

    # If no valid count tables were found, return an empty table
    if not series_list:
        return pd.DataFrame()

    # Combine all sample count Series into one matrix and fill missing values with zero
    matrix = pd.concat(series_list, axis=1).fillna(0)

    # Sort feature IDs alphabetically for reproducible output
    matrix = matrix.sort_index()

    # Set the matrix index name to the feature ID column name
    matrix.index.name = id_col

    # Make sure all matrix values are numeric
    for col in matrix.columns:
        matrix[col] = pd.to_numeric(matrix[col], errors="coerce").fillna(0)

    # If all values are whole numbers, save them as integers instead of floats
    if (matrix % 1 == 0).all().all():
        matrix = matrix.astype("int64")

    # Return the final count matrix
    return matrix


# Write all raw count matrices and, if requested, their CPM-normalised versions
def write_count_matrices(samples, matrices_dir, normalise_cpm=False, strict=False, time_stamps=False):
    """
    Write raw count matrices and optional CPM matrices.
    """
    # Store the raw matrices in memory so they can also be used later for the HTML report
    raw_matrices = {}
    # Store the CPM matrices if CPM normalisation is requested
    cpm_matrices = {}

    # Build one matrix for each FungiFunk count table type
    for table_key, table_info in COUNT_TABLES.items():
        # Build the raw count matrix for this table type
        matrix = build_count_matrix(samples, table_key, strict=strict)

        # Save the raw matrix in the dictionary using the table type as key
        raw_matrices[table_key] = matrix

        # Define the output path for the raw matrix
        raw_path = matrices_dir / table_info["raw_filename"]

        # Write the raw matrix to a TSV file
        matrix.to_csv(raw_path, sep="\t")

        # Tell the user where the raw matrix was writen
        log_msg(f"    (Written to: {raw_path})", time_stamps)

        # If requested, also calculate and write a CPM-normalised matrix
        if normalise_cpm:
            # Convert raw counts to CPM values within this count table
            cpm_matrix = calculate_cpm_matrix(matrix)

            # Save the CPM matrix in the dictionary using the table type as key
            cpm_matrices[table_key] = cpm_matrix

            # Define the output path for the CPM matrix
            cpm_path = matrices_dir / table_info["cpm_filename"]

            # Write the CPM matrix to a TSV file
            cpm_matrix.to_csv(cpm_path, sep="\t")

            # Tell the user where the CPM-normalised matrix was writen
            log_msg(f"    (Written to: {cpm_path})", time_stamps)

    # Return the matrices so they can be reused by other parts of the workflow
    return raw_matrices, cpm_matrices


# Convert raw count matrices to CPM-normalised matrices
def calculate_cpm_matrix(matrix):
    """
    Convert a raw count matrix to CPM values per sample column.
    """
    # If the matrix is empty, return an empty copy without further processing
    if matrix.empty:
        return matrix.copy()

    # Convert counts to floats because CPM values can contain decimals
    cpm_matrix = matrix.astype("float64").copy()

    # Calculate CPM independently for each sample column
    for sample in cpm_matrix.columns:
        # Calculate the total number of counts for this sample in this table
        total_count = cpm_matrix[sample].sum()

        # If the sample has zero total counts, keep the whole column as zero
        if total_count == 0:
            cpm_matrix[sample] = 0.0
        # Otherwise, scale each feature count to counts per million table counts
        else:
            cpm_matrix[sample] = (cpm_matrix[sample] / total_count) * 1_000_000

    # Round CPM values to keep the output readable
    cpm_matrix = cpm_matrix.round(6)

    # Return the CPM-normalised matrix
    return cpm_matrix


# Build one non-redundant protein metadata table from all protein count tables
def build_protein_metadata(samples, strict=False):
    """
    Build a protein metadata table from all protein read count tables.
    """
    # Store protein metadata tables read from each sample
    metadata_frames = []

    # Process samples in natural alphabetical order for reproducible output
    for sample in sorted(samples, key=natural_sort_key):
        # Get the protein read count table for this sample
        path = samples[sample].get("protein")

        # Skip samples without a protein count table
        if path is None:
            continue

        # Try to read only the header first to check which metadata columns are available
        try:
            header = read_header(path)
        # Handle unreadable protein table headers
        except Exception as exc:
            # In strict mode, stop immediately if the header cannot be read
            if strict:
                error_exit(f"Could not read protein table header for sample '{sample}': {exc}")
            # In non-strict mode, warn the user and skip this sample for protein metadata
            warn_msg(f"Skipping protein metadata for sample '{sample}': {exc}")
            continue

        # Keep only the expected metadata columns that are actually present in this file
        available_cols = [col for col in PROTEIN_METADATA_COLUMNS if col in header]

        # protein_id is required because it is used to merge metadata across samples
        if "protein_id" not in available_cols:
            # In strict mode, stop if protein_id is missing
            if strict:
                error_exit(f"Protein table for sample '{sample}' does not contain protein_id.")
            # In non-strict mode, warn the user and skip this table
            warn_msg(f"Protein table for sample '{sample}' does not contain protein_id.")
            continue

        # Try to read only the available metadata columns to reduce memory usage
        try:
            df = pd.read_csv(
                path,
                sep="\t",
                usecols=available_cols,
                dtype="string",
            )
        # Handle malformed or unreadable protein metadata tables
        except Exception as exc:
            # In strict mode, stop immediately if the table cannot be read
            if strict:
                error_exit(f"Could not read protein metadata for sample '{sample}': {exc}")
            # In non-strict mode, warn the user and skip this sample
            warn_msg(f"Skipping protein metadata for sample '{sample}': {exc}")
            continue

        # Add missing metadata columns as "-" so all samples have the same column structure
        for col in PROTEIN_METADATA_COLUMNS:
            if col not in df.columns:
                df[col] = "-"

        # Reorder columns to the standard protein metadata order
        df = df[PROTEIN_METADATA_COLUMNS]

        # Store this sample's protein metadata table
        metadata_frames.append(df)

    # If no valid protein metadata tables were found, return an empty table with standard columns
    if not metadata_frames:
        return pd.DataFrame(columns=PROTEIN_METADATA_COLUMNS)

    # Combine protein metadata from all samples into one table
    combined = pd.concat(metadata_frames, ignore_index=True)

    # Clean protein IDs by converting them to text and removing extra spaces
    combined["protein_id"] = combined["protein_id"].astype("string").str.strip()

    # Remove rows where protein_id is missing
    combined = combined[~combined["protein_id"].isna()]
    # Remove rows where protein_id is one of the predefined missing-value labels
    combined = combined[~combined["protein_id"].isin(MISSING_VALUES)]

    # Collapse duplicated protein IDs and keep the first valid metadata value for each column
    metadata = (
        combined.groupby("protein_id", as_index=False)
        .agg({col: first_valid_value for col in PROTEIN_METADATA_COLUMNS if col != "protein_id"})
    )

    # Reorder columns to the standard protein metadata order
    metadata = metadata[PROTEIN_METADATA_COLUMNS]

    # Return the final non-redundant protein metadata table
    return metadata


# Convert one input read count value into an integer or a missing value
def parse_input_read_count_value(value):
    """
    Convert one input read count value into an integer.

    Missing or invalid values are returned as pd.NA.
    """

    # Return a missing value if the input should be treated as missing
    if is_missing_value(value):
        return pd.NA

    # Convert the value to text and remove extra spaces
    value = str(value).strip()

    # Try to convert the value to an integer
    try:
        return int(float(value))

    # Return missing if the value cannot be converted
    except ValueError:
        return pd.NA


# Read input read counts for one sample
def read_input_read_counts_for_sample(samples, sample, strict=False):
    """
    Read the input read count table for one sample.
    """

    # Define default values used when the input read count table is missing
    default_values = {
        "R1_input_reads": pd.NA,
        "R2_input_reads": pd.NA,
        "total_input_reads": pd.NA,
    }

    # Get the input read count table path for this sample
    path = samples[sample].get("input_read_counts")

    # Return missing values if the input read count table was not detected
    if path is None:
        return default_values

    # Try to read the input read count table
    try:
        table = pd.read_csv(path, sep="\t", dtype=str, na_filter=False)

    # Handle files that cannot be read
    except Exception as exc:

        # Stop immediately in strict mode
        if strict:
            error_exit(
                f"Could not read input read count table for sample '{sample}': {exc}"
            )

        # Warn and return missing values in non-strict mode
        warn_msg(
            f"Could not read input read count table for sample '{sample}': {exc}"
        )
        return default_values

    # Check whether the table contains at least one row
    if table.empty:
        if strict:
            error_exit(f"Input read count table is empty for sample '{sample}': {path}")

        warn_msg(f"Input read count table is empty for sample '{sample}': {path}")
        return default_values

    # Define the columns expected in the input read count table
    required_columns = [
        "R1_input_reads",
        "R2_input_reads",
        "total_input_reads",
    ]

    # Find required columns that are missing from the table
    missing_columns = [
        column for column in required_columns if column not in table.columns
    ]

    # Stop or warn if the table does not contain the expected columns
    if missing_columns:
        message = (
            f"Input read count table for sample '{sample}' is missing column(s): "
            f"{', '.join(missing_columns)}"
        )

        if strict:
            error_exit(message)

        warn_msg(message)
        return default_values

    # Use the first row because each input read count table should contain one sample
    row = table.iloc[0]

    # Return the parsed input read count values
    return {
        "R1_input_reads": parse_input_read_count_value(row["R1_input_reads"]),
        "R2_input_reads": parse_input_read_count_value(row["R2_input_reads"]),
        "total_input_reads": parse_input_read_count_value(row["total_input_reads"]),
    }


# Divide two count values while preserving missing values
def divide_count_values(numerator, denominator, digits=6):
    """
    Divide two count values and return a rounded fraction.

    If the numerator or denominator is missing, or if the denominator is zero,
    return pd.NA.
    """

    # Return missing if one of the values is not available
    if pd.isna(numerator) or pd.isna(denominator):
        return pd.NA

    # Try to convert both values to numbers before division
    try:
        numerator = float(numerator)
        denominator = float(denominator)

    # Return missing if conversion fails
    except (TypeError, ValueError):
        return pd.NA

    # Avoid division by zero
    if denominator == 0:
        return pd.NA

    # Return the rounded fraction
    return round(numerator / denominator, digits)


# Add two count values while preserving missing values
def add_count_values(first_count, second_count):
    """
    Add two count values.

    If either count is missing, return pd.NA.
    """

    # Return missing if one of the values is not available
    if pd.isna(first_count) or pd.isna(second_count):
        return pd.NA

    # Return the integer sum
    return int(first_count) + int(second_count)


# Summarise one unannotated hits table without loading the full file into memory
def summarise_unannotated_hits_file(path, chunksize=100_000):
    """
    Summarise one unannotated hits table.

    The table is processed in chunks to reduce memory usage.
    """

    # Read only the header to check whether the mate column is available
    header = read_header(path)

    # If the mate column is missing, only the total number of unannotated hits can be counted
    if "mate" not in header:
        return {
            "unannotated_hit_count": count_data_lines(path),
            "R1_unannotated_hit_count": pd.NA,
            "R2_unannotated_hit_count": pd.NA,
        }

    # Store running counters that can be updated chunk by chunk
    counters = {
        "unannotated_hit_count": 0,
        "R1_unannotated_hit_count": 0,
        "R2_unannotated_hit_count": 0,
    }

    # Read only the mate column in chunks to reduce peak RAM usage
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=["mate"],
        chunksize=chunksize,
        dtype=str,
        na_filter=False,
    ):

        # Count how many unannotated hits are present in this chunk
        counters["unannotated_hit_count"] += len(chunk)

        # Count how many unannotated hits came from R1 and R2
        mate_counts = chunk["mate"].value_counts()
        counters["R1_unannotated_hit_count"] += int(mate_counts.get("R1", 0))
        counters["R2_unannotated_hit_count"] += int(mate_counts.get("R2", 0))

    # Return all unannotated hit counts
    return counters


# Summarise unannotated hits for one sample
def summarise_unannotated_hits_for_sample(samples, sample, strict=False):
    """
    Summarise unannotated hits for one sample.
    """

    # Define default values used when the unannotated hits table is missing
    default_values = {
        "unannotated_hit_count": pd.NA,
        "R1_unannotated_hit_count": pd.NA,
        "R2_unannotated_hit_count": pd.NA,
    }

    # Get the unannotated hits table for this sample
    path = samples[sample].get("unannotated")

    # Return missing values if the unannotated hits table is not available
    if path is None:
        return default_values

    # Try to summarise the unannotated hits table
    try:
        return summarise_unannotated_hits_file(path)

    # Handle unannotated hit tables that cannot be read
    except Exception as exc:

        # In strict mode, stop immediately if summarising fails
        if strict:
            error_exit(
                f"Could not summarise unannotated hits for sample '{sample}': {exc}"
            )

        # In non-strict mode, warn the user and return missing values
        warn_msg(f"Could not summarise unannotated hits for sample '{sample}': {exc}")
        return default_values


# Summarise one FungiFunk annotation table without loading the full file into memory
def summarise_annotation_file(path, chunksize=100_000):
    """
    Summarise one default FungiFunk annotation table.

    The table is processed in chunks to reduce memory usage.
    """

    # Define the annotation columns needed for the lightweight summary statistics
    wanted_cols = [
        "diamond_evalue",
        "bitscore",
        "mate",
    ]

    # Read only the header to check which expected columns are present
    header = read_header(path)

    # Keep only wanted columns that actually exist in this file
    usecols = [col for col in wanted_cols if col in header]

    # Store running counters that can be updated chunk by chunk
    counters = {
        "annotated_hit_count": 0,
        "bitscore_sum": 0.0,
        "bitscore_count": 0,
        "evalue_sum": 0.0,
        "evalue_count": 0,
        "R1_hit_count": 0,
        "R2_hit_count": 0,
    }

    # If none of the expected columns are present, count rows only and return missing values for statistics
    if not usecols:
        return {
            "annotated_hit_count": count_data_lines(path),
            "R1_hit_count": pd.NA,
            "R2_hit_count": pd.NA,
            "mean_bitscore": pd.NA,
            "mean_diamond_evalue": pd.NA,
        }

    # Read the annotation table in chunks to reduce peak RAM usage
    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        chunksize=chunksize,
        dtype=str,
        na_filter=False,
    ):

        # Count how many annotated hits are present in this chunk
        counters["annotated_hit_count"] += len(chunk)

        # Count how many hits came from R1 and R2
        if "mate" in chunk.columns:
            mate_counts = chunk["mate"].value_counts()
            counters["R1_hit_count"] += int(mate_counts.get("R1", 0))
            counters["R2_hit_count"] += int(mate_counts.get("R2", 0))

        # Convert bitscore values to numeric values and update bitscore statistics
        if "bitscore" in chunk.columns:
            bitscore = pd.to_numeric(chunk["bitscore"], errors="coerce").dropna()
            counters["bitscore_sum"] += float(bitscore.sum())
            counters["bitscore_count"] += int(bitscore.shape[0])

        # Convert DIAMOND e-values to numeric values and update e-value statistics
        if "diamond_evalue" in chunk.columns:
            evalues = pd.to_numeric(chunk["diamond_evalue"], errors="coerce").dropna()
            counters["evalue_sum"] += float(evalues.sum())
            counters["evalue_count"] += int(evalues.shape[0])

    # Calculate the mean bitscore if valid bitscore values were found
    mean_bitscore = (
        counters["bitscore_sum"] / counters["bitscore_count"]
        if counters["bitscore_count"] > 0
        else pd.NA
    )

    # Calculate the mean DIAMOND e-value if valid e-values were found
    mean_evalue = (
        counters["evalue_sum"] / counters["evalue_count"]
        if counters["evalue_count"] > 0
        else pd.NA
    )

    # Return all summary statistics as a dictionary
    return {
        "annotated_hit_count": counters["annotated_hit_count"],
        "R1_hit_count": counters["R1_hit_count"],
        "R2_hit_count": counters["R2_hit_count"],
        "mean_bitscore": mean_bitscore,
        "mean_diamond_evalue": mean_evalue,
    }


# Build annotation summary statistics for every detected sample
def build_annotation_summary(samples, strict=False):
    """
    Build annotation summary statistics for all samples.
    """

    # Store one annotation summary row per sample
    rows = []

    # Process samples in alphabetical order for reproducible output
    for sample in sorted(samples, key=natural_sort_key):

        # Read input read counts for this sample, if available
        input_read_counts = read_input_read_counts_for_sample(
            samples=samples,
            sample=sample,
            strict=strict,
        )

        # Get the default annotation table for this sample
        path = samples[sample].get("annotation")

        # Summarise unannotated hits from the unannotated hits table, if available
        unannotated_summary = summarise_unannotated_hits_for_sample(
            samples=samples,
            sample=sample,
            strict=strict,
        )

        # If the annotation table is missing, use missing values for annotation statistics
        if path is None:
            summary = {
                "annotated_hit_count": pd.NA,
                "R1_hit_count": pd.NA,
                "R2_hit_count": pd.NA,
                "mean_bitscore": pd.NA,
                "mean_diamond_evalue": pd.NA,
            }

        else:
            # Try to summarise the annotation table for this sample
            try:
                summary = summarise_annotation_file(path)

            # Handle annotation tables that cannot be read or summarised
            except Exception as exc:

                # In strict mode, stop immediately if summarising fails
                if strict:
                    error_exit(
                        f"Could not summarise annotation table for sample '{sample}': {exc}"
                    )

                # In non-strict mode, warn the user and keep missing values for this sample
                warn_msg(
                    f"Could not summarise annotation table for sample '{sample}': {exc}"
                )

                # Use missing values for annotation statistics if summarising failed
                summary = {
                    "annotated_hit_count": pd.NA,
                    "R1_hit_count": pd.NA,
                    "R2_hit_count": pd.NA,
                    "mean_bitscore": pd.NA,
                    "mean_diamond_evalue": pd.NA,
                }

        # Get input read counts
        total_input_reads = input_read_counts["total_input_reads"]
        R1_input_reads = input_read_counts["R1_input_reads"]
        R2_input_reads = input_read_counts["R2_input_reads"]

        # Get annotated hit counts
        annotated_hit_count = summary.get("annotated_hit_count", pd.NA)
        R1_annotated_hit_count = summary.get("R1_hit_count", pd.NA)
        R2_annotated_hit_count = summary.get("R2_hit_count", pd.NA)

        # Get unannotated hit counts
        unannotated_hit_count = unannotated_summary["unannotated_hit_count"]
        R1_unannotated_hit_count = unannotated_summary["R1_unannotated_hit_count"]
        R2_unannotated_hit_count = unannotated_summary["R2_unannotated_hit_count"]

        # Calculate the total number of DIAMOND hits retained by FungiFunk
        total_diamond_hit_count = add_count_values(
            annotated_hit_count,
            unannotated_hit_count,
        )

        # Calculate the number of R1 DIAMOND hits retained by FungiFunk
        R1_diamond_hit_count = add_count_values(
            R1_annotated_hit_count,
            R1_unannotated_hit_count,
        )

        # Calculate the number of R2 DIAMOND hits retained by FungiFunk
        R2_diamond_hit_count = add_count_values(
            R2_annotated_hit_count,
            R2_unannotated_hit_count,
        )

        # Calculate retained DIAMOND hits per input read
        diamond_hits_per_input_read = divide_count_values(
            total_diamond_hit_count,
            total_input_reads,
        )

        # Calculate retained DIAMOND hits per R1 input read
        R1_diamond_hits_per_input_read = divide_count_values(
            R1_diamond_hit_count,
            R1_input_reads,
        )

        # Calculate retained DIAMOND hits per R2 input read
        R2_diamond_hits_per_input_read = divide_count_values(
            R2_diamond_hit_count,
            R2_input_reads,
        )

        # Calculate the fraction of retained DIAMOND hits that were annotated
        annotated_hit_fraction = divide_count_values(
            annotated_hit_count,
            total_diamond_hit_count,
        )

        # Calculate the fraction of retained DIAMOND hits that were unannotated
        unannotated_hit_fraction = divide_count_values(
            unannotated_hit_count,
            total_diamond_hit_count,
        )

        # Build the final row in an interpretation-friendly column order
        row = {
            "sample": sample,
            "total_input_reads": total_input_reads,
            "total_diamond_hit_count": total_diamond_hit_count,
            "diamond_hits_per_input_read": diamond_hits_per_input_read,
            "annotated_hit_count": annotated_hit_count,
            "annotated_hit_fraction": annotated_hit_fraction,
            "unannotated_hit_count": unannotated_hit_count,
            "unannotated_hit_fraction": unannotated_hit_fraction,
            "R1_input_reads": R1_input_reads,
            "R1_diamond_hit_count": R1_diamond_hit_count,
            "R1_diamond_hits_per_input_read": R1_diamond_hits_per_input_read,
            "R2_input_reads": R2_input_reads,
            "R2_diamond_hit_count": R2_diamond_hit_count,
            "R2_diamond_hits_per_input_read": R2_diamond_hits_per_input_read,
            "mean_bitscore": summary.get("mean_bitscore", pd.NA),
            "mean_diamond_evalue": summary.get("mean_diamond_evalue", pd.NA),
        }

        # Store the completed summary row
        rows.append(row)

    # Convert all summary rows into a pandas table with a fixed column order
    return pd.DataFrame(
        rows,
        columns=[
            "sample",
            "total_input_reads",
            "total_diamond_hit_count",
            "diamond_hits_per_input_read",
            "annotated_hit_count",
            "annotated_hit_fraction",
            "unannotated_hit_count",
            "unannotated_hit_fraction",
            "R1_input_reads",
            "R1_diamond_hit_count",
            "R1_diamond_hits_per_input_read",
            "R2_input_reads",
            "R2_diamond_hit_count",
            "R2_diamond_hits_per_input_read",
            "mean_bitscore",
            "mean_diamond_evalue",
        ],
    )


# Build collapsible diagnostic plot sections for the annotation summary
def build_annotation_diagnostic_sections(annotation_summary):
    """
    Build per-sample diagnostic plots from the annotation summary table.
    """

    # Return no diagnostic sections if the annotation summary is empty
    if annotation_summary is None or annotation_summary.empty:
        return []

    # Define the diagnostic plots shown below the annotation summary table
    diagnostic_configs = [
        {
            "plot_type": "single",
            "id": "annotation-diagnostic-total-input-reads",
            "title": "Total input reads per sample",
            "description": (
                "Shows the total number of input reads provided to FungiFunk for each sample."
            ),
            "y_col": "total_input_reads",
            "plot_title": "Total input reads per sample",
        },
        {
            "plot_type": "single",
            "id": "annotation-diagnostic-diamond-hits-per-input-read",
            "title": "DIAMOND hits per input read",
            "description": (
                "Shows the number of retained DIAMOND hits relative to the total number of input reads."
            ),
            "y_col": "diamond_hits_per_input_read",
            "plot_title": "DIAMOND hits per input read",
        },
        {
            "plot_type": "single",
            "id": "annotation-diagnostic-annotated-hit-fraction",
            "title": "Annotated hit fraction",
            "description": (
                "Shows the fraction of retained DIAMOND hits that received functional annotation."
            ),
            "y_col": "annotated_hit_fraction",
            "plot_title": "Annotated hit fraction",
        },
        {
            "plot_type": "grouped",
            "id": "annotation-diagnostic-r1-r2-input-reads",
            "title": "R1 and R2 input reads per sample",
            "description": (
                "Shows the number of input reads separately for R1 and R2 in each sample."
            ),
            "y_cols": ["R1_input_reads", "R2_input_reads"],
            "y_labels": {
                "R1_input_reads": "R1 input reads",
                "R2_input_reads": "R2 input reads",
            },
            "color_map": {
                "R1_input_reads": "#69b3e7",
                "R2_input_reads": "#f28e2b",
            },
            "y_axis_title": "Input reads",
            "plot_title": "R1 and R2 input reads per sample",
        },
        {
            "plot_type": "grouped",
            "id": "annotation-diagnostic-r1-r2-diamond-hits-per-input-read",
            "title": "R1 and R2 DIAMOND hits per input read",
            "description": (
                "Shows the number of retained DIAMOND hits per input read separately for R1 and R2 in each sample."
            ),
            "y_cols": [
                "R1_diamond_hits_per_input_read",
                "R2_diamond_hits_per_input_read",
            ],
            "y_labels": {
                "R1_diamond_hits_per_input_read": "R1 DIAMOND hits per input read",
                "R2_diamond_hits_per_input_read": "R2 DIAMOND hits per input read",
            },
            "color_map": {
                "R1_diamond_hits_per_input_read": "#69b3e7",
                "R2_diamond_hits_per_input_read": "#f28e2b",
            },
            "y_axis_title": "DIAMOND hits per input read",
            "plot_title": "R1 and R2 DIAMOND hits per input read",
        },
    ]

    # Store all diagnostic sections that can be created
    diagnostic_sections = []

    # Build each diagnostic section
    for config in diagnostic_configs:

        # Build simple one-metric barplots
        if config["plot_type"] == "single":

            # Skip this plot if the required columns are not available
            if (
                "sample" not in annotation_summary.columns
                or config["y_col"] not in annotation_summary.columns
            ):
                continue

            # Keep only the sample column and the metric used in this diagnostic plot
            plot_data = annotation_summary[["sample", config["y_col"]]].copy()

            # Sort samples in natural alphanumeric order before plotting
            plot_data = sort_table_by_sample(plot_data)

            # Convert the plotted metric to numeric values
            plot_data[config["y_col"]] = pd.to_numeric(
                plot_data[config["y_col"]],
                errors="coerce",
            )

            # Remove samples without a valid value for this metric
            plot_data = plot_data.dropna(subset=[config["y_col"]])

            # Skip this plot if no valid values remain
            if plot_data.empty:
                continue

            # Create the barplot HTML using the existing plotting helper
            plot_html = make_total_barplot_html(
                plot_data,
                title="",
                x_col="sample",
                y_col=config["y_col"],
            )

        # Build grouped two-metric barplots
        elif config["plot_type"] == "grouped":

            # Define the columns required for this grouped diagnostic plot
            required_columns = ["sample"] + config["y_cols"]

            # Skip this plot if any required column is missing
            if any(column not in annotation_summary.columns for column in required_columns):
                continue

            # Keep only the columns needed for this grouped diagnostic plot
            plot_data = annotation_summary[required_columns].copy()

            # Create the grouped barplot HTML
            plot_html = make_grouped_barplot_html(
                plot_data,
                title="",
                x_col="sample",
                y_cols=config["y_cols"],
                y_labels=config["y_labels"],
                y_axis_title=config["y_axis_title"],
                color_map=config["color_map"],
            )

            # Skip this plot if no valid grouped plot could be created
            if not plot_html:
                continue

        # Skip unknown plot types
        else:
            continue

        # Store the completed diagnostic section
        diagnostic_sections.append(
            {
                "id": config["id"],
                "title": config["title"],
                "description": config["description"],
                "plot_title": config["plot_title"],
                "plot_html": plot_html,
            }
        )

    # Return all diagnostic sections
    return diagnostic_sections


# Build a compact table with the number of detected features per sample
def build_detected_features_summary(raw_matrices, samples):
    """
    Build a compact table with the number of detected features per sample.
    """

    # Store one detected-feature summary row per sample
    rows = []

    # Process samples in alphabetical order for reproducible output
    for sample in sorted(samples, key=natural_sort_key):

        # Count detected features from the already-created count matrices
        row = {
            "sample": sample,
            "proteins": count_detected_features_from_matrix(raw_matrices.get("protein"), sample),
            "KEGG_ko": count_detected_features_from_matrix(raw_matrices.get("ko"), sample),
            "CAZy": count_detected_features_from_matrix(raw_matrices.get("cazy"), sample),
            "PFAMs": count_detected_features_from_matrix(raw_matrices.get("pfam"), sample),
            "eggNOG_OGs": count_detected_features_from_matrix(raw_matrices.get("eggnog_ogs"), sample),
        }

        # Store the completed row
        rows.append(row)

    # Convert all detected-feature summary rows into a pandas table
    return pd.DataFrame(rows)


# Count how many features are detected in one sample from a count matrix
def count_detected_features_from_matrix(matrix, sample):
    """
    Count how many features have a non-zero count for one sample in one matrix.
    """
    # Return missing if the matrix is unavailable, empty, or does not contain this sample
    if matrix is None or matrix.empty or sample not in matrix.columns:
        return pd.NA

    # Count features with counts greater than zero for this sample
    return int((matrix[sample] > 0).sum())


# Build collapsible diagnostic plot sections for the detected features summary
def build_detected_features_diagnostic_sections(detected_features_summary):
    """
    Build per-sample diagnostic plots from the detected features summary table.
    """

    # Return no diagnostic sections if the detected features summary is empty
    if detected_features_summary is None or detected_features_summary.empty:
        return []

    # Define the diagnostic plots shown below the detected features summary table
    diagnostic_configs = [
        {
            "id": "detected-features-diagnostic-proteins",
            "title": "Detected proteins per sample",
            "description": (
                "Shows the number of detected protein features in each sample."
            ),
            "y_col": "proteins",
            "plot_title": "Detected proteins per sample",
        },
        {
            "id": "detected-features-diagnostic-kegg-ko",
            "title": "Detected KEGG orthologs per sample",
            "description": (
                "Shows the number of detected KEGG ortholog features in each sample."
            ),
            "y_col": "KEGG_ko",
            "plot_title": "Detected KEGG orthologs per sample",
        },
        {
            "id": "detected-features-diagnostic-cazy",
            "title": "Detected CAZy families per sample",
            "description": (
                "Shows the number of detected CAZy family features in each sample."
            ),
            "y_col": "CAZy",
            "plot_title": "Detected CAZy families per sample",
        },
        {
            "id": "detected-features-diagnostic-pfams",
            "title": "Detected PFAM features per sample",
            "description": (
                "Shows the number of detected PFAM features in each sample."
            ),
            "y_col": "PFAMs",
            "plot_title": "Detected PFAM features per sample",
        },
        {
            "id": "detected-features-diagnostic-eggnog-ogs",
            "title": "Detected eggNOG OGs per sample",
            "description": (
                "Shows the number of detected eggNOG orthologous group features in each sample."
            ),
            "y_col": "eggNOG_OGs",
            "plot_title": "Detected eggNOG OGs per sample",
        },
    ]

    # Store all diagnostic sections that can be created
    diagnostic_sections = []

    # Build each diagnostic section
    for config in diagnostic_configs:

        # Skip this plot if the required columns are not available
        if (
            "sample" not in detected_features_summary.columns
            or config["y_col"] not in detected_features_summary.columns
        ):
            continue

        # Keep only the sample column and the metric used in this diagnostic plot
        plot_data = detected_features_summary[["sample", config["y_col"]]].copy()

        # Sort samples in natural alphanumeric order before plotting
        plot_data = sort_table_by_sample(plot_data)

        # Convert the plotted metric to numeric values
        plot_data[config["y_col"]] = pd.to_numeric(
            plot_data[config["y_col"]],
            errors="coerce",
        )

        # Remove samples without a valid value for this metric
        plot_data = plot_data.dropna(subset=[config["y_col"]])

        # Skip this plot if no valid values remain
        if plot_data.empty:
            continue

        # Create the barplot HTML using the existing plotting helper
        plot_html = make_total_barplot_html(
            plot_data,
            title="",
            x_col="sample",
            y_col=config["y_col"],
        )

        # Store the completed diagnostic section
        diagnostic_sections.append(
            {
                "id": config["id"],
                "title": config["title"],
                "description": config["description"],
                "plot_title": config["plot_title"],
                "plot_html": plot_html,
            }
        )

    # Return all diagnostic sections
    return diagnostic_sections


# Select the most abundant functional features for display in the HTML report
def select_top_count_matrix(matrix, top_n=50, no_top=False):
    """
    Select top functional features by total count across samples.
    """
    # Return an empty table if the matrix is missing or has no data
    if matrix is None or matrix.empty:
        return pd.DataFrame()

    # If no_top is enabled, keep the full matrix without filtering
    if no_top:
        return matrix

    # Calculate the total count of each feature across all samples
    totals = matrix.sum(axis=1).sort_values(ascending=False)

    # Select the feature IDs with the highest total counts
    selected_ids = totals.head(top_n).index

    # Return the matrix subset containing only the selected top features
    return matrix.loc[selected_ids]


# Extract species and accession information from a FungiFunk protein ID
def parse_protein_id_metadata(protein_id):
    """
    Extract species and accession from a FungiFunk protein ID.
    """

    # Return missing values if the protein ID is missing
    if is_missing_value(protein_id):
        return "-", "-"

    # Convert the protein ID to text
    protein_id = str(protein_id)

    # Split the protein ID using the FungiFunk separator
    parts = protein_id.split("|")

    # Extract species if the expected protein ID structure is available
    species = parts[1].replace("_", " ") if len(parts) >= 2 else "-"

    # Extract the final accession if the expected protein ID structure is available
    accession = parts[-1] if len(parts) >= 3 else "-"

    # Return the extracted species and accession
    return species, accession


# Extract the broad CAZy class from a CAZy family ID
def parse_cazy_class(cazy_id):
    """
    Extract the broad CAZy class from a CAZy family ID.
    """

    # Return missing if the CAZy ID is missing
    if is_missing_value(cazy_id):
        return "-"

    # Convert the CAZy ID to text and remove extra spaces
    cazy_id = str(cazy_id).strip()

    # Check the longest prefix first because CBM has three letters
    for prefix in ["CBM", "GH", "GT", "PL", "CE", "AA"]:

        # Return the matching CAZy class prefix
        if cazy_id.startswith(prefix):
            return prefix

    # Return missing if no known CAZy class was detected
    return "-"


# Join associated annotation values into a compact text field
def join_limited_values(values, max_values=MAX_ASSOCIATED_VALUES):
    """
    Join unique values into a compact comma-separated text field.
    """

    # Remove missing values and repeated values while preserving sorted output
    clean_values = sorted(
        {
            str(value).strip()
            for value in values
            if not is_missing_value(value)
        }
    )

    # Return the standard missing-value label if no real values are available
    if not clean_values:
        return "-"

    # Keep only a limited number of values to avoid very long cells
    selected_values = clean_values[:max_values]

    # Add a short suffix if some values were not shown
    if len(clean_values) > max_values:
        selected_values.append(f"...(+{len(clean_values) - max_values} more)")

    # Return the compact text field
    return ", ".join(selected_values)


# Select feature IDs that will be displayed in the HTML report
def get_html_feature_ids(matrix, html_top_n=50, html_no_top=False):
    """
    Select feature IDs that need metadata for the HTML report.
    """

    # Return an empty set if the matrix is missing or empty
    if matrix is None or matrix.empty:
        return set()

    # Use all feature IDs only when top filtering is disabled
    if html_no_top:
        return set(matrix.index.astype(str))

    # Rank features by total count across samples
    totals = matrix.sum(axis=1).sort_values(ascending=False)

    # Return only the top feature IDs displayed in the HTML report
    return set(totals.head(html_top_n).index.astype(str))


# Build simple metadata for selected features when no extra context is requested
def build_simple_feature_metadata(id_col, target_features):
    """
    Build a minimal metadata table containing only selected feature IDs.
    """

    # Return an empty table with the expected ID column if no features are available
    if not target_features:
        return pd.DataFrame(columns=[id_col])

    # Store selected feature IDs in sorted order
    rows = [{id_col: feature} for feature in sorted(target_features)]

    # Return the simple feature metadata table
    return pd.DataFrame(rows, columns=[id_col])


# Build human-readable metadata for protein features shown in the HTML report
def build_protein_feature_metadata(protein_metadata, target_proteins):
    """
    Build metadata for selected protein features.
    """

    # Define the output columns
    output_columns = [
        "protein_id",
        "species",
        "accession",
        "Preferred_name",
        "Description",
    ]

    # Return an empty table if no metadata or target proteins are available
    if protein_metadata is None or protein_metadata.empty or not target_proteins:
        return pd.DataFrame(columns=output_columns)

    # Keep only selected proteins
    df = protein_metadata[
        protein_metadata["protein_id"].astype(str).isin(target_proteins)
    ].copy()

    # Add missing metadata columns if needed
    for col in ["Preferred_name", "Description"]:
        if col not in df.columns:
            df[col] = "-"

    # Extract species and accession from protein IDs
    parsed_values = df["protein_id"].apply(parse_protein_id_metadata)

    # Store parsed species
    df["species"] = parsed_values.apply(lambda value: value[0])

    # Store parsed accession
    df["accession"] = parsed_values.apply(lambda value: value[1])

    # Keep the final clean output columns
    df = df[output_columns].fillna("-")

    # Return metadata sorted by protein ID
    return df.sort_values("protein_id")


# Build contextual metadata for selected KEGG KO features
def build_ko_feature_metadata(protein_metadata, target_kos):
    """
    Build metadata for selected KEGG KO features using protein metadata.
    """

    # Define the output columns
    output_columns = [
        "KEGG_ko",
        "representative_preferred_name",
        "representative_description",
        "associated_EC",
        "associated_KEGG_Pathway",
    ]

    # Return an empty table if no metadata or target KOs are available
    if protein_metadata is None or protein_metadata.empty or not target_kos:
        return pd.DataFrame(columns=output_columns)

    # Store preferred names and descriptions associated with each KO
    context_counts = defaultdict(lambda: defaultdict(int))

    # Store EC and KEGG pathway values associated with each KO
    associated_ec = defaultdict(set)
    associated_pathway = defaultdict(set)

    # Add missing columns if needed
    for col in ["KEGG_ko", "Preferred_name", "Description", "EC", "KEGG_Pathway"]:
        if col not in protein_metadata.columns:
            protein_metadata[col] = "-"

    # Process protein metadata row by row to keep RAM usage low
    for row in protein_metadata[
        ["KEGG_ko", "Preferred_name", "Description", "EC", "KEGG_Pathway"]
    ].itertuples(index=False):

        # Extract values from the current row
        ko_values, preferred_name, description, ec_values, pathway_values = row

        # Split multi-value KO annotations
        for ko in split_annotation_values(ko_values):

            # Keep only KOs selected for the HTML report
            if ko not in target_kos:
                continue

            # Use preferred name and description as the context
            context = (
                first_valid_value([preferred_name]),
                first_valid_value([description]),
            )

            # Count how often this context is associated with this KO
            context_counts[ko][context] += 1

            # Store associated EC values
            associated_ec[ko].update(split_annotation_values(ec_values))

            # Store associated KEGG pathway values
            associated_pathway[ko].update(split_annotation_values(pathway_values))

    # Store one metadata row per KO
    rows = []

    # Build output rows in sorted order
    for ko in sorted(target_kos):

        # Select the most frequent associated context
        if context_counts[ko]:
            top_context, _ = max(
                context_counts[ko].items(),
                key=lambda item: item[1],
            )
            representative_preferred_name, representative_description = top_context
        else:
            representative_preferred_name = "-"
            representative_description = "-"

        # Store the KO metadata row
        rows.append(
            {
                "KEGG_ko": ko,
                "representative_preferred_name": representative_preferred_name,
                "representative_description": representative_description,
                "associated_EC": join_limited_values(associated_ec[ko]),
                "associated_KEGG_Pathway": join_limited_values(associated_pathway[ko]),
            }
        )

    # Return the KO metadata table
    return pd.DataFrame(rows, columns=output_columns)


# Build metadata for selected CAZy features
def build_cazy_feature_metadata(target_cazys):
    """
    Build metadata for selected CAZy features.
    """

    # Define the output columns
    output_columns = [
        "CAZy",
        "CAZy_class",
        "CAZy_class_description",
    ]

    # Store one metadata row per CAZy feature
    rows = []

    # Process selected CAZy IDs in sorted order
    for cazy_id in sorted(target_cazys):

        # Extract the broad CAZy class
        cazy_class = parse_cazy_class(cazy_id)

        # Get the class description
        cazy_class_description = CAZY_CLASS_DESCRIPTIONS.get(cazy_class, "-")

        # Store the CAZy metadata row
        rows.append(
            {
                "CAZy": cazy_id,
                "CAZy_class": cazy_class,
                "CAZy_class_description": cazy_class_description,
            }
        )

    # Return the CAZy metadata table
    return pd.DataFrame(rows, columns=output_columns)


# Build derived context metadata for selected PFAM and eggNOG OG features in one annotation-table pass
def build_associated_context_metadata_tables(
    samples,
    feature_targets,
    strict=False,
    chunksize=100_000,
):
    """
    Build derived context metadata for selected PFAM and eggNOG OG features.

    This is not official database metadata. It reports the most frequent protein
    annotation associated with each selected feature in this FungiFunk run.

    The annotation tables are scanned only once to reduce runtime.
    """

    # Map internal table keys to annotation-table column names
    feature_columns = {
        "pfam": "PFAMs",
        "eggnog_ogs": "eggNOG_OGs",
    }

    # Define output columns for each context metadata table
    output_columns = {
        "pfam": [
            "PFAMs",
            "top_associated_preferred_name",
            "top_associated_description",
            "top_associated_hit_count",
            "top_associated_hit_fraction",
        ],
        "eggnog_ogs": [
            "eggNOG_OGs",
            "top_associated_preferred_name",
            "top_associated_description",
            "top_associated_hit_count",
            "top_associated_hit_fraction",
        ],
    }

    # Count total feature-associated hits for each feature type
    associated_hit_counts = {
        table_key: defaultdict(int)
        for table_key in feature_columns
    }

    # Count associated protein annotation contexts for each feature type
    context_counts = {
        table_key: defaultdict(lambda: defaultdict(int))
        for table_key in feature_columns
    }

    # Keep only feature types that have selected target features
    active_feature_targets = {
        table_key: set(targets)
        for table_key, targets in feature_targets.items()
        if targets
    }

    # Return empty metadata tables if no target features were selected
    if not active_feature_targets:
        return {
            table_key: pd.DataFrame(columns=output_columns[table_key])
            for table_key in feature_columns
        }

    # Process annotation files sample by sample to keep RAM usage low
    for sample in sorted(samples, key=natural_sort_key):

        # Get the default annotation table for this sample
        path = samples[sample].get("annotation")

        # Skip samples without annotation table
        if path is None:
            continue

        # Read only the header first
        header = read_header(path)

        # Define columns needed from this annotation file
        wanted_cols = ["Preferred_name", "Description"]

        # Add available feature columns that have selected target features
        for table_key in active_feature_targets:
            feature_col = feature_columns[table_key]

            if feature_col in header:
                wanted_cols.append(feature_col)

        # Keep only columns that actually exist in the file
        usecols = [col for col in wanted_cols if col in header]

        # Skip this file if no selected feature column is available
        if not any(feature_columns[key] in usecols for key in active_feature_targets):
            continue

        # Try to process this annotation table in chunks
        try:
            for chunk in pd.read_csv(
                path,
                sep="\t",
                usecols=usecols,
                chunksize=chunksize,
                dtype=str,
                na_filter=False,
            ):

                # Add missing context columns if they are absent
                for col in ["Preferred_name", "Description"]:
                    if col not in chunk.columns:
                        chunk[col] = "-"

                # Process each active feature type from the same chunk
                for table_key, target_features in active_feature_targets.items():

                    # Get the annotation column for this feature type
                    feature_col = feature_columns[table_key]

                    # Skip this feature type if the column is not present in this chunk
                    if feature_col not in chunk.columns:
                        continue

                    # Keep rows where the feature column contains at least one value
                    feature_chunk = chunk[
                        ~chunk[feature_col].apply(is_missing_value)
                    ]

                    # Process only the required columns row by row
                    for feature_values, preferred_name, description in feature_chunk[
                        [feature_col, "Preferred_name", "Description"]
                    ].itertuples(index=False, name=None):

                        # Clean the associated context values
                        preferred_name = first_valid_value([preferred_name])
                        description = first_valid_value([description])

                        # Use preferred name and description as the associated context
                        context = (preferred_name, description)

                        # Split multi-value feature fields into clean annotation values
                        features = split_annotation_values(feature_values)

                        # For eggNOG OGs, keep only the most specific valid OG
                        if table_key == "eggnog_ogs" and features:
                            features = [features[-1]]

                        # Process the selected feature values
                        for feature in features:

                            # Keep only features selected for the HTML report
                            if feature not in target_features:
                                continue

                            # Count one associated hit for this feature
                            associated_hit_counts[table_key][feature] += 1

                            # Count how often this context is associated with this feature
                            context_counts[table_key][feature][context] += 1

        # Handle annotation tables that cannot be processed
        except Exception as exc:

            # In strict mode, stop immediately
            if strict:
                error_exit(
                    f"Could not build feature context metadata from annotation table for sample '{sample}': {exc}"
                )

            # In non-strict mode, warn and continue
            warn_msg(
                f"Could not build feature context metadata from annotation table for sample '{sample}': {exc}"
            )

    # Store final metadata tables
    metadata_tables = {}

    # Build one metadata table per feature type
    for table_key, feature_col in feature_columns.items():

        # Store one metadata row per selected feature
        rows = []

        # Get selected target features for this table type
        target_features = feature_targets.get(table_key, set())

        # Build output rows in sorted order
        for feature in sorted(target_features):

            # Get the total number of associated hits for this feature
            associated_hit_count = associated_hit_counts[table_key].get(feature, 0)

            # Select the most frequent associated context
            if context_counts[table_key][feature]:
                top_context, top_count = max(
                    context_counts[table_key][feature].items(),
                    key=lambda item: item[1],
                )
                top_preferred_name, top_description = top_context
            else:
                top_count = 0
                top_preferred_name = "-"
                top_description = "-"

            # Calculate the fraction of hits represented by the top context
            top_fraction = (
                round(top_count / associated_hit_count, 4)
                if associated_hit_count > 0
                else pd.NA
            )

            # Store the metadata row
            rows.append(
                {
                    feature_col: feature,
                    "top_associated_preferred_name": top_preferred_name,
                    "top_associated_description": top_description,
                    "top_associated_hit_count": top_count,
                    "top_associated_hit_fraction": top_fraction,
                }
            )

        # Store this metadata table
        metadata_tables[table_key] = pd.DataFrame(
            rows,
            columns=output_columns[table_key],
        )

    # Return all context metadata tables
    return metadata_tables


# Select all feature IDs that need metadata from raw and CPM matrices
def get_html_feature_ids_from_raw_and_cpm(
    raw_matrix,
    cpm_matrix=None,
    html_top_n=50,
    html_no_top=False,
):
    """
    Select feature IDs that need metadata for raw and CPM HTML tables.
    """

    # Select feature IDs from the raw matrix
    feature_ids = get_html_feature_ids(
        raw_matrix,
        html_top_n=html_top_n,
        html_no_top=html_no_top,
    )

    # Add feature IDs from the CPM matrix, if available
    if cpm_matrix is not None:
        feature_ids.update(
            get_html_feature_ids(
                cpm_matrix,
                html_top_n=html_top_n,
                html_no_top=html_no_top,
            )
        )

    # Return the combined feature ID set
    return feature_ids


# Build and write feature metadata tables used to enrich the HTML report
def build_feature_metadata_tables(
    samples,
    raw_matrices,
    protein_metadata,
    feature_metadata_dir,
    cpm_matrices=None,
    html_top_n=50,
    html_no_top=False,
    calculate_pfam_eggnog_context=False,
    strict=False,
    time_stamps=False,
):
    """
    Build metadata tables for functional features shown in the HTML report.
    """

    # Create the feature metadata output directory
    feature_metadata_dir.mkdir(parents=True, exist_ok=True)

    # Store metadata tables by feature type
    metadata_tables = {}

    # Use an empty dictionary if CPM matrices are not available
    if cpm_matrices is None:
        cpm_matrices = {}

    # Select top features from raw and CPM matrices
    target_proteins = get_html_feature_ids_from_raw_and_cpm(
        raw_matrices.get("protein"),
        cpm_matrices.get("protein"),
        html_top_n=html_top_n,
        html_no_top=html_no_top,
    )

    target_kos = get_html_feature_ids_from_raw_and_cpm(
        raw_matrices.get("ko"),
        cpm_matrices.get("ko"),
        html_top_n=html_top_n,
        html_no_top=html_no_top,
    )

    target_cazys = get_html_feature_ids_from_raw_and_cpm(
        raw_matrices.get("cazy"),
        cpm_matrices.get("cazy"),
        html_top_n=html_top_n,
        html_no_top=html_no_top,
    )

    target_pfams = get_html_feature_ids_from_raw_and_cpm(
        raw_matrices.get("pfam"),
        cpm_matrices.get("pfam"),
        html_top_n=html_top_n,
        html_no_top=html_no_top,
    )

    target_eggnog_ogs = get_html_feature_ids_from_raw_and_cpm(
        raw_matrices.get("eggnog_ogs"),
        cpm_matrices.get("eggnog_ogs"),
        html_top_n=html_top_n,
        html_no_top=html_no_top,
    )

    # Build protein feature metadata
    metadata_tables["protein"] = build_protein_feature_metadata(
        protein_metadata=protein_metadata,
        target_proteins=target_proteins,
    )

    # Write protein feature metadata
    protein_feature_metadata_path = feature_metadata_dir / "protein_feature_metadata.tsv"
    metadata_tables["protein"].to_csv(protein_feature_metadata_path, sep="\t", index=False)
    log_msg(f"    (Written to: {protein_feature_metadata_path})", time_stamps)

    # Build KO feature metadata
    metadata_tables["ko"] = build_ko_feature_metadata(
        protein_metadata=protein_metadata,
        target_kos=target_kos,
    )

    # Write KO feature metadata
    ko_metadata_path = feature_metadata_dir / "ko_metadata.tsv"
    metadata_tables["ko"].to_csv(ko_metadata_path, sep="\t", index=False)
    log_msg(f"    (Written to: {ko_metadata_path})", time_stamps)

    # Build CAZy feature metadata
    metadata_tables["cazy"] = build_cazy_feature_metadata(
        target_cazys=target_cazys,
    )

    # Write CAZy feature metadata
    cazy_metadata_path = feature_metadata_dir / "cazy_metadata.tsv"
    metadata_tables["cazy"].to_csv(cazy_metadata_path, sep="\t", index=False)
    log_msg(f"    (Written to: {cazy_metadata_path})", time_stamps)

    # Calculate PFAM and eggNOG OG context only when requested by the user
    if calculate_pfam_eggnog_context:

        # Tell the user that the optional context step is being run
        log_msg(
            "Calculating PFAM and eggNOG OG context from annotation tables...",
            time_stamps,
        )

        # Build PFAM and eggNOG OG context metadata in a single annotation-table pass
        context_metadata_tables = build_associated_context_metadata_tables(
            samples=samples,
            feature_targets={
                "pfam": target_pfams,
                "eggnog_ogs": target_eggnog_ogs,
            },
            strict=strict,
        )

        # Store PFAM context metadata
        metadata_tables["pfam"] = context_metadata_tables["pfam"]

        # Store eggNOG OG context metadata
        metadata_tables["eggnog_ogs"] = context_metadata_tables["eggnog_ogs"]

    # Otherwise, write minimal PFAM and eggNOG OG metadata without scanning annotation tables
    else:

        # Tell the user that the optional context step is being skipped
        log_msg(
            "Skipping PFAM and eggNOG OG context calculation. Use --calculate_pfam_eggnog_context to enable it.",
            time_stamps,
        )

        # Build minimal PFAM metadata containing only selected PFAM IDs
        metadata_tables["pfam"] = build_simple_feature_metadata(
            id_col="PFAMs",
            target_features=target_pfams,
        )

        # Build minimal eggNOG OG metadata containing only selected eggNOG OG IDs
        metadata_tables["eggnog_ogs"] = build_simple_feature_metadata(
            id_col="eggNOG_OGs",
            target_features=target_eggnog_ogs,
        )

    # Write PFAM feature metadata
    pfam_metadata_path = feature_metadata_dir / "pfam_metadata.tsv"
    metadata_tables["pfam"].to_csv(pfam_metadata_path, sep="\t", index=False)
    log_msg(f"    (Written to: {pfam_metadata_path})", time_stamps)

    # Write eggNOG OG feature metadata
    eggnog_metadata_path = feature_metadata_dir / "eggnog_ogs_metadata.tsv"
    metadata_tables["eggnog_ogs"].to_csv(eggnog_metadata_path, sep="\t", index=False)
    log_msg(f"    (Written to: {eggnog_metadata_path})", time_stamps)

    # Return metadata tables so they can be reused in the HTML report
    return metadata_tables


# Convert a pandas table into HTML so it can be inserted into the summary report
def dataframe_to_html(df, max_rows=None):
    """
    Convert a DataFrame to an HTML table.
    """
    # If the table is missing or empty, show a simple message in the report
    if df is None or df.empty:
        return "<p>No data available.</p>"

    # Work on a copy so the original table is not modified
    display_df = df.copy()

    # If a row limit is provided, keep only the first rows for display
    if max_rows is not None and len(display_df) > max_rows:
        display_df = display_df.head(max_rows)

    # Convert the table to HTML using a CSS class for report styling
    return display_df.to_html(
        index=False,
        escape=False,
        classes=["summary-table", "paginated-table"],
        border=0,
    )


# Create a small summary table for selected functional features shown in the HTML report
def matrix_top_table(
    matrix,
    top_matrix,
    table_key,
    feature_metadata=None,
    total_col="total_count",
):
    """
    Create a small summary table for selected top features.
    """

    # Return an empty table if the matrix is missing or empty
    if matrix is None or matrix.empty or top_matrix is None or top_matrix.empty:
        return pd.DataFrame()

    # Get the feature ID column for this count table
    id_col = COUNT_TABLES[table_key]["id_col"]

    # Calculate total values across samples for selected features
    totals = matrix.loc[top_matrix.index].sum(axis=1)

    # Create the base top-feature table
    top_table = totals.reset_index()

    # Rename columns to make the table clearer
    top_table.columns = [id_col, total_col]

    # Get metadata for this feature type, if available
    metadata = None
    if feature_metadata is not None:
        metadata = feature_metadata.get(table_key)

    # Add metadata columns if metadata is available
    if metadata is not None and not metadata.empty and id_col in metadata.columns:
        top_table = top_table.merge(metadata, on=id_col, how="left")

    # Sort by the selected total column
    top_table = top_table.sort_values(total_col, ascending=False)

    # Return the enriched top-feature table
    return top_table


# Create an interactive Plotly barplot for the HTML report
def make_total_barplot_html(table, title, x_col, y_col):
    """
    Create a Plotly barplot and return it as HTML.
    """
    # If the table is missing or empty, show a simple message instead of a plot
    if table is None or table.empty:
        return "<p>No data available for this barplot.</p>"

    # Import Plotly only when a barplot is needed
    try:
        import plotly.express as px
    # If Plotly is not available, keep the report running and show a clear message
    except ImportError:
        return "<p>Plotly is not installed; barplot could not be generated.</p>"

    # Build the barplot using the selected table columns
    fig = px.bar(
        table,
        x=x_col,
        y=y_col,
        title=title,
    )

    # Add plot titles, size, and margins
    fig.update_layout(
        xaxis_title=x_col,
        yaxis_title=y_col,
        height=500,
        margin=dict(l=80, r=40, t=80, b=180),
    )

    # Rotate x-axis labels to make long feature names easier to read
    fig.update_xaxes(tickangle=45)

    # Return only the plot HTML so it can be embedded inside the full report
    return fig.to_html(full_html=False, include_plotlyjs=False)


# Build a grouped barplot for comparing multiple metrics per sample
def make_grouped_barplot_html(
    table,
    title,
    x_col,
    y_cols,
    y_labels,
    y_axis_title,
    color_map=None,
):
    """
    Build a grouped barplot as an HTML Plotly figure.

    This is useful when each sample should show more than one bar.
    """

    # Return an empty plot if the input table is missing or empty
    if table is None or table.empty:
        return ""

    # Keep only the columns needed for plotting
    plot_data = table[[x_col] + y_cols].copy()

    # Sort samples in natural alphanumeric order when plotting by sample
    if x_col == "sample":
        plot_data = sort_table_by_sample(plot_data)

    # Convert plotted columns to numeric values
    for y_col in y_cols:
        plot_data[y_col] = pd.to_numeric(
            plot_data[y_col],
            errors="coerce",
        )

    # Remove samples where all plotted values are missing
    plot_data = plot_data.dropna(subset=y_cols, how="all")

    # Return an empty plot if no valid values remain
    if plot_data.empty:
        return ""

    # Create an empty Plotly figure
    fig = go.Figure()

    # Add one bar trace per plotted metric
    for y_col in y_cols:

        # Define optional marker settings
        marker_settings = {}

        # Use custom colours if they were provided
        if color_map is not None and y_col in color_map:
            marker_settings["color"] = color_map[y_col]

        # Add the current metric as one grouped bar series
        fig.add_trace(
            go.Bar(
                x=plot_data[x_col],
                y=plot_data[y_col],
                name=y_labels.get(y_col, y_col),
                marker=marker_settings,
            )
        )

    # Format the grouped barplot
    fig.update_layout(
        title=title,
        xaxis_title="Sample" if x_col == "sample" else x_col,
        yaxis_title=y_axis_title,
        barmode="group",
        legend_title_text="",
        template="plotly_white",
        margin=dict(l=60, r=30, t=70, b=90),
    )

    # Rotate sample names to keep them readable
    fig.update_xaxes(tickangle=45)

    # Return the figure as an HTML fragment
    return fig.to_html(
        full_html=False,
        include_plotlyjs=False,
    )


# Build a boxplot showing value distributions across samples
def make_matrix_boxplot_html(
    matrix,
    title,
    feature_col=None,
    value_label="Value",
    feature_label="Feature",
):
    """
    Build a vertical boxplot as an HTML Plotly figure.

    Each box shows the distribution of one selected feature or taxon across samples.
    Only outlier points are shown to keep the report compact.
    """

    # Return an empty plot if the matrix is missing or empty
    if matrix is None or matrix.empty:
        return ""

    # Work on a copy so the original matrix is not modified
    plot_matrix = matrix.copy()

    # If no feature column name was provided, use the matrix index name or a generic name
    if feature_col is None:
        feature_col = plot_matrix.index.name or feature_label.lower()

    # If the feature column is not present, move the matrix index into a regular column
    if feature_col not in plot_matrix.columns:

        # Move the index into the first column
        plot_matrix = plot_matrix.reset_index()

        # Rename the new index column to the requested feature column name
        plot_matrix = plot_matrix.rename(
            columns={plot_matrix.columns[0]: feature_col}
        )

    # Return an empty plot if the feature column is still not available
    if feature_col not in plot_matrix.columns:
        return ""

    # Use all columns except the feature column as sample columns
    sample_columns = [
        column for column in plot_matrix.columns
        if column != feature_col
    ]

    # A boxplot needs at least two samples to show a distribution
    if len(sample_columns) < 2:
        return ""

    # Keep only the feature column and sample value columns
    plot_matrix = plot_matrix[[feature_col] + sample_columns].copy()

    # Convert the matrix from wide format to long format for plotting
    plot_data = plot_matrix.melt(
        id_vars=feature_col,
        value_vars=sample_columns,
        var_name="sample",
        value_name="value",
    )

    # Convert plotted values to numeric values
    plot_data["value"] = pd.to_numeric(
        plot_data["value"],
        errors="coerce",
    )

    # Remove missing values before plotting
    plot_data = plot_data.dropna(subset=["value"])

    # Return an empty plot if no valid values remain
    if plot_data.empty:
        return ""

    # Convert feature or taxon names to text for plotting
    plot_matrix[feature_col] = plot_matrix[feature_col].astype(str)
    plot_data[feature_col] = plot_data[feature_col].astype(str)

    # Keep the feature or taxon order from the input matrix
    feature_order = plot_matrix[feature_col].tolist()

    # Create the vertical boxplot
    fig = go.Figure()

    fig.add_trace(
        go.Box(
            x=plot_data[feature_col],
            y=plot_data["value"],
            boxpoints="outliers",
            customdata=plot_data[["sample"]],
            hovertemplate=(
                f"<b>{feature_label}: %{{x}}</b><br>"
                "Sample: %{customdata[0]}<br>"
                f"{value_label}: %{{y}}"
                "<extra></extra>"
            ),
            name="",
        )
    )

    # Scale plot height for a vertical boxplot
    plot_height = 600

    # Format the boxplot
    fig.update_layout(
        title=title,
        xaxis_title=feature_label,
        yaxis_title=value_label,
        template="plotly_white",
        height=plot_height,
        showlegend=False,
        margin=dict(l=80, r=30, t=70, b=180),
    )

    # Keep the selected features or taxa in the same order as the table/barplot
    fig.update_xaxes(
        categoryorder="array",
        categoryarray=feature_order,
        tickangle=45,
    )

    # Return the figure as an HTML fragment
    return fig.to_html(
        full_html=False,
        include_plotlyjs=False,
    )


# Create an interactive clustered heatmap with a sample dendrogram
def make_clustermap_html(
    matrix,
    title,
    colorbar_title,
    row_label="Feature",
):
    """
    Create a clustered heatmap and return it as HTML.

    The matrix is log1p-transformed before clustering and plotting. Samples are
    shown with a top dendrogram. Rows are reordered by clustering, but their
    dendrogram is not displayed.
    """

    # Return an empty plot if the matrix is missing or empty
    if matrix is None or matrix.empty:
        return ""

    # Import clustering and plotting dependencies only when needed
    try:
        import numpy as np
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
        from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
        from scipy.spatial.distance import pdist

    # Return a clear message if required dependencies are unavailable
    except ImportError:
        return "<p>Clustered heatmap could not be generated because scipy and/or plotly are not installed.</p>"

    # Convert all matrix values to numeric values
    numeric_matrix = matrix.apply(
        pd.to_numeric,
        errors="coerce",
    )

    # Replace infinite values with missing values
    numeric_matrix = numeric_matrix.replace(
        [np.inf, -np.inf],
        np.nan,
    )

    # Remove empty rows and columns
    numeric_matrix = numeric_matrix.dropna(axis=0, how="all")
    numeric_matrix = numeric_matrix.dropna(axis=1, how="all")

    # Replace remaining missing values with zero
    numeric_matrix = numeric_matrix.fillna(0)

    # A clustered heatmap needs at least two rows and two columns
    if numeric_matrix.shape[0] < 2 or numeric_matrix.shape[1] < 2:
        return ""

    # Negative values are not expected, but clip them before log transformation
    numeric_matrix = numeric_matrix.clip(lower=0)

    # Skip the plot if all values are zero
    if numeric_matrix.to_numpy().sum() == 0:
        return ""

    # Transform values to reduce the influence of highly abundant features
    log_matrix = np.log1p(numeric_matrix)

    # Calculate hierarchical clustering for rows and columns
    row_linkage = linkage(
        pdist(log_matrix.values, metric="euclidean"),
        method="average",
    )

    col_linkage = linkage(
        pdist(log_matrix.values.T, metric="euclidean"),
        method="average",
    )

    # Get the clustered row order without drawing the row dendrogram
    row_order = leaves_list(row_linkage)

    # Get the column dendrogram structure for plotting
    col_dendrogram = dendrogram(
        col_linkage,
        no_plot=True,
    )

    # Get the clustered column order
    col_order = col_dendrogram["leaves"]

    # Reorder the original and transformed matrices
    ordered_original_matrix = numeric_matrix.iloc[row_order, col_order]
    ordered_log_matrix = log_matrix.iloc[row_order, col_order]

    # Store ordered row and column labels
    row_names = [str(index) for index in ordered_log_matrix.index]
    column_names = [str(column) for column in ordered_log_matrix.columns]

    # Create numeric positions that match scipy dendrogram leaf positions
    row_positions = [5 + 10 * index for index in range(len(row_names))]
    column_positions = [5 + 10 * index for index in range(len(column_names))]

    # Store hover information for every heatmap cell
    row_grid = np.repeat(
        np.array(row_names, dtype=object)[:, None],
        len(column_names),
        axis=1,
    )

    column_grid = np.repeat(
        np.array(column_names, dtype=object)[None, :],
        len(row_names),
        axis=0,
    )

    customdata = np.empty(
        (len(row_names), len(column_names), 3),
        dtype=object,
    )

    customdata[:, :, 0] = row_grid
    customdata[:, :, 1] = column_grid
    customdata[:, :, 2] = ordered_original_matrix.values

    # Create a subplot layout with a top dendrogram and a heatmap
    fig = make_subplots(
        rows=2,
        cols=1,
        row_heights=[0.22, 0.78],
        vertical_spacing=0.03,
        specs=[
            [{"type": "scatter"}],
            [{"type": "heatmap"}],
        ],
    )

    # Add the column dendrogram above the heatmap
    for x_coordinates, y_coordinates in zip(
        col_dendrogram["icoord"],
        col_dendrogram["dcoord"],
    ):
        fig.add_trace(
            go.Scatter(
                x=x_coordinates,
                y=y_coordinates,
                mode="lines",
                line=dict(
                    color="black",
                    width=1,
                ),
                hoverinfo="skip",
                showlegend=False,
            ),
            row=1,
            col=1,
        )

    # Add the log1p-transformed heatmap
    fig.add_trace(
        go.Heatmap(
            z=ordered_log_matrix.values,
            x=column_positions,
            y=row_positions,
            customdata=customdata,
            colorscale="Viridis",
            colorbar=dict(
                title=f"log1p({colorbar_title})",
            ),
            hovertemplate=(
                f"{row_label}: %{{customdata[0]}}<br>"
                "Sample: %{customdata[1]}<br>"
                f"Original {colorbar_title}: %{{customdata[2]}}<br>"
                "log1p value: %{z}"
                "<extra></extra>"
            ),
        ),
        row=2,
        col=1,
    )

    # Scale plot height based on the number of rows
    plot_height = max(
        550,
        min(
            1800,
            280 + 25 * len(row_names),
        ),
    )

    # Format the full clustered heatmap
    fig.update_layout(
        title=title,
        template="plotly_white",
        height=plot_height,
        showlegend=False,
        margin=dict(l=220, r=120, t=80, b=130),
    )

    # Format the column dendrogram axes
    fig.update_xaxes(
        range=[0, 10 * len(column_names)],
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        row=1,
        col=1,
    )

    fig.update_yaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        row=1,
        col=1,
    )

    # Format the heatmap sample axis
    fig.update_xaxes(
        title_text="Sample",
        range=[0, 10 * len(column_names)],
        tickmode="array",
        tickvals=column_positions,
        ticktext=column_names,
        tickangle=45,
        showgrid=False,
        zeroline=False,
        row=2,
        col=1,
    )

    # Format the heatmap row axis
    fig.update_yaxes(
        title_text=row_label,
        tickmode="array",
        tickvals=row_positions,
        ticktext=row_names,
        range=[10 * len(row_names), 0],
        side="left",
        showgrid=False,
        zeroline=False,
        row=2,
        col=1,
    )

    # Return only the plot HTML so it can be embedded inside the full report
    return fig.to_html(
        full_html=False,
        include_plotlyjs=False,
    )


# Build the HTML note shown above functional feature boxplots
def get_functional_boxplot_note_html():
    """
    Return the HTML explanation shown above functional feature boxplots.
    """

    # Explain what the functional feature boxplots show
    return """
<p class="small-text">
    Boxplots show the distribution of selected features across samples.
    Visible outlier points can be inspected by hovering over them.
</p>
"""


# Build the HTML note shown above functional feature clustered heatmaps
def get_functional_clustermap_note_html():
    """
    Return the HTML explanation shown above functional feature clustered heatmaps.
    """

    # Explain how functional clustered heatmaps are generated
    return """
<p class="small-text">
    The clustered heatmap uses log1p-transformed values to reduce the influence of highly abundant features.
    Samples are clustered with hierarchical clustering using Euclidean distance and average linkage, and are shown with a top dendrogram.
    Features are also reordered by abundance-profile similarity, but their row dendrogram is not shown.
</p>
"""


# Create a relative file path that can be used as a link inside the HTML report
def relative_link(path, outdir):
    """
    Return a path relative to the output directory for HTML links.
    """
    # Try to express the file path relative to the summary output directory
    try:
        return path.relative_to(outdir).as_posix()
    # If the path is outside the output directory, fall back to the full path
    except ValueError:
        return str(path)


# Build the HTML note that explains the annotation summary columns
def get_annotation_summary_note_html():
    """
    Return the HTML explanation block for the annotation summary table.
    """

    # Return a collapsible HTML block using the same style as feature metadata notes
    return """
<details class="feature-info column-guide-details">
    <summary>Column guide</summary>

    <div class="column-guide">
        <dl>
            <dt><code>sample</code></dt>
            <dd>Sample name parsed from the FungiFunk output filenames.</dd>

            <dt><code>total_input_reads</code></dt>
            <dd>Total number of input reads provided to FungiFunk.</dd>

            <dt><code>total_diamond_hit_count</code></dt>
            <dd><code>annotated_hit_count + unannotated_hit_count</code>.</dd>

            <dt><code>diamond_hits_per_input_read</code></dt>
            <dd><code>total_diamond_hit_count / total_input_reads</code>.</dd>

            <dt><code>annotated_hit_count</code></dt>
            <dd>Number of retained DIAMOND hits with EggNOG annotation.</dd>

            <dt><code>annotated_hit_fraction</code></dt>
            <dd><code>annotated_hit_count / total_diamond_hit_count</code>.</dd>

            <dt><code>unannotated_hit_count</code></dt>
            <dd>Number of retained DIAMOND hits without EggNOG annotation.</dd>

            <dt><code>unannotated_hit_fraction</code></dt>
            <dd><code>unannotated_hit_count / total_diamond_hit_count</code>.</dd>

            <dt><code>R1_input_reads</code></dt>
            <dd>Number of input reads in the R1 FASTQ file.</dd>

            <dt><code>R1_diamond_hit_count</code></dt>
            <dd>Total retained DIAMOND hits from R1 reads.</dd>

            <dt><code>R1_diamond_hits_per_input_read</code></dt>
            <dd><code>R1_diamond_hit_count / R1_input_reads</code>.</dd>

            <dt><code>R2_input_reads</code></dt>
            <dd>Number of input reads in the R2 FASTQ file.</dd>

            <dt><code>R2_diamond_hit_count</code></dt>
            <dd>Total retained DIAMOND hits from R2 reads.</dd>

            <dt><code>R2_diamond_hits_per_input_read</code></dt>
            <dd><code>R2_diamond_hit_count / R2_input_reads</code>.</dd>

            <dt><code>mean_bitscore</code></dt>
            <dd>Mean DIAMOND bitscore across annotated hits.</dd>

            <dt><code>mean_diamond_evalue</code></dt>
            <dd>Mean DIAMOND e-value across annotated hits.</dd>
        </dl>
    </div>
</details>
"""


# Build the short HTML description shown below each functional feature title
def get_feature_description_html(table_key, calculate_pfam_eggnog_context=False):
    """
    Return the short description for one functional feature type.
    """

    # Define PFAM context note depending on whether optional context was calculated
    if calculate_pfam_eggnog_context:
        pfam_context_note = """
<p class="small-text">
    The top-associated PFAM columns are dataset-derived context.
    They should not be interpreted as official PFAM database descriptions.
</p>
"""
    else:
        pfam_context_note = """
<p class="small-text">
    PFAM context columns were not calculated. Use <code>--calculate_pfam_eggnog_context</code>
    to add dataset-derived top-associated protein context.
</p>
"""

    # Define eggNOG OG context note depending on whether optional context was calculated
    if calculate_pfam_eggnog_context:
        eggnog_context_note = """
<p class="small-text">
    The top-associated eggNOG OG columns are dataset-derived context.
    They should not be interpreted as official eggNOG database descriptions.
</p>
"""
    else:
        eggnog_context_note = """
<p class="small-text">
    eggNOG OG context columns were not calculated. Use <code>--calculate_pfam_eggnog_context</code>
    to add dataset-derived top-associated protein context.
</p>
"""

    # Define feature-specific descriptions shown outside the column guide box
    descriptions = {
        "protein": """
<p class="feature-description">
    Protein features are the reference proteins reported by FungiFunk.
    Counts show how many reads were assigned to each protein across samples.
</p>
""",
        "ko": """
<p class="feature-description">
    KEGG orthologs are functional gene or protein groups.
    They provide functional context for proteins with similar molecular roles.
</p>
""",
        "cazy": """
<p class="feature-description">
    CAZy features describe carbohydrate-active enzyme families.
    They are useful for identifying enzymes involved in carbohydrate degradation, modification, or binding.
</p>
""",
        "pfam": f"""
<p class="feature-description">
    PFAM features describe conserved protein families or domains.
    A single protein can contain more than one PFAM feature, so PFAM counts represent domain or family assignments rather than unique proteins.
</p>

{pfam_context_note}
""",
        "eggnog_ogs": f"""
<p class="feature-description">
    eggNOG OGs are orthologous groups of evolutionarily related genes or proteins.
    Counts use only the most specific valid eggNOG OG reported for each annotated hit.
    The complete comma-separated eggNOG OG hierarchy remains preserved in the per-sample annotation tables.
</p>

{eggnog_context_note}
""",
    }

    # Return the matching description, or an empty string if the table type is unknown
    return descriptions.get(table_key, "")


# Build the collapsible column guide for one functional feature table
def get_feature_column_guide_html(
    table_key,
    count_view="raw",
    calculate_pfam_eggnog_context=False,
):
    """
    Return a collapsible HTML column guide for one feature type and count view.
    """

    # Select the total column name and explanation depending on the count view
    if count_view == "cpm":
        total_column = "total_CPM"
        total_description = (
            "CPM-normalised values summed across samples."
        )
        guide_title = "Column guide"
    else:
        total_column = "total_count"
        total_description = "Raw read counts summed across samples."
        guide_title = "Column guide"

    # Define PFAM context columns only when optional context was calculated
    if calculate_pfam_eggnog_context:
        pfam_context_columns = """
            <dt><code>top_associated_preferred_name</code></dt>
            <dd>Most frequent protein preferred name observed with this PFAM in the annotation tables.</dd>

            <dt><code>top_associated_description</code></dt>
            <dd>Most frequent protein description observed with this PFAM in the annotation tables.</dd>

            <dt><code>top_associated_hit_count</code></dt>
            <dd>Number of annotation-table hits supporting the top associated context.</dd>

            <dt><code>top_associated_hit_fraction</code></dt>
            <dd>Fraction of hits assigned to this eggNOG OG as their most specific OG that support the top associated context.</dd>
"""
    else:
        pfam_context_columns = ""

    # Define eggNOG OG context columns only when optional context was calculated
    if calculate_pfam_eggnog_context:
        eggnog_context_columns = """
            <dt><code>top_associated_preferred_name</code></dt>
            <dd>Most frequent protein preferred name observed with this eggNOG OG in the annotation tables.</dd>

            <dt><code>top_associated_description</code></dt>
            <dd>Most frequent protein description observed with this eggNOG OG in the annotation tables.</dd>

            <dt><code>top_associated_hit_count</code></dt>
            <dd>Number of annotation-table hits supporting the top associated context.</dd>

            <dt><code>top_associated_hit_fraction</code></dt>
            <dd>Fraction of hits containing this eggNOG OG that support the top associated context.</dd>
"""
    else:
        eggnog_context_columns = ""

    # Define feature-specific metadata columns
    metadata_columns = {
        "protein": """
            <dt><code>protein_id</code></dt>
            <dd>Matched protein identifier.</dd>

            <dt><code>species</code></dt>
            <dd>Species parsed from the protein ID.</dd>

            <dt><code>accession</code></dt>
            <dd>Protein accession parsed from the protein ID.</dd>

            <dt><code>Preferred_name</code></dt>
            <dd>Short protein name from the FungiFunk annotation output.</dd>

            <dt><code>Description</code></dt>
            <dd>Protein description from the FungiFunk annotation output.</dd>
""",
        "ko": """
            <dt><code>KEGG_ko</code></dt>
            <dd>Detected KEGG ortholog ID.</dd>

            <dt><code>representative_preferred_name</code></dt>
            <dd>Most frequent protein preferred name associated with this KO in the FungiFunk outputs.</dd>

            <dt><code>representative_description</code></dt>
            <dd>Most frequent protein description associated with this KO in the FungiFunk outputs.</dd>

            <dt><code>associated_EC</code></dt>
            <dd>EC numbers observed in proteins associated with this KO.</dd>

            <dt><code>associated_KEGG_Pathway</code></dt>
            <dd>KEGG pathways observed in proteins associated with this KO.</dd>
""",
        "cazy": """
            <dt><code>CAZy</code></dt>
            <dd>Detected CAZy family.</dd>

            <dt><code>CAZy_class</code></dt>
            <dd>Broad CAZy class extracted from the family prefix, such as GH, GT, PL, CE, AA, or CBM.</dd>

            <dt><code>CAZy_class_description</code></dt>
            <dd>Short description of the broad CAZy class.</dd>
""",
        "pfam": f"""
            <dt><code>PFAMs</code></dt>
            <dd>Detected PFAM feature.</dd>

{pfam_context_columns}
""",
        "eggnog_ogs": f"""
            <dt><code>eggNOG_OGs</code></dt>
            <dd>Most specific eggNOG orthologous group counted for annotated hits.</dd>

{eggnog_context_columns}
""",
    }

    # Get metadata columns for this feature type
    feature_columns = metadata_columns.get(table_key, "")

    # Return the collapsible column guide
    return f"""
<details class="feature-info column-guide-details">
    <summary>{guide_title}</summary>

    <div class="column-guide">
        <dl>
{feature_columns}
            <dt><code>{total_column}</code></dt>
            <dd>{total_description}</dd>
        </dl>
    </div>
</details>
"""


# Convert one detected file value into a compact True/False HTML label
def format_detected_file_status_for_html(value):
    """
    Format one detected file value for the HTML sample summary table.
    """

    # Treat real boolean values directly
    if isinstance(value, bool):

        # Return True if the boolean says the file was detected
        if value:
            return '<span class="status-true">True</span>'

        # Return False if the boolean says the file was not detected
        return '<span class="status-false">False</span>'

    # Return False for missing values
    if is_missing_value(value):
        return '<span class="status-false">False</span>'

    # Convert the value to clean text
    value_text = str(value).strip()

    # Treat common false-like text values as missing files
    if value_text.lower() in {"false", "no", "none", "nan", "na", "n/a", "-"}:
        return '<span class="status-false">False</span>'

    # Treat true-like text values as detected files without a path tooltip
    if value_text.lower() in {"true", "yes"}:
        return '<span class="status-true">True</span>'

    # Escape the path before adding it to the HTML tooltip
    escaped_path = escape(value_text, quote=True)

    # Return True with the detected path available as a tooltip
    return f'<span class="status-true" title="{escaped_path}">True</span>'


# Build a compact HTML-only version of the sample summary table
def build_sample_summary_html_display(sample_summary):
    """
    Build a compact sample summary table for the HTML report.

    The full sample_summary.tsv file keeps the complete paths.
    This HTML-only table shows True/False values and keeps paths as tooltips.
    """

    # Return an empty table if the sample summary is missing
    if sample_summary is None or sample_summary.empty:
        return pd.DataFrame()

    # Define the output columns shown in the HTML report
    html_columns = [
        (
            "Found input read counts?",
            ["input_read_counts_path", "has_input_read_counts", "input_read_counts"],
        ),
        (
            "Found annotation table?",
            ["annotation_path", "annotation_file", "annotation"],
        ),
        (
            "Found protein counts?",
            ["protein_read_counts_path", "protein_path", "protein_read_counts", "protein"],
        ),
        (
            "Found KO counts?",
            ["ko_read_counts_path", "ko_path", "ko_read_counts", "ko"],
        ),
        (
            "Found eggNOG OG counts?",
            [
                "eggnog_ogs_read_counts_path",
                "eggnog_ogs_path",
                "eggNOG_OGs_path",
                "eggnog_ogs_read_counts",
                "eggnog_ogs",
                "eggNOG_OGs",
            ],
        ),
        (
            "Found CAZy counts?",
            ["cazy_read_counts_path", "cazy_path", "CAZy_path", "cazy_read_counts", "cazy", "CAZy"],
        ),
        (
            "Found PFAM counts?",
            ["pfam_read_counts_path", "pfam_path", "PFAMs_path", "pfam_read_counts", "pfam", "PFAMs"],
        ),
        (
            "Found unannotated hits?",
            ["unannotated_hits_path", "unannotated_path", "unannotated_hits", "unannotated"],
        ),
    ]

    # Detect the sample-name column
    sample_col = "sample" if "sample" in sample_summary.columns else sample_summary.columns[0]

    # Store one compact HTML row per sample
    rows = []

    # Process each sample summary row
    for _, row in sample_summary.iterrows():

        # Store the sample name as escaped text
        display_row = {
            "sample": escape(str(row[sample_col]), quote=False)
        }

        # Convert each detected-file column to True/False HTML
        for output_label, candidate_columns in html_columns:

            # Find the first matching column available in the real sample summary table
            source_col = next(
                (col for col in candidate_columns if col in sample_summary.columns),
                None,
            )

            # Use missing if no matching source column was found
            value = row[source_col] if source_col is not None else pd.NA

            # Store the compact HTML value
            display_row[output_label] = format_detected_file_status_for_html(value)

        # Store the completed row
        rows.append(display_row)

    # Return the compact HTML display table
    return pd.DataFrame(rows)


# Check whether the compact HTML sample summary contains at least one missing file
def sample_summary_html_has_false(sample_summary_html_display):
    """
    Check whether the HTML sample summary contains at least one False value.
    """

    # Return False if the display table is missing or empty
    if sample_summary_html_display is None or sample_summary_html_display.empty:
        return False

    # Check all columns except the sample column
    for col in sample_summary_html_display.columns:

        # Skip the sample-name column
        if col == "sample":
            continue

        # Return True as soon as one missing-file status is found
        if sample_summary_html_display[col].astype(str).str.contains("status-false").any():
            return True

    # Return False if no missing-file status was found
    return False


# Read report details from the run information log
def read_report_details_from_log(run_info_log, outdir):
    """
    Read report details from the FungiFunk_summary run information log.
    """

    # Define the default report details used if the log is missing or incomplete
    details = {
        "tool": "FungiFunk_summary",
        "version": "not available",
        "run_info_log": "not available",
        "command": "not available",
        "sections": [],
    }

    # Use the standard log path if no path was provided
    if run_info_log is None:
        run_info_log = Path(outdir) / "FungiFunk_summary_run_info.log"

    # Convert the log path to a pathlib object
    run_info_log = Path(run_info_log)

    # Store the log path even if the file cannot be read
    details["run_info_log"] = str(run_info_log)

    # Return default details if the log file does not exist
    if not run_info_log.exists():
        details["sections"] = [
            {
                "title": "",
                "rows": [
                    (
                        "Run information log",
                        f"Run information log could not be found: {run_info_log}",
                    )
                ],
            }
        ]
        return details

    # Store rows that appear before named sections such as Input or Run options
    current_section = {
        "title": "",
        "rows": [],
    }

    # Store all report detail sections
    sections = [current_section]

    # Read the run information log
    try:
        with open(run_info_log, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:

                # Remove the line break from the current line
                line = line.rstrip("\n")

                # Skip empty lines because spacing is handled by the HTML layout
                if not line.strip():
                    continue

                # Remove leading and trailing spaces for parsing
                stripped_line = line.strip()

                # Skip the tool name because it is already shown in the report introduction
                if stripped_line.startswith("Tool:"):
                    details["tool"] = stripped_line.split(":", 1)[1].strip()
                    continue

                # Reformat the version line for clearer display
                if stripped_line.startswith("Version:"):
                    version = stripped_line.split(":", 1)[1].strip()
                    details["version"] = version
                    current_section["rows"].append(
                        ("FungiFunk_summary version", version)
                    )
                    continue

                # Store the command separately so it can be shown in its own code block
                if stripped_line.startswith("Command:"):
                    details["command"] = stripped_line.split(":", 1)[1].strip()
                    continue

                # Start a new section when a section header is found
                if stripped_line.endswith(":") and ":" not in stripped_line[:-1]:
                    current_section = {
                        "title": stripped_line,
                        "rows": [],
                    }
                    sections.append(current_section)
                    continue

                # Store regular key-value rows
                if ":" in stripped_line:
                    key, value = stripped_line.split(":", 1)
                    current_section["rows"].append(
                        (key.strip(), value.strip())
                    )
                    continue

                # Keep unusual lines if they ever appear in the log
                current_section["rows"].append((stripped_line, ""))

    # Keep default details if the log cannot be read
    except Exception as exc:
        details["sections"] = [
            {
                "title": "",
                "rows": [
                    ("Run information log", f"Could not be read: {run_info_log}"),
                    ("Reason", str(exc)),
                ],
            }
        ]
        return details

    # Add the run information log path as an extra useful report detail
    sections[0]["rows"].append(("Run information log", str(run_info_log)))

    # Remove empty sections
    details["sections"] = [
        section for section in sections if section["rows"] or section["title"]
    ]

    # Return the collected report details
    return details


# Build the HTML details block shown near the top of the report
def build_report_details_html(report_details):
    """
    Build a collapsible HTML block with report metadata.
    """

    # Store rendered HTML for each report detail section
    section_html_blocks = []

    # Build one formatted block per report detail section
    for section in report_details.get("sections", []):

        # Escape the section title before inserting it into HTML
        section_title = escape(str(section.get("title", "")), quote=True)

        # Store all rows for this section
        row_html_blocks = []

        # Build one two-column row per metadata entry
        for key, value in section.get("rows", []):

            # Escape key and value before inserting them into HTML
            key = escape(str(key), quote=True)
            value = escape(str(value), quote=True)

            # Add the formatted metadata row
            row_html_blocks.append(
                f"""
            <dt>{key}:</dt>
            <dd>{value}</dd>
"""
            )

        # Join all rows for this section
        rows_html = "\n".join(row_html_blocks)

        # Add a section title only for named sections
        if section_title:
            section_heading_html = f"<h4>{section_title}</h4>"
        else:
            section_heading_html = ""

        # Add the section block
        section_html_blocks.append(
            f"""
        <div class="report-details-section">
            {section_heading_html}
            <dl>
{rows_html}
            </dl>
        </div>
"""
        )

    # Join all rendered report detail sections
    sections_html = "\n".join(section_html_blocks)

    # Escape the command before inserting it into the HTML code block
    command = escape(
        str(report_details.get("command", "not available")),
        quote=True,
    )

    # Return the collapsible report details block
    return f"""
<details class="report-details">
    <summary>Report details</summary>

    <div class="report-details-content">
{sections_html}

        <p><strong>Command used:</strong></p>
        <pre class="command-block"><code>{command}</code></pre>
    </div>
</details>
"""


# Generate the interactive HTML report from summary tables and functional matrices
def generate_html_report(
    outdir,
    matrices_dir,
    summaries_dir,
    sample_summary,
    annotation_summary,
    detected_features_summary,
    raw_matrices,
    cpm_matrices,
    feature_metadata=None,
    html_top_n=50,
    html_no_top=False,
    normalise_cpm=False,
    calculate_pfam_eggnog_context=False,
    time_stamps=False,
    run_info_log=None,
):
    """
    Generate the FungiFunk_summary HTML report.
    """
    # Import HTML-report dependencies only when the report is generated
    try:
        from jinja2 import Template
        import plotly.offline as plotly_offline
    # If dependencies are missing, skip the HTML report but keep the rest of the workflow running
    except ImportError:
        warn_msg(
            "HTML report could not be generated because jinja2 and/or plotly are not installed. "
            "Install them or run with --skip_html_report."
        )
        return

    # Get Plotly JavaScript so interactive plots work inside the standalone HTML report
    plotly_js = plotly_offline.get_plotlyjs()

    # Build a compact HTML-only sample summary table
    sample_summary_html_display = build_sample_summary_html_display(sample_summary)

    # Check whether any expected FungiFunk output file was not detected
    sample_summary_has_false = sample_summary_html_has_false(sample_summary_html_display)

    # Define the sample summary TSV link for the HTML report
    sample_summary_link = relative_link(summaries_dir / "sample_summary.tsv", outdir)

    # Define the annotation summary TSV link for the HTML report
    annotation_summary_link = relative_link(summaries_dir / "annotation_summary.tsv", outdir)

    # Define the detected features summary TSV link for the HTML report
    detected_features_summary_link = relative_link(summaries_dir / "detected_features_summary.tsv", outdir)

    # Read report metadata from the run information log
    report_details = read_report_details_from_log(run_info_log=run_info_log, outdir=outdir)

    # Build the collapsible report details block
    report_details_html = build_report_details_html(report_details)

    # Build the functional feature boxplot explanation used in the HTML report
    functional_boxplot_note_html = get_functional_boxplot_note_html()

    # Build the functional feature clustered heatmap explanation used in the HTML report
    functional_clustermap_note_html = get_functional_clustermap_note_html()

    # Get tool name and version for the report introduction
    report_tool = report_details.get("tool", "FungiFunk_summary")
    report_version = report_details.get("version", "not available")

    # Create a warning note only when at least one expected file was not detected
    sample_summary_footer_html = ""

    if sample_summary_has_false:
        sample_summary_footer_html = """
    <div class="warning sample-summary-warning">
        <p>
            <strong>Warning:</strong> At least one expected FungiFunk output file was not detected.
            FungiFunk_summary searches recursively inside the input directory for files with the expected suffixes:
        </p>

        <ul>
            <li><code>*_input_read_counts.tsv</code></li>
            <li><code>*_annotation.tsv</code></li>
            <li><code>*_protein_read_counts.tsv</code></li>
            <li><code>*_ko_read_counts.tsv</code></li>
            <li><code>*_eggnog_ogs_read_counts.tsv</code></li>
            <li><code>*_cazy_read_counts.tsv</code></li>
            <li><code>*_pfam_read_counts.tsv</code></li>
            <li><code>*_unannotated_hits.tsv</code></li>
        </ul>

        <p>
            A <strong>False</strong> value means that no matching file was detected for that sample and output type.
            This may happen if the file is missing, renamed, stored outside the input directory,
            or does not follow the expected filename suffix.
        </p>
    </div>
    """

    # Build annotation summary diagnostic plot sections
    annotation_diagnostic_sections = build_annotation_diagnostic_sections(annotation_summary)

    # Build detected features summary diagnostic plot sections
    detected_features_diagnostic_sections = build_detected_features_diagnostic_sections(detected_features_summary)

    summary_tables = [
        {
            "title": "Sample summary",
            "description": (
                "Quick overview of detected FungiFunk output files per sample. "
                "Hover over <strong>True</strong> to see the detected file path. <br> "
                f"Full paths are also available in <a href='{sample_summary_link}'>summaries/sample_summary.tsv</a>."
            ),
            "note_html": "",
            "html": dataframe_to_html(sample_summary_html_display),
            "footer_html": sample_summary_footer_html,
        },
        {
            "title": "Annotation summary",
            "description": (
                "This table summarises the per-sample progression from input reads, "
                "to retained DIAMOND hits, to functionally annotated and unannotated hits. "
                "Fraction columns are reported as values from 0 to 1, not as percentages. <br>"
                f"Full table is also available in <a href='{annotation_summary_link}'>summaries/annotation_summary.tsv</a>."
            ),
            "note_html": get_annotation_summary_note_html(),
            "html": dataframe_to_html(annotation_summary),
            "footer_html": "",
        },
        {
            "title": "Detected features summary",
            "description": (
                "Number of detected functional features per sample. <br>"
                f"Full table is also available in <a href='{detected_features_summary_link}'>summaries/detected_features_summary.tsv</a>."
            ),
            "note_html": "",
            "html": dataframe_to_html(detected_features_summary),
            "footer_html": "",
        },
    ]

    # Prepare functional feature sections for the HTML report
    feature_sections = []

    # Process each functional matrix type
    for table_key, matrix in raw_matrices.items():

        # Skip missing or empty matrices
        if matrix is None or matrix.empty:
            continue

        # Get information for this count table
        info = COUNT_TABLES[table_key]

        # Get a readable title for this feature type
        section_title = info.get("title", table_key)

        # Get the feature ID column for this table type
        id_col = info["id_col"]

        # Select top features from the raw count matrix
        selected_matrix = select_top_count_matrix(
            matrix,
            top_n=html_top_n,
            no_top=html_no_top,
        )

        # Create the enriched raw top-feature table
        raw_top_table = matrix_top_table(
            matrix=matrix,
            top_matrix=selected_matrix,
            table_key=table_key,
            feature_metadata=feature_metadata,
            total_col="total_count",
        )

        # Define the raw matrix path
        raw_matrix_path = matrices_dir / info["raw_filename"]

        # Define the CPM matrix path
        cpm_matrix_path = matrices_dir / info["cpm_filename"]

        # Create the CPM matrix link only when CPM normalisation was requested and the file exists
        cpm_matrix_link = (
            relative_link(cpm_matrix_path, outdir)
            if normalise_cpm and cpm_matrix_path.exists()
            else None
        )

        # Prepare CPM-specific HTML elements
        cpm_top_table_html = None
        cpm_boxplot_title = None
        cpm_boxplot_html = None
        cpm_clustermap_title = None
        cpm_clustermap_html = None

        # Add CPM tables and plots only when CPM normalisation was requested and matrices are available
        if normalise_cpm and cpm_matrices is not None and table_key in cpm_matrices:

            # Get the CPM matrix for this feature type
            cpm_matrix = cpm_matrices[table_key]

            # Continue only if the CPM matrix is available and not empty
            if cpm_matrix is not None and not cpm_matrix.empty:

                # Select top features from the CPM matrix
                selected_cpm_matrix = select_top_count_matrix(
                    cpm_matrix,
                    top_n=html_top_n,
                    no_top=html_no_top,
                )

                # Create the enriched CPM top-feature table
                cpm_top_table = matrix_top_table(
                    matrix=cpm_matrix,
                    top_matrix=selected_cpm_matrix,
                    table_key=table_key,
                    feature_metadata=feature_metadata,
                    total_col="total_CPM",
                )

                # Convert the CPM top-feature table to HTML
                cpm_top_table_html = dataframe_to_html(cpm_top_table)

                # Define CPM plot titles shown in the HTML report
                cpm_boxplot_title = f"{section_title}: CPM distribution per feature"
                cpm_clustermap_title = f"{section_title}: clustered CPM abundance patterns"

                # Create the CPM boxplot without an internal Plotly title
                cpm_boxplot_html = make_matrix_boxplot_html(
                    selected_cpm_matrix,
                    title="",
                    feature_col=id_col,
                    value_label="CPM",
                    feature_label="Feature",
                )

                # Create the CPM clustered heatmap without an internal Plotly title
                cpm_clustermap_html = make_clustermap_html(
                    selected_cpm_matrix,
                    title="",
                    colorbar_title="CPM",
                    row_label="Feature",
                )

        # Define raw-count plot titles shown in the HTML report
        raw_boxplot_title = f"{section_title}: raw read count distribution per feature"
        raw_clustermap_title = f"{section_title}: clustered raw read count patterns"

        # Create the raw-count boxplot without an internal Plotly title
        raw_boxplot_html = make_matrix_boxplot_html(
            selected_matrix,
            title="",
            feature_col=id_col,
            value_label="Raw read count",
            feature_label="Feature",
        )

        # Create the raw-count clustered heatmap without an internal Plotly title
        raw_clustermap_html = make_clustermap_html(
            selected_matrix,
            title="",
            colorbar_title="read count",
            row_label="Feature",
        )

        # Store all HTML elements for this functional feature section
        section = {
            "title": section_title,
            "raw_matrix_link": relative_link(raw_matrix_path, outdir),
            "cpm_matrix_link": cpm_matrix_link,
            "description_html": get_feature_description_html(
                table_key,
                calculate_pfam_eggnog_context=calculate_pfam_eggnog_context,
            ),
            "raw_column_guide_html": get_feature_column_guide_html(
                table_key=table_key,
                count_view="raw",
                calculate_pfam_eggnog_context=calculate_pfam_eggnog_context,
            ),
            "cpm_column_guide_html": get_feature_column_guide_html(
                table_key=table_key,
                count_view="cpm",
                calculate_pfam_eggnog_context=calculate_pfam_eggnog_context,
            ),
            "raw_top_table_html": dataframe_to_html(raw_top_table),
            "raw_boxplot_title": raw_boxplot_title,
            "raw_boxplot_html": raw_boxplot_html,
            "raw_clustermap_title": raw_clustermap_title,
            "raw_clustermap_html": raw_clustermap_html,
            "cpm_top_table_html": cpm_top_table_html,
            "cpm_boxplot_title": cpm_boxplot_title,
            "cpm_boxplot_html": cpm_boxplot_html,
            "cpm_clustermap_title": cpm_clustermap_title,
            "cpm_clustermap_html": cpm_clustermap_html,
        }

        # Store this functional feature section
        feature_sections.append(section)

    # Define the HTML report template
    html_template = Template(
        """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <title>FungiFunk summary report</title>
    <script type="text/javascript">{{ plotly_js | safe }}</script>
    <style>
        body {
            font-family: Arial, sans-serif;
            margin: 30px;
            color: #222;
            line-height: 1.5;
        }

        h1 {
            border-bottom: 2px solid #333;
            padding-bottom: 8px;
        }

        h2 {
            margin-top: 40px;
            border-bottom: 1px solid #aaa;
            padding-bottom: 5px;
        }

        h3 {
            margin-top: 30px;
        }

        .note {
            background: #f5f5f5;
            border-left: 4px solid #777;
            padding: 10px;
            margin: 15px 0;
        }

        .warning {
            background: #fff3cd;
            border-left: 4px solid #cc9900;
            padding: 10px;
            margin: 15px 0;
        }

        .summary-table {
            border-collapse: collapse;
            font-size: 13px;
            margin: 15px 0;
            width: 100%;
        }

        .summary-table th {
            background: #eeeeee;
            border: 1px solid #cccccc;
            padding: 6px;
            text-align: left;
        }

        .summary-table td {
            border: 1px solid #dddddd;
            padding: 6px;
            vertical-align: top;
        }

        .scroll-box {
            overflow-x: auto;
            max-width: 100%;
        }

        .small-text {
            font-size: 12px;
            color: #555;
        }

        .status-true {
            font-weight: bold;
            color: #1f7a1f;
            cursor: help;
        }

        .status-false {
            font-weight: bold;
            color: #a33a3a;
        }

        .sample-summary-warning {
            margin-top: 10px;
        }

        .sample-summary-warning p {
            margin: 8px 0;
        }

        .sample-summary-warning ul {
            margin: 8px 0 8px 24px;
        }

        .sample-summary-warning li {
            margin: 3px 0;
        }

        .feature-info {
            background: #f8f9fb;
            border: 1px solid #dddddd;
            border-left: 4px solid #4f6f8f;
            padding: 12px 14px;
            margin: 12px 0 18px 0;
        }

        .feature-description {
            margin: 0 0 12px 0;
        }

        .column-guide {
            margin-top: 8px;
        }

        .column-guide h4 {
            margin: 8px 0;
        }

        .column-guide dl {
            display: grid;
            grid-template-columns: 260px 1fr;
            column-gap: 16px;
            row-gap: 6px;
            margin: 0;
        }

        .column-guide dt {
            font-weight: bold;
        }

        .column-guide dd {
            margin: 0;
        }

        .column-guide-details summary {
            font-weight: bold;
            cursor: pointer;
            margin: 0;
        }

        .column-guide-details[open] summary {
            margin-bottom: 10px;
        }

        .column-guide-details .column-guide {
            margin-top: 8px;
        }

        .feature-count-view-details {
            margin: 6px 0 8px 32px;
            padding-bottom: 0;
        }

        .feature-count-view-details > summary {
            display: block;
            font-size: 16px;
            font-weight: bold;
            cursor: pointer;
            margin: 4px 0 6px 0;
        }

        .feature-count-view-details > summary::before {
            content: "▶ ";
            font-size: 13px;
        }

        .feature-count-view-details[open] > summary::before {
            content: "▼ ";
        }

        .feature-count-view-details > summary::-webkit-details-marker {
            display: none;
        }

        .feature-count-view-content {
            margin-top: 6px;
        }

        .analysis-note {
            background: #f5f5f5;
            border-left: 4px solid #aaaaaa;
            padding: 8px 10px;
            margin: 8px 0 15px 0;
            font-size: 13px;
        }

        .pagination-controls {
            display: flex;
            align-items: center;
            gap: 8px;
            margin: 8px 0 20px 0;
            font-size: 13px;
        }

        .pagination-controls button {
            border: 1px solid #cccccc;
            background: #f5f5f5;
            padding: 4px 8px;
            cursor: pointer;
        }

        .pagination-controls button:disabled {
            opacity: 0.5;
            cursor: not-allowed;
        }

        .pagination-info {
            color: #555555;
        }

        a {
            color: #005ea8;
        }

        .report-toc {
            margin: 18px 0 28px 0;
        }

        .report-toc > summary {
            display: list-item;
            font-size: 1.5em;
            font-weight: bold;
            cursor: pointer;
            margin-top: 40px;
            border-bottom: 1px solid #aaa;
            padding-bottom: 5px;
        }

        .report-toc-content {
            margin-top: 12px;
        }

        .report-toc a {
            color: inherit;
            text-decoration: none;
        }

        .report-toc a:hover {
            text-decoration: underline;
        }

        .report-toc ol {
            list-style-type: none;
            margin-top: 6px;
            margin-bottom: 6px;
            padding-left: 22px;
        }

        .report-toc ol ol {
            padding-left: 24px;
        }

        .report-toc li {
            margin: 4px 0;
        }

        .report-toc .toc-views {
            margin-top: 4px;
            margin-bottom: 8px;
            padding-left: 22px;
            list-style-type: none;
        }

        .report-toc .toc-views li {
            margin: 3px 0;
        }

        .toc-help {
            margin-top: 0;
        }

        .summary-table-details {
            margin: 22px 0 34px 0;
            padding-bottom: 8px;
        }

        .summary-table-details > summary {
            display: block;
            font-size: 18px;
            font-weight: bold;
            cursor: pointer;
            margin: 18px 0 14px 0;
        }

        .summary-table-details > summary::before {
            content: "▼ ";
            font-size: 14px;
        }

        .summary-table-details:not([open]) > summary::before {
            content: "▶ ";
        }

        .summary-table-details > summary::-webkit-details-marker {
            display: none;
        }

        .summary-table-content {
            margin-top: 10px;
        }

        .report-details {
            background: #f8f9fb;
            border: 1px solid #dddddd;
            border-left: 4px solid #4f6f8f;
            padding: 12px 16px;
            margin: 16px 0 24px 0;
        }

        .report-details > summary {
            display: list-item;
            font-weight: bold;
            cursor: pointer;
            font-size: 15px;
            margin: 0;
        }

        .report-details-content {
            margin-top: 12px;
            font-size: 13px;
        }

        .report-details-section {
            margin-bottom: 12px;
        }

        .report-details-section h4 {
            margin: 12px 0 6px 0;
            font-size: 13px;
        }

        .report-details dl {
            display: grid;
            grid-template-columns: 220px 1fr;
            column-gap: 14px;
            row-gap: 5px;
            margin: 0;
        }

        .report-details dt {
            font-weight: bold;
        }

        .report-details dd {
            margin: 0;
        }

        .report-details-section h4 + dl {
            margin-left: 22px;
        }

        .command-block {
            white-space: pre-wrap;
            overflow-x: auto;
            background: #eeeeee;
            padding: 10px;
            border-radius: 4px;
        }

        .feature-section-details {
            margin: 24px 0 34px 0;
            padding-bottom: 8px;
        }

        .feature-section-details > summary {
            display: block;
            font-size: 18px;
            font-weight: bold;
            cursor: pointer;
            margin: 18px 0 14px 0;
        }

        .feature-section-details > summary::before {
            content: "▼ ";
            font-size: 14px;
        }

        .feature-section-details:not([open]) > summary::before {
            content: "▶ ";
        }

        .feature-section-details > summary::-webkit-details-marker {
            display: none;
        }

        .feature-section-content {
            margin-top: 10px;
        }

        .intro-divider {
            border: 0;
            border-top: 1px solid #aaaaaa;
            margin: 24px 0 28px 0;
        }

        .summary-table-details > summary,
        .feature-section-details > summary,
        .feature-count-view-details > summary {
            list-style: none;
        }

        .summary-table-details > summary::marker,
        .feature-section-details > summary::marker,
        .feature-count-view-details > summary::marker {
            content: "";
            display: none;
        }

        .plot-block {
            margin-top: 2.6rem;
            margin-bottom: 3.4rem;
            padding-top: 0;
        }

        .plot-block + .plot-block {
            border-top: 2px solid #d5dbe3;
            margin-top: 4rem;
            padding-top: 2.4rem;
        }

        .plot-title {
            font-family: "Open Sans", verdana, arial, sans-serif;
            font-size: 17px;
            font-weight: 400;
            color: #2a3f5f;
            text-align: left;
            margin: 0 0 0.65rem 80px;
            line-height: 1.3;
        }

        .plot-block .small-text {
            margin-top: 0;
            margin-bottom: 1.2rem;
            margin-left: 80px;
            max-width: 900px;
        }

        @media (max-width: 800px) {
            .plot-title,
            .plot-block .small-text {
                margin-left: 0;
            }
        }

    </style>
</head>

<body>

    {% macro plot_block(title, note_html, plot_html) %}
        {% if plot_html %}
            <div class="plot-block">
                {% if title %}
                    <div class="plot-title">{{ title }}</div>
                {% endif %}

                {% if note_html %}
                    {{ note_html | safe }}
                {% endif %}

                {{ plot_html | safe }}
            </div>
        {% endif %}
    {% endmacro %}

    <h1>FungiFunk summary report</h1>

    <p>
        This report was generated by <strong>{{ report_tool }}{% if report_version and report_version != "not available" %} v{{ report_version }}{% endif %}</strong>. <br>
        It summarises multiple FungiFunk runs by showing sample-level output checks,
        annotation summaries, and detected feature summaries. It also compares detected
        functional features across samples, including proteins, KEGG orthologs, eggNOG
        orthologous groups, CAZy families, and PFAM domains. <br>
        Complete TSV matrices and summary tables are written to the output directory. <br>
        To keep the HTML report responsive, functional tables and plots show the top 50
        features by default. This limit can be changed with <code>--html_top_n</code>
        or disabled with <code>--html_no_top</code>.
    </p>

    {% if report_details_html %}
        {{ report_details_html | safe }}
    {% endif %}

    {% if html_no_top %}
    <div class="warning">
        <strong>Warning:</strong> Functional feature top filtering was disabled with <code>--html_no_top</code>.
        This HTML report may be large and slow to open.
    </div>
    {% endif %}

    <details class="report-toc" open>
        <summary>Contents</summary>

        <nav class="report-toc-content">
            <p class="toc-help">
                Click a section name to jump directly to that part of the report.
            </p>

            <ol>
                <li>
                    <a href="#summary-tables">1. Summary tables</a>
                    <ol>
                        {% for table in summary_tables %}
                            <li>
                                <a href="#summary-table-{{ loop.index }}">
                                    1.{{ loop.index }}. {{ table.title }}
                                </a>

                                {% if table.title == "Annotation summary" and annotation_diagnostic_sections %}
                                    <ul class="toc-views">
                                        {% for diagnostic in annotation_diagnostic_sections %}
                                            <li>
                                                <a href="#{{ diagnostic.id }}">
                                                    {{ diagnostic.title }}
                                                </a>
                                            </li>
                                        {% endfor %}
                                    </ul>
                                {% endif %}

                                {% if table.title == "Detected features summary" and detected_features_diagnostic_sections %}
                                    <ul class="toc-views">
                                        {% for diagnostic in detected_features_diagnostic_sections %}
                                            <li>
                                                <a href="#{{ diagnostic.id }}">
                                                    {{ diagnostic.title }}
                                                </a>
                                            </li>
                                        {% endfor %}
                                    </ul>
                                {% endif %}
                            </li>
                        {% endfor %}
                    </ol>
                </li>

                <li>
                    <a href="#functional-features">2. Functional features</a>
                    <ol>
                        {% for section in feature_sections %}
                            <li>
                                <a href="#functional-feature-{{ loop.index }}">
                                    2.{{ loop.index }}. {{ section.title }}
                                </a>

                                <ul class="toc-views">
                                    <li>
                                        <a href="#functional-feature-{{ loop.index }}-raw">
                                            Raw read counts
                                        </a>
                                    </li>

                                    {% if section.cpm_top_table_html %}
                                    <li>
                                        <a href="#functional-feature-{{ loop.index }}-cpm">
                                            CPM-normalised counts
                                        </a>
                                    </li>
                                    {% endif %}
                                </ul>
                            </li>
                        {% endfor %}
                    </ol>
                </li>
            </ol>
        </nav>
    </details>

    <h2 id="summary-tables">1. Summary tables</h2>

    {% for table in summary_tables %}
        <details id="summary-table-{{ loop.index }}" class="summary-table-details" open>
            <summary>1.{{ loop.index }}. {{ table.title }}</summary>

            <div class="summary-table-content">
                {% if table.description %}
                    <p>{{ table.description | safe }}</p>
                {% endif %}

                {% if table.note_html %}
                    {{ table.note_html | safe }}
                {% endif %}

                <div class="scroll-box">
                    {{ table.html | safe }}
                </div>

                {% if table.title == "Annotation summary" and annotation_diagnostic_sections %}
                    {% for diagnostic in annotation_diagnostic_sections %}
                        <details id="{{ diagnostic.id }}" class="annotation-diagnostic-details feature-count-view-details">
                            <summary>{{ diagnostic.title }}</summary>

                            <div class="feature-count-view-content">
                                <p>
                                    {{ diagnostic.description }}
                                </p>

                                {{ plot_block(diagnostic.plot_title, "", diagnostic.plot_html) }}
                            </div>
                        </details>
                    {% endfor %}
                {% endif %}

                {% if table.title == "Detected features summary" and detected_features_diagnostic_sections %}
                    {% for diagnostic in detected_features_diagnostic_sections %}
                        <details id="{{ diagnostic.id }}" class="detected-features-diagnostic-details feature-count-view-details">
                            <summary>{{ diagnostic.title }}</summary>

                            <div class="feature-count-view-content">
                                <p>
                                    {{ diagnostic.description }}
                                </p>

                                {{ plot_block(diagnostic.plot_title, "", diagnostic.plot_html) }}
                            </div>
                        </details>
                    {% endfor %}
                {% endif %}

            {% if table.footer_html %}
                {{ table.footer_html | safe }}
            {% endif %}
            </div>
        </details>
    {% endfor %}

    <h2 id="functional-features">2. Functional features</h2>

    <p>
        Functional features include matched proteins, KEGG orthologs, eggNOG orthologous groups,
        CAZy families, and PFAM domains. This section compares their abundance across samples
        using FungiFunk read-count tables. For each feature type, the report shows selected
        features as tables, boxplots, and clustered heatmaps based on raw read counts. <br>
        By default, the HTML report displays the top 50 features by total abundance across samples, unless
        <code>--html_top_n</code> or <code>--html_no_top</code> was used. <br>
        If <code>--normalise_cpm</code> was used, CPM-normalised views are also included and selected
        independently by total CPM across samples.<br>
        Feature metadata shown in the functional tables is generated from FungiFunk outputs.
    </p>

    {% for section in feature_sections %}
        <details id="functional-feature-{{ loop.index }}" class="feature-section-details" open>
            <summary>2.{{ loop.index }}. {{ section.title }}</summary>

            <div class="feature-section-content">
                {% if section.description_html %}
                    {{ section.description_html | safe }}
                {% endif %}

                <details id="functional-feature-{{ loop.index }}-raw" class="feature-count-view-details">
                    <summary>Raw read counts</summary>

                    <div class="feature-count-view-content">
                        <p>
                            Features are selected by total raw count summed across samples according to
                            <code>--html_top_n</code> (default: 50), unless <code>--html_no_top</code> was used.
                        </p>

                        {% if section.raw_column_guide_html %}
                            {{ section.raw_column_guide_html | safe }}
                        {% endif %}

                        <p>
                            Full raw matrix:
                            <a href="{{ section.raw_matrix_link }}">{{ section.raw_matrix_link }}</a>
                        </p>

                        <div class="scroll-box">
                            {{ section.raw_top_table_html | safe }}
                        </div>

                        {{ plot_block(section.raw_boxplot_title, functional_boxplot_note_html, section.raw_boxplot_html) }}

                        {{ plot_block(section.raw_clustermap_title, functional_clustermap_note_html, section.raw_clustermap_html) }}
                    </div>
                </details>

                {% if section.cpm_top_table_html %}
                <details id="functional-feature-{{ loop.index }}-cpm" class="feature-count-view-details">
                    <summary>CPM-normalised counts</summary>

                    <div class="feature-count-view-content">
                        <p>
                            Features are selected independently by total CPM summed across samples according to
                            <code>--html_top_n</code> (default: 50), unless <code>--html_no_top</code> was used.
                        </p>

                        {% if section.cpm_column_guide_html %}
                            {{ section.cpm_column_guide_html | safe }}
                        {% endif %}

                        {% if section.cpm_matrix_link %}
                        <p>
                            Full CPM matrix:
                            <a href="{{ section.cpm_matrix_link }}">{{ section.cpm_matrix_link }}</a>
                        </p>
                        {% endif %}

                        <div class="scroll-box">
                            {{ section.cpm_top_table_html | safe }}
                        </div>

                        {{ plot_block(section.cpm_boxplot_title, functional_boxplot_note_html, section.cpm_boxplot_html) }}

                        {{ plot_block(section.cpm_clustermap_title, functional_clustermap_note_html, section.cpm_clustermap_html) }}
                    </div>
                </details>
                {% endif %}
            </div>
        </details>
    {% endfor %}

    <script type="text/javascript">
        document.addEventListener("DOMContentLoaded", function () {
            const rowsPerPage = 10;

            document.querySelectorAll("table.paginated-table").forEach(function (table) {
                const tbody = table.querySelector("tbody");

                if (!tbody) {
                    return;
                }

                const rows = Array.from(tbody.querySelectorAll("tr"));

                if (rows.length <= rowsPerPage) {
                    return;
                }

                let currentPage = 1;
                const totalPages = Math.ceil(rows.length / rowsPerPage);

                const controls = document.createElement("div");
                controls.className = "pagination-controls";

                const previousButton = document.createElement("button");
                previousButton.textContent = "Previous";

                const nextButton = document.createElement("button");
                nextButton.textContent = "Next";

                const pageInfo = document.createElement("span");
                pageInfo.className = "pagination-info";

                controls.appendChild(previousButton);
                controls.appendChild(pageInfo);
                controls.appendChild(nextButton);

                table.insertAdjacentElement("afterend", controls);

                function showPage(page) {
                    currentPage = page;

                    const start = (currentPage - 1) * rowsPerPage;
                    const end = start + rowsPerPage;

                    rows.forEach(function (row, index) {
                        if (index >= start && index < end) {
                            row.style.display = "";
                        } else {
                            row.style.display = "none";
                        }
                    });

                    pageInfo.textContent = "Page " + currentPage + " of " + totalPages + " (" + rows.length + " rows)";
                    previousButton.disabled = currentPage === 1;
                    nextButton.disabled = currentPage === totalPages;
                }

                previousButton.addEventListener("click", function () {
                    if (currentPage > 1) {
                        showPage(currentPage - 1);
                    }
                });

                nextButton.addEventListener("click", function () {
                    if (currentPage < totalPages) {
                        showPage(currentPage + 1);
                    }
                });

                showPage(1);
            });
        });
    </script>

    <script>
        // Open collapsed parent sections when clicking internal report links
        document.querySelectorAll('a[href^="#"]').forEach(function(linkElement) {
            linkElement.addEventListener("click", function(event) {
                var targetId = linkElement.getAttribute("href").substring(1);
                var targetElement = document.getElementById(targetId);

                if (!targetElement) {
                    return;
                }

                event.preventDefault();

                var parentDetails = targetElement.closest("details");

                while (parentDetails) {
                    parentDetails.open = true;
                    parentDetails = parentDetails.parentElement.closest("details");
                }

                targetElement.scrollIntoView({
                    behavior: "smooth",
                    block: "start"
                });

                if (window.history && window.history.pushState) {
                    window.history.pushState(null, "", "#" + targetId);
                }

                if (window.Plotly) {
                    targetElement.querySelectorAll(".plotly-graph-div").forEach(function(plotElement) {
                        Plotly.Plots.resize(plotElement);
                    });
                }
            });
        });

        // Resize Plotly figures when a collapsed report section is opened
        document.querySelectorAll(
            ".report-toc, .report-details, .summary-table-details, .feature-section-details, .feature-count-view-details"
        ).forEach(function(detailsElement) {
            detailsElement.addEventListener("toggle", function() {
                if (detailsElement.open && window.Plotly) {
                    detailsElement.querySelectorAll(".plotly-graph-div").forEach(function(plotElement) {
                        Plotly.Plots.resize(plotElement);
                    });
                }
            });
        });
    </script>

</body>
</html>
        """
    )

    # Fill the HTML template with generated tables, plots, links, and warning flags
    html = html_template.render(
        plotly_js=plotly_js,
        summary_tables=summary_tables,
        annotation_diagnostic_sections=annotation_diagnostic_sections,
        detected_features_diagnostic_sections=detected_features_diagnostic_sections,
        feature_sections=feature_sections,
        functional_boxplot_note_html=functional_boxplot_note_html,
        functional_clustermap_note_html=functional_clustermap_note_html,
        html_no_top=html_no_top,
        report_tool=report_tool,
        report_version=report_version,
        report_details_html=report_details_html,
    )

    # Define the final HTML report output path
    html_path = outdir / "FungiFunk_summary_report.html"

    # Write the completed HTML report to disk
    with open(html_path, "w", encoding="utf-8") as handle:
        handle.write(html)
    
    # Tell the user where the html was writen
    log_msg(f"    (Written to: {html_path})", time_stamps)



# Main function that runs the full FungiFunk_summary workflow
def main():
    """
    Run FungiFunk_summary.
    """
    # Read command-line arguments passed by the Bash wrapper
    args = parse_args()

    # Convert input and output directories to absolute paths
    input_dir = Path(args.input_dir).resolve()
    outdir = Path(args.outdir).resolve()

    # Define standard subdirectories inside the summary output directory
    matrices_dir = outdir / "matrices"
    summaries_dir = outdir / "summaries"

    # Stop if the input directory does not exist
    if not input_dir.exists():
        error_exit(f"--input_dir does not exist: {input_dir}")

    # Create output directories if they do not already exist
    outdir.mkdir(parents=True, exist_ok=True)
    matrices_dir.mkdir(parents=True, exist_ok=True)
    summaries_dir.mkdir(parents=True, exist_ok=True)

    # Warn the user if all functional features will be included in the HTML report
    if args.html_no_top:
        warn_msg(
            "--html_no_top was selected. The HTML report may become very large and slow to open."
        )


    print(f"################################ PROCESSING FungiFunk OUTPUT DATA ################################", flush=True)

##### Checking FungiFunk output data #####
    # Search for FungiFunk output files inside the input directory
    log_msg("Checking FungiFunk outputs...", args.timestamps)

    samples, inventory_df, duplicate_records = discover_fungifunk_outputs(
        input_dir,
        strict=args.strict,
    )


    # Check whether the detected outputs are complete enough to continue
    validate_discovered_outputs(samples, strict=args.strict)


    # Write an inventory table with all detected FungiFunk output files
    log_msg("Writing detected output files table...", args.timestamps)

    inventory_path = summaries_dir / "detected_output_files.tsv"

    inventory_df.to_csv(inventory_path, sep="\t", index=False)

    log_msg(f"    (Written to: {inventory_path})", args.timestamps)


    # Report how many samples were detected
    log_msg(f"Number of samples detected in output files: {len(samples)}.", args.timestamps)


    # Build and write the per-sample file detection summary
    log_msg("Writing sample summary...", args.timestamps)

    sample_summary = build_sample_summary(samples)

    sample_summary_path = summaries_dir / "sample_summary.tsv"

    sample_summary.to_csv(
        sample_summary_path,
        sep="\t",
        index=False,
    )

    log_msg(f"    (Written to: {sample_summary_path})", args.timestamps)


##### Count matrices #####
    # Build and write raw count matrices and optional CPM matrices
    log_msg("Writing matrices based on Cazy, EggnogOG, KO, PFAM and protein count tables...", args.timestamps)

    raw_matrices, cpm_matrices = write_count_matrices(
        samples=samples,
        matrices_dir=matrices_dir,
        normalise_cpm=args.normalise_cpm,
        strict=args.strict,
        time_stamps=args.timestamps,
    )


##### Protein metadata table #####
    # Build and write a non-redundant protein metadata table
    log_msg("Writing protein metadata table...", args.timestamps)

    protein_metadata = build_protein_metadata(
        samples=samples,
        strict=args.strict,
    )

    protein_metadata_path = summaries_dir / "protein_metadata.tsv"

    protein_metadata.to_csv(
        protein_metadata_path,
        sep="\t",
        index=False,
    )

    log_msg(f"    (Written to: {protein_metadata_path})", args.timestamps)


##### Summary of default annotation tables #####
    # Summarise default annotation tables across all samples
    log_msg("Writing summary of FungiFunk default annotation tables...", args.timestamps)

    annotation_summary = build_annotation_summary(
        samples=samples,
        strict=args.strict,
    )

    annotation_summary_path = summaries_dir / "annotation_summary.tsv"

    annotation_summary.to_csv(
        annotation_summary_path,
        sep="\t",
        index=False,
    )

    log_msg(f"    (Written to: {annotation_summary_path})", args.timestamps)


##### Summary of detected features #####
    # Build and write the compact detected-features summary table
    log_msg("Writing summary of detected features...", args.timestamps)

    detected_features_summary = build_detected_features_summary(
        raw_matrices=raw_matrices,
        samples=samples,
    )

    detected_features_summary_path = summaries_dir / "detected_features_summary.tsv"

    detected_features_summary.to_csv(
        detected_features_summary_path,
        sep="\t",
        index=False,
    )

    log_msg(f"    (Written to: {detected_features_summary_path})", args.timestamps)

    ##### Feature metadata #####
    # Build metadata tables to make HTML feature tables easier to interpret
    log_msg("Writing feature metadata tables...", args.timestamps)

    # Define the feature metadata output directory
    feature_metadata_dir = summaries_dir / "feature_metadata"

    # Build and write feature metadata tables
    feature_metadata = build_feature_metadata_tables(
        samples=samples,
        raw_matrices=raw_matrices,
        cpm_matrices=cpm_matrices,
        protein_metadata=protein_metadata,
        feature_metadata_dir=feature_metadata_dir,
        html_top_n=args.html_top_n,
        html_no_top=args.html_no_top,
        calculate_pfam_eggnog_context=args.calculate_pfam_eggnog_context,
        strict=args.strict,
        time_stamps=args.timestamps,
    )


    print(f"########################################## HTML REPORT ##########################################", flush=True)

    # Generate the HTML report unless the user requested to skip it
    if not args.skip_html_report:
        log_msg("Generating HTML report...", args.timestamps)

        generate_html_report(
            outdir=outdir,
            matrices_dir=matrices_dir,
            summaries_dir=summaries_dir,
            sample_summary=sample_summary,
            annotation_summary=annotation_summary,
            detected_features_summary=detected_features_summary,
            raw_matrices=raw_matrices,
            cpm_matrices=cpm_matrices,
            feature_metadata=feature_metadata,
            html_top_n=args.html_top_n,
            html_no_top=args.html_no_top,
            normalise_cpm=args.normalise_cpm,
            calculate_pfam_eggnog_context=args.calculate_pfam_eggnog_context,
            time_stamps=args.timestamps,
            run_info_log=args.run_info_log,
        )

    # Print final success message
    log_msg("Done.", args.timestamps)

    print(f"############################ FungiFunk_summary finished successfully #############################", flush=True)


if __name__ == "__main__":
    main()