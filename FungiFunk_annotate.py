#!/usr/bin/env python3

# Allows modern type hints without problems in some Python versions
from __future__ import annotations

# Lets us define and read command-line arguments
import argparse

# Gives us a cleaner and safer way to handle file paths
from pathlib import Path

# Imports pandas, which we use to read, merge, and save tables
import pandas as pd

# Import datetime to keep track of the time for the timestamps
from datetime import datetime

# Import counter to create the counts tables in a more RAM efficient way
from collections import Counter


# Define a function that reads and returns the command-line arguments
def parse_args() -> argparse.Namespace:

    # Create the argument parser and add a description shown in --help
    parser = argparse.ArgumentParser(
        description=(
            "Merge DIAMOND outputs from R1 and R2, match protein hits against "
            "EggNOG annotations, and generate a full functional annotation table."
        )
    )

    # Add the required argument --r1 for the DIAMOND output of read 1
    parser.add_argument(
        "--r1",
        required=True,
        help="Path to the DIAMOND output file for R1.",
    )

    # Add the required argument --r2 for the DIAMOND output of read 2
    parser.add_argument(
        "--r2",
        required=True,
        help="Path to the DIAMOND output file for R2.",
    )

    # Add the required argument --eggnog_db for the EggNOG annotation table
    parser.add_argument(
        "--eggnog_db",
        required=True,
        help="Path to the EggNOG .emapper.annotations file.",
    )

    # Add the required argument --prefix, which will be used in output filenames
    parser.add_argument(
        "--prefix",
        required=True,
        help="Output prefix.",
    )

    # Add the required argument --outdir, where output files will be written
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )

    # Add the optional argument --write_full_annotation to check if the user used it
    parser.add_argument(
        "--write_full_annotation",
        action="store_true",
        help="Write the full extended annotation table.",
    )

    # Add the optional argument --skip_unannotated_hits to check if the user used it
    parser.add_argument(
        "--skip_unannotated_hits",
        action="store_true",
        help="Skip writing DIAMOND hits without EggNOG annotation.",
    )

    # Add the optional argument --skip_count_tables to check if the user used it
    parser.add_argument(
        "--skip_count_tables",
        action="store_true",
        help="Skip writing all read count tables.",
    )

    # Add the optional argument --skip_protein_read_counts to check if the user used it
    parser.add_argument(
        "--skip_protein_read_counts",
        action="store_true",
        help="Skip writing read counts per reference protein.",
    )

    # Add the optional argument --skip_ko_read_counts to check if the user used it
    parser.add_argument(
        "--skip_ko_read_counts",
        action="store_true",
        help="Skip writing read counts per KEGG ortholog (KO).",
    )

    # Add the optional argument --skip_eggnog_ogs_read_counts to check if the user used it
    parser.add_argument(
        "--skip_eggnog_ogs_read_counts",
        action="store_true",
        help="Skip writing read counts per eggNOG orthologous group.",
    )

    # Add the optional argument --skip_cazy_read_counts to check if the user used it
    parser.add_argument(
        "--skip_cazy_read_counts",
        action="store_true",
        help="Skip writing read counts per CAZy family.",
    )

    # Add the optional argument --skip_pfam_read_counts to check if the user used it
    parser.add_argument(
        "--skip_pfam_read_counts",
        action="store_true",
        help="Skip writing read counts per PFAM domain.",
    )

    # Add the optional argument --compress_ouputs to check if the user used it
    parser.add_argument(
        "--compress_outputs",
        action="store_true",
        help="Compress large annotation outputs as .tsv.gz.",
    )

    # Add the optional arguent --timestamps to check if the user used it
    parser.add_argument(
        "--timestamps",
        action="store_true",
        help="Add timestamps to progress messages"
    )

    # Parse the arguments provided by the user and return them
    return parser.parse_args()


# Define a function to read one DIAMOND tabular output file
def read_diamond_tab(path: str | Path) -> pd.DataFrame:

    # This docstring explains what the function expects and what each DIAMOND column means
    """
    Read a DIAMOND tabular output file with 12 standard BLAST-like columns.

    Column meanings:
    1. query_id   = read ID
    2. protein_id = matched protein ID from the FungiFunk reference database
    3. pident     = percent identity
    4. aln_length = alignment length
    5. mismatch   = number of mismatches
    6. gapopen    = number of gap openings
    7. qstart     = query start
    8. qend       = query end
    9. sstart     = subject start
    10. send      = subject end
    11. evalue    = E-value
    12. bitscore  = bit score
    """

    # Create a list with the column names that we want to assign to the 12 DIAMOND columns
    diamond_cols = [
        "query_id",
        "protein_id",
        "pident",
        "aln_length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
    ]

    # Read the DIAMOND file as a tab-separated table
    # header=None means the file has no header row
    # names=diamond_cols assigns the column names we defined above
    # dtype forces query_id and protein_id to be read as strings, which is important for safe merging
    df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=diamond_cols,
        dtype={
            "query_id": "string",
            "protein_id": "string",
        },
    )

    # Return the DIAMOND dataframe
    return df


# Define a function to read the EggNOG annotation file
def read_eggnog_annotations(path: str | Path) -> pd.DataFrame:

    # This docstring explains the structure of the EggNOG file and why we need special handling
    """
    Read an EggNOG .emapper.annotations file.

    The file starts with commented metadata lines beginning with '##' and one
    header line beginning with '#query'. We strip the leading '#' from the header
    so that the first column becomes 'query'.

    We then rename 'query' to 'protein_id' so that it matches the DIAMOND table
    and can be merged directly.
    """

    # Read the EggNOG file as a tab-separated table
    # comment="#" tells pandas to ignore all lines starting with '#'
    # That means pandas will only read the real data rows, not the metadata/header comments
    # header=None means pandas should not assume that one of those data rows is the header
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
    )

    # Open the EggNOG file manually so we can recover the real header line
    with open(path, "r", encoding="utf-8") as handle:

        # Start by assuming we have not found the header line yet
        header_line = None

        # Read the file line by line.
        for line in handle:

            # Check whether the current line is the real EggNOG header line
            if line.startswith("#query"):

                # Save the header line without the final newline character
                header_line = line.strip()

                # Stop reading once we have found the header
                break

    # If we never found a line starting with #query, the file is not in the expected format
    if header_line is None:
        # Stop the script with a clear error message
        raise ValueError(
            "Could not find the EggNOG header line starting with '#query'."
        )

    # Remove the leading '#' from the header line and split it into column names
    eggnog_cols = header_line.lstrip("#").split("\t")

    # Check that the number of columns in the data matches the number of column names in the header
    if df.shape[1] != len(eggnog_cols):
        # Stop with an informative error if the file structure is inconsistent
        raise ValueError(
            f"EggNOG column mismatch: data has {df.shape[1]} columns but header has {len(eggnog_cols)}."
        )

    # Assign the recovered EggNOG column names to the dataframe
    df.columns = eggnog_cols

    # Rename the EggNOG column 'query' to 'protein_id' so that it matches the DIAMOND dataframe
    # This makes the later merge much cleaner
    df = df.rename(columns={"query": "protein_id"})

    # Force the protein_id column to be string type for safer merging
    df["protein_id"] = df["protein_id"].astype("string")

    # Return the EggNOG dataframe
    return df


# Define a function to count how many rows support each value in a given column
def make_simple_count_table(df: pd.DataFrame, column: str, count_col: str = "read_count") -> pd.DataFrame:
    """
    Count how many rows support each value in a given column.

    Rows with missing values are excluded.
    """
    count_df = (
        df[df[column].notna()]
        .groupby(column)
        .size()
        .reset_index(name=count_col)
        .sort_values(count_col, ascending=False)
    )
    return count_df


# Define a funtion to count how many rows support each individual annotation in a column containing comma-separated values (VERY RAM CONSUMING)
# def make_split_count_table(
#     df: pd.DataFrame,
#     column: str,
#     sep: str = ",",
#     count_col: str = "read_count"
# ) -> pd.DataFrame:
#     """
#     Count how many rows support each individual annotation value in a column
#     containing comma-separated entries.

#     Example:
#     'ko:K00001,ko:K00002' will be split into two separate rows before counting.

#     Rows with missing values or '-' are excluded.
#     """
#     # Keep only rows where the annotation exists and is not just a placeholder
#     sub_df = df[df[column].notna()].copy()
#     sub_df = sub_df[sub_df[column] != "-"]
#     sub_df = sub_df[sub_df[column] != "--"]

#     # Split comma-separated values into Python lists
#     sub_df[column] = sub_df[column].astype(str).str.split(sep)

#     # Expand the lists so each value gets its own row
#     sub_df = sub_df.explode(column)

#     # Remove surrounding spaces
#     sub_df[column] = sub_df[column].str.strip()

#     # Drop empty values that may appear after splitting
#     sub_df = sub_df[sub_df[column] != ""]

#     # Count how many reads support each annotation
#     count_df = (
#         sub_df.groupby(column)
#         .size()
#         .reset_index(name=count_col)
#         .sort_values(count_col, ascending=False)
#     )
#     return count_df

# Define a funtion to count how many rows support each individual annotation in a column containing comma-separated values (RAM friendly)
def make_split_count_table(
    df: pd.DataFrame,
    column: str,
    sep: str = ",",
    count_col: str = "read_count"
) -> pd.DataFrame:
    """
    Count how many rows support each individual annotation value in a column
    containing comma-separated entries.

    Example:
    'ko:K00001,ko:K00002' will be split and counted as two separate values.

    Rows with missing values, empty values, '-' or '--' are excluded.

    This low-memory version avoids copying the full dataframe and avoids
    creating a large exploded intermediate table.
    """
    # Store annotation counts directly, without creating an exploded dataframe
    counter = Counter()

    # Values that should not be counted as real annotations
    invalid_values = {"", "-", "--"}

    # Iterate only over the target annotation column
    for value in df[column]:

        # Skip missing values such as NaN
        if pd.isna(value):
            continue

        # Convert the value to string before splitting
        value = str(value)

        # Skip placeholder values
        if value in invalid_values:
            continue

        # Split comma-separated annotations and count each individual item
        for item in value.split(sep):

            # Remove surrounding whitespace
            item = item.strip()

            # Skip empty or placeholder values after splitting
            if item and item not in invalid_values:
                counter[item] += 1

    # Return an empty table with the expected columns if no counts were found
    if not counter:
        return pd.DataFrame(columns=[column, count_col])

    # Convert the counter into a dataframe
    count_df = pd.DataFrame(
        counter.items(),
        columns=[column, count_col]
    )

    # Sort annotations from most to least supported
    count_df = (
        count_df
        .sort_values(count_col, ascending=False)
        .reset_index(drop=True)
    )

    return count_df


# Extract the last valid annotation value from a separated annotation field
def get_last_valid_annotation_value(value, sep=","):
    """
    Extract the last valid annotation value from a separated annotation field.

    This is used for eggNOG OGs, where the comma-separated values represent
    increasingly specific orthologous groups. Only the most specific valid value
    should be counted.
    """

    # Values that should not be counted as real annotations
    invalid_values = {"", "-", "--", "nan", "NA", "N/A", "None", "none"}

    # Return missing if the value itself is missing
    if pd.isna(value):
        return None

    # Convert the value to text
    value = str(value).strip()

    # Return missing if the complete value is not valid
    if value in invalid_values:
        return None

    # Inspect values from right to left so the most specific valid value is used
    for item in reversed(value.split(sep)):

        # Clean the current item
        item = item.strip()

        # Return the first valid item found from the end
        if item and item not in invalid_values:
            return item

    # Return missing if no valid item was found
    return None


# Count only the most specific annotation value in a separated annotation field
def make_last_split_count_table(
    df: pd.DataFrame,
    column: str,
    sep: str = ",",
    count_col: str = "read_count",
) -> pd.DataFrame:
    """
    Count the last valid separated annotation value in each row.

    This is intended for eggNOG OGs, where the comma-separated values represent
    increasingly specific orthologous groups. Counting only the last valid value
    avoids inflating broad groups such as root or Eukaryota.
    """

    # Store annotation counts directly, without creating an exploded dataframe
    counter = Counter()

    # Iterate only over the target annotation column
    for value in df[column]:

        # Extract the most specific valid annotation value
        item = get_last_valid_annotation_value(
            value,
            sep=sep,
        )

        # Count the value if one was found
        if item is not None:
            counter[item] += 1

    # Return an empty table with the expected columns if no counts were found
    if not counter:
        return pd.DataFrame(columns=[column, count_col])

    # Convert the counter into a dataframe
    count_df = pd.DataFrame(
        counter.items(),
        columns=[column, count_col],
    )

    # Sort annotations from most to least supported
    count_df = (
        count_df
        .sort_values(count_col, ascending=False)
        .reset_index(drop=True)
    )

    return count_df


# Define a function to count how many reads support each protein_id and keep a few useful annotation columns
def make_protein_read_count_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Count reads per protein_id and keep a few useful annotation columns.
    """
    count_df = (
        df.groupby("protein_id")
        .size()
        .reset_index(name="read_count")
    )

    # Select a few annotation columns and keep only one annotation row per protein_id
    annot_cols = ["protein_id", "Description", "Preferred_name", "EC", "KEGG_ko", "KEGG_Pathway"]
    annot_df = df[annot_cols].drop_duplicates(subset=["protein_id"])

    # Merge the read counts per protein with the annotation columns and sort them by descending "read_count"
    out_df = count_df.merge(annot_df, on="protein_id", how="left")
    out_df = out_df.sort_values("read_count", ascending=False)

    # Return the merged and sorted table
    return out_df


# Define a function to write as tsv or tsv.gz depending on the --compress_ouputs flag
def write_output_table(df: pd.DataFrame, path: Path, compress: bool = False) -> Path:
    """
    Write a dataframe as TSV or TSV.GZ depending on the compress flag.

    Returns the final output path.
    """
    if compress:
        compressed_path = Path(f"{path}.gz")
        df.to_csv(compressed_path, sep="\t", index=False, compression="gzip")
        return compressed_path

    df.to_csv(path, sep="\t", index=False)
    return path


# Define a function to add a time stamp to the prints depending ont he --time_stamp flag
def log_msg(message, timestamps=False):
    """
    Print progress messages.

    If timestamps is True, messages are printed with the current time.
    Otherwise, messages are printed normally.
    """
    if timestamps:
        print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)
    else:
        print(message, flush=True)


# Define the main function that runs the whole workflow
def main() -> None:

    # Read all command-line arguments provided by the user
    args = parse_args()

    # Convert the output directory string into a Path object
    outdir = Path(args.outdir)

    # Create the output directory if it does not already exist
    # parents=True allows creation of parent folders if needed
    # exist_ok=True avoids an error if the folder already exists
    outdir.mkdir(parents=True, exist_ok=True)

##### DIAMOND R1 #####
    # Print a progress message so the user knows what the script is doing
    log_msg("Reading DIAMOND output for R1...", args.timestamps)

    # Read the R1 DIAMOND output into a dataframe
    r1_df = read_diamond_tab(args.r1)

    # Add a new column called 'mate' to record that these rows came from R1
    r1_df["mate"] = "R1"


##### DIAMOND R2 #####
    # Print the next progress message
    log_msg("Reading DIAMOND output for R2...", args.timestamps)

    # Read the R2 DIAMOND output into a dataframe
    r2_df = read_diamond_tab(args.r2)

    # Add a new column called 'mate' to record that these rows came from R2
    r2_df["mate"] = "R2"


##### R1 + R2 #####
    # Tell the user that both DIAMOND outputs will now be combined
    log_msg("Concatenating R1 and R2 hits...", args.timestamps)

    # Concatenate the R1 and R2 dataframes vertically into one single dataframe
    # ignore_index=True resets the row numbering in the combined table
    diamond_df = pd.concat([r1_df, r2_df], ignore_index=True)

    # Count how many DIAMOND hits are present before any filtering is applied
    total_hits_before_filtering = len(diamond_df)


##### EGGNOG database annotations #####
    # Print the next progress message
    log_msg("Reading EggNOG database annotations...", args.timestamps)

    # Read the EggNOG annotation file into a dataframe
    eggnog_df = read_eggnog_annotations(args.eggnog_db)


##### DIAMOND + EGGNOG #####
    # Tell the user that the merge step is starting.
    log_msg("Merging DIAMOND hits with EggNOG annotations...", args.timestamps)

    # Merge the DIAMOND table with the EggNOG table using 'protein_id' as the key
    # on="protein_id" means both tables are matched through that column
    # how="left" means we keep all DIAMOND hits even if some do not find an EggNOG annotation
    # validate="many_to_one" checks that many DIAMOND rows may point to one EggNOG protein annotation
    merged_df = diamond_df.merge(
        eggnog_df,
        on="protein_id",
        how="left",
        validate="many_to_one",
    )

    # Rename columns to specify where they come from (diamond or eggnog tsv file)
    merged_df = merged_df.rename(
        columns={
            "evalue_x": "diamond_evalue",
            "evalue_y": "eggnog_evalue",
            "score": "eggnog_score",
        }
    )


##### Output paths #####

    #print(f"", flush=True)
    print(f"##################################### WRITING OUTPUTS #############################################", flush=True)

    # Create the ouput path for the default table.
    default_annotation_path = outdir / f"{args.prefix}_annotation.tsv"

    # Create the output path for the full merged annotation table
    full_annotation_path = outdir / f"{args.prefix}_full_annotation.tsv"

    # Create the output path for the DIAMOND hits that did not match EggNOG
    unannotated_path = outdir / f"{args.prefix}_unannotated_hits.tsv"

    # Create the output path for the protein read count table
    protein_read_counts_path = outdir / f"{args.prefix}_protein_read_counts.tsv"

    # Create the ouput path for the ko read counts table
    ko_read_counts_path = outdir / f"{args.prefix}_ko_read_counts.tsv"

    # Create the output path for the eggNOG_OGs read counts table
    eggnog_ogs_read_counts_path = outdir / f"{args.prefix}_eggnog_ogs_read_counts.tsv"

    # Create the output path for the CAZY read counts table
    cazy_read_counts_path = outdir / f"{args.prefix}_cazy_read_counts.tsv"

    # Create the output path for the PFAM read counts table
    pfam_read_counts_path = outdir / f"{args.prefix}_pfam_read_counts.tsv"


##### Default annotation table #####
    # Create a list of the default columns that constitute the default ouput table
    default_annotation_cols = [
        "query_id",
        "protein_id",
        "diamond_evalue",
        "bitscore",
        "mate",
        "seed_ortholog",
        "eggNOG_OGs",
        "max_annot_lvl",
        "Description",
        "Preferred_name",
        "EC",
        "KEGG_ko",
        "KEGG_Pathway",
        "CAZy",
        "PFAMs",
    ]

    # Create the deafult table that only contains the default coulmns
    default_annotation_df = merged_df[default_annotation_cols].copy()

    # Tell the user that the default annotation table will be written
    log_msg(f"Writing default annotation table...", args.timestamps)

    # Save default columns to a tsv or tsv.gz file
    written_default_path = write_output_table(
        default_annotation_df,
        default_annotation_path,
        compress=args.compress_outputs,
    )

    # Tell the user where the annotation table was written
    log_msg(f"    (Written to: {written_default_path})", args.timestamps)


##### Full annotation table #####
    # If the user has indicated the "--write_full_annotation" then:
    if args.write_full_annotation:
        # Tell the user that the full annotation table will be written
        log_msg(f"Writing full annotation table...", args.timestamps)
        #Save the full annotation table as a tsv or tsv.gz file
        written_full_path = write_output_table(
            merged_df,
            full_annotation_path,
            compress=args.compress_outputs,
        )
        # Tell the user where the full annotation table was wtitten
        log_msg(f"    (Written to: {written_full_path})", args.timestamps)
    else:
        log_msg("Skipping full annotation table.", args.timestamps)


##### Unannotated hits #####
    # Check if the user has requested --skip_unannotated_hits
    if args.skip_unannotated_hits:
        log_msg("Skipping unannotated hits table.", args.timestamps)
    else:
        # Print a progress message for the next output file
        log_msg(f"Writing unannotated hits table...", args.timestamps)
        # Select only the rows where EggNOG annotation is missing
        # We use 'seed_ortholog' as a simple indicator that no annotation was assigned because that column would be NA
        # .copy() makes an explicit copy of the filtered dataframe
        unannotated_df = merged_df[merged_df["seed_ortholog"].isna()].copy()
        # Save the unannotated hits to a tsv or tsv.gz file
        written_unannotated_path = write_output_table(
            unannotated_df,
            unannotated_path,
            compress=args.compress_outputs,
        )
        # Tell the user where the unannotated hist table was written
        log_msg(f"    (Written to: {written_unannotated_path})", args.timestamps)


##### Protein read counts #####
    # Check if the user has requested --skip_count_tables or --skip_protein_read_counts
    if args.skip_count_tables or args.skip_protein_read_counts:
        log_msg("Skipping protein read counts table", args.timestamps)
    else:
        # Print a progress message for the count read count table
        log_msg(f"Writing protein read counts table...", args.timestamps)
        # Use the make_protein_count_table function on the merged dataframe
        protein_read_counts_df = make_protein_read_count_table(merged_df)
        # Save the protein read count table as a tab-separated file
        protein_read_counts_df.to_csv(protein_read_counts_path, sep="\t", index=False)
        # Tell the user where the protein counts table was written
        log_msg(f"    (Written to: {protein_read_counts_path})", args.timestamps)


##### KEGG KO read count #####
    # Check if the user has requested --skip_count_tables or --skip_ko_read_counts
    if args.skip_count_tables or args.skip_ko_read_counts:
        log_msg("Skipping KEGG KO read counts table.", args.timestamps)
    else:
        # Print a progress message for the KEGG KO read count table
        log_msg(f"Writing KEGG KO read count table...", args.timestamps)
        # Use the make_split_count_table function on the merged dataframe for "KEGG_ko" column
        ko_counts_df = make_split_count_table(merged_df, "KEGG_ko")
        # Save the KEGG KO read table as a tab-separated file
        ko_counts_df.to_csv(ko_read_counts_path, sep="\t", index=False)
        # Tell the user where the KEGG KO counts table was written
        log_msg(f"    (Written to: {ko_read_counts_path})", args.timestamps)


##### eggOG_OGs read count #####
    # Check if the user has requested --skip_count_tables or --skip_eggnog_ogs_read_counts
    if args.skip_count_tables or args.skip_eggnog_ogs_read_counts:
        log_msg("Skipping eggNOG OG read counts table.", args.timestamps)
    else:
        # Print a progress message for the eggOG_OGs read count table
        log_msg("Writing eggNOG OG read count table...", args.timestamps)
        # Use the make_last_split_count_table function on the merged dataframe for "eggNOG_OGs" column
        eggnog_ogs_counts_df = make_last_split_count_table(merged_df,"eggNOG_OGs")
        # Save the eggNOG_OGs read table as a tab-separated file
        eggnog_ogs_counts_df.to_csv(eggnog_ogs_read_counts_path, sep="\t", index=False)
        # Tell the user where the eggNOG_OGs counts table was written
        log_msg(f"    (Written to: {eggnog_ogs_read_counts_path})", args.timestamps)


##### CAZy read count #####
    # Check if the user has requested --skip_count_tables or --skip_cazy_read_counts
    if args.skip_count_tables or args.skip_cazy_read_counts:
        log_msg("Skipping CAZy read counts table.", args.timestamps)
    else:
        # Print a progress message for the CAZy read count table
        log_msg("Writing CAZy read count table...", args.timestamps)
        # Use the make_split_count_table function on the merged dataframe for "CAZy" column
        cazy_counts_df = make_split_count_table(merged_df, "CAZy")
        # Save the CAZy read table as a tab-separated file
        cazy_counts_df.to_csv(cazy_read_counts_path, sep="\t", index=False)
        # Tell the user where the CAZy counts table was written
        log_msg(f"    (Written to: {cazy_read_counts_path})", args.timestamps)


##### PFAM read count #####
    # Check if the user has requested --skip_count_tables or --skip_pfam_read_counts
    if args.skip_count_tables or args.skip_pfam_read_counts:
        log_msg("Skipping PFAM read counts table.", args.timestamps)
    else:
        # Print a progress message for the PFAM read count table
        log_msg("Writing PFAM read count table...", args.timestamps)
        # Use the make_split_count_table function on the merged dataframe for "PFAM" column
        pfam_counts_df = make_split_count_table(merged_df, "PFAMs")
        # Save the PFAM read table as a tab-separated file
        pfam_counts_df.to_csv(pfam_read_counts_path, sep="\t", index=False)
        # Tell the user where the PFAM counts table was written
        log_msg(f"    (Written to: {pfam_read_counts_path})", args.timestamps)


##### Final prints #####
    # Print a final message indicating that the script finished successfully
    log_msg("Done.", args.timestamps)

    #print(f"", flush=True)
    print(f"############################################# SUMMARY ############################################", flush=True)

    # Print the total number of combined DIAMOND hits from R1 and R2
    print(f"Total DIAMOND hits: {len(diamond_df)}", flush=True)

    # Count and print how many merged rows received an EggNOG annotation
    # isna() identifies missing values
    # ~ inverts that, so we count non-missing values
    # .sum() counts how many True values there are
    print(f"EggNOG-annotated hits: {(~merged_df['seed_ortholog'].isna()).sum()}", flush=True)

    # Count and print how many merged rows did not receive an EggNOG annotation
    print(f"EggNOG-unannotated hits: {merged_df['seed_ortholog'].isna().sum()}", flush=True)


# Check whether the script is being run directly from the command line
if __name__ == "__main__":

    # If yes, run the main workflow
    main()
