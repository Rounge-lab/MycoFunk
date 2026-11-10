#!/usr/bin/env bash

# Treat unset variables as errors
# This helps catch typos in variable names and prevents silent failures
set -u

# Define the path to the location of the script
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

# Define the name and version of the tool
TOOL_NAME="FungiFunk_batch"
VERSION="0.1.9"

# Preserve the original command-line arguments exactly as supplied
ORIGINAL_ARGS=("$@")

# Save a shell-escaped representation of the command for logging
COMMAND_RUN="$(printf '%q ' "$0" "${ORIGINAL_ARGS[@]}")"
COMMAND_RUN="${COMMAND_RUN% }"


########################################################################################################
# Help
########################################################################################################

Help() {
    echo "                                                                                  "
    echo " _____                  _ _____            _        _           _       _         "
    echo "|  ___|   _ _ __   __ _(_)  ___|   _ _ __ | | __   | |__   __ _| |_ ___| |__      "
    echo "| |_ | | | |  _ \ / _  | | |_ | | | |  _ \| |/ /   | '_ \ / _  | __/ __|  _ \     "
    echo "|  _|| |_| | | | | (_| | |  _|| |_| | | | |   <    | |_) | (_| | || (__| | | |    "
    echo "|_|   \__,_|_| |_|\__, |_|_|   \__,_|_| |_|_|\_\___|_.__/ \__,_|\__\___|_| |_|    "
    echo "                  |___/                                                           "
    echo "                                                                                  "
    echo "${TOOL_NAME} v${VERSION}: batch runner for FungiFunk"
    echo ""
    echo "Runs FungiFunk.sh sequentially on multiple paired-end samples listed in a TSV samplesheet."
    echo ""
    echo "##########################"
    echo "         USAGE            "
    echo "##########################"
    echo "  FungiFunk_batch.sh -s <samplesheet.tsv> -o <outdir> [optional arguments]"
    echo ""
    echo "##########################"
    echo "       SAMPLE SHEET       "
    echo "##########################"
    echo "The samplesheet must be tab-separated and contain this header:"
    echo ""
    echo "  prefix<TAB>reads1<TAB>reads2"
    echo ""
    echo "Example:"
    echo ""
    echo "  prefix      reads1                         reads2"
    echo "  sampleA     /path/sampleA_R1.fastq.gz      /path/sampleA_R2.fastq.gz"
    echo "  sampleB     /path/sampleB_R1.fastq.gz      /path/sampleB_R2.fastq.gz"
    echo ""
    echo "##########################"
    echo "        ARGUMENTS         "
    echo "##########################"
    echo "Required:"
    echo "  -s, --samplesheet            Path to samplesheet TSV"
    echo "  -o, --outdir                 Output directory for the full batch"
    echo ""
    echo "Optional:"
    echo "  --fungifunk_script              Path to FungiFunk.sh"
    echo "                                  [default: script_dir/FungiFunk.sh]"
    echo "  -d, --diamond_db                Path to DIAMOND (.dmnd) fungal database file"
    echo "  -e, --eggnog_db                 Path to EggNOG (.emapper.annotations) fungal database file"
    echo "  -t, --threads                   Number of threads per sample [default: 1]"
    echo "  --write_full_annotation         Write the full extended annotation table for each sample"
    echo "  --skip_unannotated_hits         Do not write DIAMOND hits without EggNOG annotation"
    echo "  --skip_count_tables             Do not write any read count tables"
    echo "  --skip_protein_read_counts      Do not write read counts per reference protein"
    echo "  --skip_ko_read_counts           Do not write read counts per KEGG ortholog (KO)"
    echo "  --skip_eggnog_ogs_read_counts   Do not write read counts per eggNOG orthologous group"
    echo "  --skip_cazy_read_counts         Do not write read counts per CAZy family"
    echo "  --skip_pfam_read_counts         Do not write read counts per PFAM domain"
    echo "  --compress_outputs              Compress large annotation outputs as .tsv.gz"
    echo "  --evalue <value>                Maximum E-value to report alignments in DIAMOND"
    echo "  --min-score <value>             Minimum bit score to report alignments in DIAMOND"
    echo "                                  (in DIAMOND, --min-score overrides --evalue)"
    echo "  --id <value>                    Minimum sequence identity percentage to report an alignment in DIAMOND"
    echo "  --query-cover <value>           Minimum query cover percentage to report an alignment in DIAMOND"
    echo "  --top <value>                   Report hits whose scores are within <value>% of the best hit for each query"
    echo "  --stop_on_error                 Stop the batch run if one sample fails"
    echo "  --timestamps                    Add timestamps to FungiFunk progress messages"
    echo "  --version                       Show FungiFunk_batch version"
    echo "  -h, --help                      Show this help message"
    echo ""
    echo "##########################"
    echo "        FILTERING         "
    echo "##########################"
    echo "These options are passed directly to DIAMOND blastx and control which alignments"
    echo "are reported before the annotation step."
    echo ""
    echo "  --evalue <value>        Maximum E-value to report alignments"
    echo "  --min-score <value>     Minimum bit score to report alignments"
    echo "                           Note: in DIAMOND, --min-score overrides --evalue"
    echo "  --id <value>            Minimum sequence identity percentage to report an alignment"
    echo "  --query-cover <value>   Minimum query cover percentage to report an alignment"
    echo "  --top <value>           Report hits whose scores are within <value>% of the best hit for each query"
    echo ""
    echo "If these options are not provided, FungiFunk does not add them to the DIAMOND command."
    echo "This allows DIAMOND to use its own default behaviour for these parameters."
    echo ""
    echo "##########################"
    echo "         OUTPUTS          "
    echo "##########################"
    echo "For each sample:"
    echo "  Default:"
    echo "      <outdir>/<prefix>/                                FungiFunk output directory for that sample"
    echo "      <outdir>/<prefix>/FungiFunk_<prefix>.batch.log    Log capturing stdout and stderr from FungiFunk.sh"
    echo "      FungiFunk_<prefix>_run_info.log                   Metadata log with version, inputs, databases, and parameters"
    echo "      <prefix>_annotation.tsv                           Basic annotation table with core columns:"
    echo "                                                        query_id, protein_id, diamond_evalue, bitscore,"
    echo "                                                        mate, seed_ortholog, eggNOG_OGs, Description,"
    echo "                                                        Preferred_name, EC, KEGG_ko, KEGG_Pathway,"
    echo "                                                        CAZy, PFAMs"
    echo "      <prefix>_unannotated_hits.tsv                     DIAMOND hits with no EggNOG annotation"
    echo "      <prefix>_protein_read_counts.tsv                  Read counts per reference protein"
    echo "      <prefix>_ko_read_counts.tsv                       Read counts per KEGG ortholog (KO)"
    echo "      <prefix>_eggnog_ogs_read_counts.tsv               Read counts per eggNOG orthologous group"
    echo "      <prefix>_cazy_read_counts.tsv                     Read counts per CAZy family"
    echo "      <prefix>_pfam_read_counts.tsv                     Read counts per PFAM domain"
    echo ""
    echo "  Optional:"
    echo "      <prefix>_full_annotation.tsv                      Full extended annotation table"
    echo ""   
    echo ""
    echo "Batch-level files:"
    echo "  <outdir>/FungiFunk_batch_run_info.log                 Metadata log with version, command, samplesheet, and shared parameters"                       
    echo "  <outdir>/FungiFunk_batch_summary.tsv                  Summary table with per-sample status, exit code, and log file paths"
    echo ""
    echo "##########################"
    echo "         NOTES            "
    echo "##########################"
    echo "  By default, the batch continues if a sample fails."
    echo "  Use --stop_on_error to stop at the first failed sample."
    echo ""
}


########################################################################################################
# Define default values for arguments
########################################################################################################

samplesheet=""
outdir=""
threads="1"
# Store the path to the main single-sample FungiFunk.sh script
# By default, assume it is in the same directory as this batch script
fungifunk_script="${script_dir}/FungiFunk.sh"
# Store an optional DIAMOND database path
# Empty means we will let FungiFunk.sh use its own default
diamond_db=""
# Store an optional EggNOG database path
# Empty means we will let FungiFunk.sh use its own default
eggnog_db=""
write_full_annotation="false"
skip_unannotated_hits="false"
skip_count_tables="false"
skip_protein_read_counts="false"
skip_ko_read_counts="false"
skip_eggnog_ogs_read_counts="false"
skip_cazy_read_counts="false"
skip_pfam_read_counts="false"
compress_outputs="false"
# DIAMOND filtering options.
# These options are passed directly to DIAMOND blastx by FungiFunk.sh if provided by the user
evalue=""
min_score=""
id=""
query_cover=""
top=""
# Store whether the batch should stop after the first failed sample
# Default is false, meaning the batch continues after failures
stop_on_error="false"
timestamps="false"


########################################################################################################
# Parse command-line arguments
########################################################################################################

# While there are still command-line arguments left, read and classify them one by one
while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--samplesheet)
            samplesheet="$2"
            shift 2
            ;;
        -o|--outdir)
            outdir="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        --fungifunk_script)
            fungifunk_script="$2"
            shift 2
            ;;
        -d|--diamond_db)
            diamond_db="$2"
            shift 2
            ;;
        -e|--eggnog_db)
            eggnog_db="$2"
            shift 2
            ;;
        --write_full_annotation)
            write_full_annotation="true"
            shift 1
            ;;
        --skip_unannotated_hits)
            skip_unannotated_hits="true"
            shift
            ;;
        --skip_count_tables)
            skip_count_tables="true"
            shift
            ;;
        --skip_protein_read_counts)
            skip_protein_read_counts="true"
            shift
            ;;
        --skip_ko_read_counts)
            skip_ko_read_counts="true"
            shift
            ;;
        --skip_eggnog_ogs_read_counts)
            skip_eggnog_ogs_read_counts="true"
            shift
            ;;
        --skip_cazy_read_counts)
            skip_cazy_read_counts="true"
            shift
            ;;
        --skip_pfam_read_counts)
            skip_pfam_read_counts="true"
            shift
            ;;
        -c|--compress_outputs)
            compress_outputs="true"
            shift 1
            ;;
        --evalue)
            evalue="$2"
            shift 2
            ;;
        --min-score)
            min_score="$2"
            shift 2
            ;;
        --id)
            id="$2"
            shift 2
            ;;
        --query-cover)
            query_cover="$2"
            shift 2
            ;;
        --top)
            top="$2"
            shift 2
            ;;
        --stop_on_error)
            stop_on_error="true"
            shift 1
            ;;
        --timestamps)
            timestamps="true"
            shift 1
            ;;
        --version)
            echo "${TOOL_NAME} v${VERSION}"
            exit 0
            ;;
        -h|--help)
            Help
            exit 0
            ;;
        *)
        # If the current argument does not match any option name, stop and return an error
            echo "ERROR: Unknown argument: $1" >&2
            Help
            exit 1
            ;;
    esac
done


########################################################################################################
# Validate arguments
########################################################################################################

if [[ -z "$samplesheet" ]]; then
    echo "ERROR: Missing required argument: -s/--samplesheet" >&2
    Help
    exit 1
fi

if [[ -z "$outdir" ]]; then
    echo "ERROR: Missing required argument: -o/--outdir" >&2
    Help
    exit 1
fi

if [[ ! -f "$samplesheet" ]]; then
    echo "ERROR: Samplesheet not found: $samplesheet" >&2
    exit 1
fi

if [[ ! -f "$fungifunk_script" ]]; then
    echo "ERROR: FungiFunk.sh not found: $fungifunk_script" >&2r.
    exit 1
fi

if [[ ! -r "$fungifunk_script" ]]; then
    echo "ERROR: FungiFunk.sh is not readable: $fungifunk_script" >&2
    exit 1
fi

if [[ -n "$diamond_db" && ! -f "$diamond_db" ]]; then
    echo "ERROR: DIAMOND database not found: $diamond_db" >&2
    exit 1
fi

if [[ -n "$eggnog_db" && ! -f "$eggnog_db" ]]; then
    echo "ERROR: EggNOG database not found: $eggnog_db" >&2
    exit 1
fi

if ! [[ "$threads" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --threads must be a positive integer: $threads" >&2
    exit 1
fi

if [[ "$threads" -lt 1 ]]; then
    echo "ERROR: --threads must be at least 1" >&2
    exit 1
fi

# Check that --evalue is a valid numeric value if it was provided
if [[ -n "$evalue" ]] && ! [[ "$evalue" =~ ^[0-9]+([.][0-9]+)?([eE]-?[0-9]+)?$ ]]; then
    echo "ERROR: --evalue must be a valid numeric value: $evalue" >&2
    exit 1
fi

# Check that --min-score is a valid numeric value if it was provided
if [[ -n "$min_score" ]] && ! [[ "$min_score" =~ ^[0-9]+([.][0-9]+)?([eE]-?[0-9]+)?$ ]]; then
    echo "ERROR: --min-score must be a valid numeric value: $min_score" >&2
    exit 1
fi

if [[ -n "$min_score" && -n "$evalue" ]]; then
    echo "[WARNING] Both --min-score and --evalue were provided."
    echo "[WARNING] In DIAMOND, --min-score overrides the E-value setting."
fi

# Check that --id is a valid numeric value if it was provided
if [[ -n "$id" ]] && ! [[ "$id" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: --id must be a valid numeric value: $id" >&2
    exit 1
fi

# Check that --id is between 0 and 100 if it was provided
if [[ -n "$id" ]] && awk "BEGIN {exit !($id < 0 || $id > 100)}"; then
    echo "ERROR: --id must be between 0 and 100: $id" >&2
    exit 1
fi

# Check that --query-cover is a valid numeric value if it was provided
if [[ -n "$query_cover" ]] && ! [[ "$query_cover" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
    echo "ERROR: --query-cover must be a valid numeric value: $query_cover" >&2
    exit 1
fi

# Check that --query-cover is between 0 and 100 if it was provided
if [[ -n "$query_cover" ]] && awk "BEGIN {exit !($query_cover < 0 || $query_cover > 100)}"; then
    echo "ERROR: --query-cover must be between 0 and 100: $query_cover" >&2
    exit 1
fi

# Check that --top is a valid numeric value if it was provided
if [[ -n "$top" ]] && ! [[ "$top" =~ ^[0-9]+([.][0-9]+)?([eE]-?[0-9]+)?$ ]]; then
    echo "ERROR: --top must be a valid numeric value: $top" >&2
    exit 1
fi

# Check that --top is between 0 and 100 if it was provided
if [[ -n "$top" ]] && awk "BEGIN {exit !($top < 0 || $top > 100)}"; then
    echo "ERROR: --top must be between 0 and 100: $top" >&2
    exit 1
fi



########################################################################################################
# Check samplesheet header
########################################################################################################

# Read the first line of the samplesheet
header="$(head -n 1 "$samplesheet")"

# Define the expected header
expected_header=$'prefix\treads1\treads2'

# Check whether the samplesheet header is exactly correct, if not print an error
if [[ "$header" != "$expected_header" ]]; then
    echo "ERROR: Invalid samplesheet header." >&2
    echo "Expected tab-separated header:" >&2
    echo -e "prefix\treads1\treads2" >&2
    echo "" >&2
    echo "Found:" >&2
    echo "$header" >&2
    exit 1
fi


########################################################################################################
# Prepare batch output directories, files and log
########################################################################################################

# Create the batch output directory if it does not already exist
mkdir -p "$outdir"

# Define the path to the batch run log
batch_log="${outdir}/FungiFunk_batch_run_info.log"

# Define the path to the batch summary table
batch_summary="${outdir}/FungiFunk_batch_summary.tsv"

# Write batch run metadata to log file
cat > "$batch_log" <<EOF
Tool name: ${TOOL_NAME}
Version: ${VERSION}
Run date: $(date '+%Y-%m-%d %H:%M:%S')
Working directory: $(pwd)
Command: ${COMMAND_RUN}

Samples information:
    Samplesheet: ${samplesheet}
    Output directory: ${outdir}
    FungiFunk script: ${fungifunk_script}

Databases:
    DIAMOND database: ${diamond_db:-FungiFunk.sh default}
    EggNOG database: ${eggnog_db:-FungiFunk.sh default}

Run options:
    Threads per sample: ${threads}
    Write full annotation table: ${write_full_annotation}
    Skip unannotated hits table: ${skip_unannotated_hits}
    Skip all count tables: ${skip_count_tables}
    Skip protein read counts table: ${skip_protein_read_counts}
    Skip KO read counts table: ${skip_ko_read_counts}
    Skip eggNOG OG read counts table: ${skip_eggnog_ogs_read_counts}
    Skip CAZy read counts table: ${skip_cazy_read_counts}
    Skip PFAM read counts table: ${skip_pfam_read_counts}
    Compress outputs: ${compress_outputs}
    Stop on error: ${stop_on_error}
    Timestamps: ${timestamps}

Filtering:
    E-value filter: ${evalue:-DIAMOND default = 0.001}
    Min-score filter: ${min_score:-not used}
    Identity filter: ${id:-not used}
    Query-cover filter: ${query_cover:-not used}
    Top filter: ${top:-not used}

EOF

# Create the batch summary table and write its header
echo -e "prefix\treads1\treads2\toutdir\tstatus\texit_code\treason\tlog_file" > "$batch_summary"



########################################################################################################
# Build common FungiFunk.sh arguments
########################################################################################################

# Create an empty bash array to store optional arguments shared by all samples
common_args=()

# If a custom DIAMOND database was provided...
if [[ -n "$diamond_db" ]]; then
    # Add the DIAMOND database argument to the common arguments array
    common_args+=( -d "$diamond_db" )
fi

if [[ -n "$eggnog_db" ]]; then
    common_args+=( -e "$eggnog_db" )
fi

 # If the user requested full annotation tables...
if [[ "$write_full_annotation" == "true" ]]; then
    # Add the flag to the common arguments array
    common_args+=( --write_full_annotation )
fi

if [[ "$skip_unannotated_hits" == "true" ]]; then
    common_args+=( --skip_unannotated_hits )
fi

if [[ "$skip_count_tables" == "true" ]]; then
    common_args+=( --skip_count_tables )
fi

if [[ "$skip_protein_read_counts" == "true" ]]; then
    common_args+=( --skip_protein_read_counts )
fi

if [[ "$skip_ko_read_counts" == "true" ]]; then
    common_args+=( --skip_ko_read_counts )
fi

if [[ "$skip_eggnog_ogs_read_counts" == "true" ]]; then
    common_args+=( --skip_eggnog_ogs_read_counts )
fi

if [[ "$skip_cazy_read_counts" == "true" ]]; then
    common_args+=( --skip_cazy_read_counts )
fi

if [[ "$skip_pfam_read_counts" == "true" ]]; then
    common_args+=( --skip_pfam_read_counts )
fi

if [[ "$compress_outputs" == "true" ]]; then
    common_args+=( --compress_outputs )
fi

if [[ -n "$evalue" ]]; then
    common_args+=( --evalue "$evalue" )
fi

if [[ -n "$min_score" ]]; then
    common_args+=( --min_score "$min_score" )
fi

if [[ -n "$id" ]]; then
    common_args+=( --id "$id" )
fi

if [[ -n "$query_cover" ]]; then
    common_args+=( --query-cover "$query_cover" )
fi

if [[ -n "$top" ]]; then
    common_args+=( --top "$top" )
fi

if [[ "$timestamps" == "true" ]]; then
    common_args+=( --timestamps )
fi

########################################################################################################
# Print batch start messages
########################################################################################################

echo "###################################################################################################"
echo " _____                  _ _____            _        _           _       _         "
echo "|  ___|   _ _ __   __ _(_)  ___|   _ _ __ | | __   | |__   __ _| |_ ___| |__      "
echo "| |_ | | | |  _ \ / _  | | |_ | | | |  _ \| |/ /   | '_ \ / _  | __/ __|  _ \     "
echo "|  _|| |_| | | | | (_| | |  _|| |_| | | | |   <    | |_) | (_| | || (__| | | |    "
echo "|_|   \__,_|_| |_|\__, |_|_|   \__,_|_| |_|_|\_\___|_.__/ \__,_|\__\___|_| |_|    "
echo "                  |___/                                                           "
echo "                                                                                  "
echo ""
echo "Running ${TOOL_NAME} v${VERSION}"
echo "Samplesheet: $samplesheet"
echo "Output directory: $outdir"
echo "Batch summary: $batch_summary"
echo "Batch run info: $batch_log"
echo ""


########################################################################################################
# Initialize counters
########################################################################################################

# Count how many sample rows were processed
total_samples=0

# Count how many samples completed successfully
successful_samples=0

# Count how many samples failed
failed_samples=0


########################################################################################################
# Run FungiFunk.sh for each sample in the samplesheet
########################################################################################################

# Read the samplesheet line by line
# IFS=$'\t' tells bash to split columns using tabs
# read -r prevents backslash interpretation
# prefix, reads1, and reads2 receive the first three columns
# extra receives any unexpected fourth column or additional content
while IFS=$'\t' read -r prefix reads1 reads2 extra; do

    # If prefix, reads1, and reads2 are all empty, treat this as an empty line
    if [[ -z "${prefix:-}" && -z "${reads1:-}" && -z "${reads2:-}" ]]; then
        # Skip empty lines
        continue
    fi

    # Increase the total sample counter by one
    total_samples=$((total_samples + 1))

    # Define the output directory for this sample
    sample_outdir="${outdir}/${prefix}"

    # Define the log file for this sample
    # This log will capture stdout and stderr from FungiFunk.sh
    sample_log="${sample_outdir}/FungiFunk_${prefix}.batch.log"

    # Create the sample output directory if it does not already exist
    mkdir -p "$sample_outdir"

    # Print sample data
    echo "---------------------------------------------------------------------------------------------------"
    echo "[INFO] Starting sample: $prefix"
    echo "[INFO] reads1: $reads1"
    echo "[INFO] reads2: $reads2"
    echo "[INFO] sample output: $sample_outdir"
    echo "[INFO] sample log: $sample_log"


    ####################################################################################################
    # Per-sample validation
    ####################################################################################################

    # Check whether the row has more than three columns
    # This usually means the samplesheet is malformed
    if [[ -n "${extra:-}" ]]; then
        # Store a simple reason for the failure
        reason="Samplesheet row has more than 3 columns"
        # Print a clear error message
        echo "[ERROR] Sample $prefix failed before running FungiFunk."
        # Print the reason.
        echo "        Reason: $reason"
        # Add a FAILED entry to the batch summary
        # exit_code is NA because FungiFunk.sh was not run
        echo -e "${prefix}\t${reads1}\t${reads2}\t${sample_outdir}\tFAILED\tNA\t${reason}\t${sample_log}" >> "$batch_summary"
        # Increase the failed sample counter by one
        failed_samples=$((failed_samples + 1))

        # If the user requested stopping on the first error...
        if [[ "$stop_on_error" == "true" ]]; then
            # Print a stopping message.
            echo "[ERROR] --stop_on_error was used. Stopping batch."
            # Stop the batch with an error code
            exit 1
        fi
        # Continue to the next sample if stop_on_error is not enabled
        continue
    fi

    # Check whether any required samplesheet field is missing
    if [[ -z "$prefix" || -z "$reads1" || -z "$reads2" ]]; then
        # Store a simple reason for the failure
        reason="Missing prefix, reads1, or reads2 in samplesheet row"
        # Print a clear error message.
        echo "[ERROR] Sample row failed before running FungiFunk."
        # Print the reason.
        echo "        Reason: $reason"
        # Add a FAILED entry to the batch summary
        # Use NA for fields that may be missing
        echo -e "${prefix:-NA}\t${reads1:-NA}\t${reads2:-NA}\t${sample_outdir:-NA}\tFAILED\tNA\t${reason}\t${sample_log:-NA}" >> "$batch_summary"
        # Increase the failed sample counter by one
        failed_samples=$((failed_samples + 1))

        # If the user requested stopping on the first error...
        if [[ "$stop_on_error" == "true" ]]; then
            # Print a stopping message.
            echo "[ERROR] --stop_on_error was used. Stopping batch."
        # Stop the batch with an error code
            exit 1
        fi
        # Continue to the next sample if stop_on_error is not enabled
        continue
    fi

    # Check whether the reads1 file exists
    if [[ ! -f "$reads1" ]]; then
        # Store a simple reason for the failure
        reason="reads1 file not found"
        # Print a clear error message
        echo "[ERROR] Sample $prefix failed before running FungiFunk."
        # Print the reason and the missing file path
        echo "        Reason: $reason: $reads1"
        # Add a FAILED entry to the batch summary
        # exit_code is NA because FungiFunk.sh was not run
        echo -e "${prefix}\t${reads1}\t${reads2}\t${sample_outdir}\tFAILED\tNA\t${reason}\t${sample_log}" >> "$batch_summary"
        # Increase the failed sample counter by one
        failed_samples=$((failed_samples + 1))

        # If the user requested stopping on the first error...
        if [[ "$stop_on_error" == "true" ]]; then
            # Print a stopping message
            echo "[ERROR] --stop_on_error was used. Stopping batch."
            # Stop the batch with an error code
            exit 1

        fi
        # Continue to the next sample if stop_on_error is not enabled
        continue
    fi

    # Check whether the reads2 file exists
    if [[ ! -f "$reads2" ]]; then
        # Store a simple reason for the failure
        reason="reads2 file not found"
        # Print a clear error message
        echo "[ERROR] Sample $prefix failed before running FungiFunk."
        # Print the reason and the missing file path
        echo "        Reason: $reason: $reads2"
        # Add a FAILED entry to the batch summary
        # exit_code is NA because FungiFunk.sh was not run
        echo -e "${prefix}\t${reads1}\t${reads2}\t${sample_outdir}\tFAILED\tNA\t${reason}\t${sample_log}" >> "$batch_summary"
        # Increase the failed sample counter by one
        failed_samples=$((failed_samples + 1))

        # If the user requested stopping on the first error...
        if [[ "$stop_on_error" == "true" ]]; then
            # Print a stopping message.
            echo "[ERROR] --stop_on_error was used. Stopping batch."
            # Stop the batch with an error code
            exit 1
        fi
        # Continue to the next sample if stop_on_error is not enabled
        continue
    fi


    ####################################################################################################
    # Run FungiFunk.sh for this sample
    ####################################################################################################

    # Run the single-sample FungiFunk.sh script
    # -1 passes reads1
    # -2 passes reads2
    # -p passes the sample prefix
    # -o passes the sample-specific output directory
    # -t passes the number of threads
    # "${common_args[@]}" passes all optional common batch settings safely
    # > "$sample_log" redirects stdout to FungiFunk_${prefix}.batch.log
    # 2>&1 redirects stderr to the same log file
    bash "$fungifunk_script" \
        -1 "$reads1" \
        -2 "$reads2" \
        -p "$prefix" \
        -o "$sample_outdir" \
        -t "$threads" \
        "${common_args[@]}" \
        > "$sample_log" 2>&1

    # Store the exit code from FungiFunk.sh
    exit_code=$?

    # Check whether FungiFunk.sh completed successfully
    if [[ "$exit_code" -eq 0 ]]; then
        # Store a success reason
        reason="Completed successfully"
        # Print a success message
        echo "[SUCCESS] Sample $prefix completed successfully."
        # Add a SUCCESS entry to the batch summary
        echo -e "${prefix}\t${reads1}\t${reads2}\t${sample_outdir}\tSUCCESS\t${exit_code}\t${reason}\t${sample_log}" >> "$batch_summary"
        # Increase the successful sample counter by one
        successful_samples=$((successful_samples + 1))
    # If FungiFunk.sh returned a non-zero exit code...
    else
        # Store a general failure reason
        reason="FungiFunk.sh exited with non-zero status"
        # Print a clear error message
        echo "[ERROR] Sample $prefix failed."
        echo "        Exit code: $exit_code"
        echo "        Reason: $reason"
        echo "        Log file: $sample_log"
        echo ""
        echo "        Last log lines:"
        echo "        ----------------------------------------------------"
        # Print the last 25 lines of the sample log
        # sed adds indentation to each line for readability
        tail -n 25 "$sample_log" | sed 's/^/        /'
        echo "        ----------------------------------------------------"
        echo ""
        # Add a FAILED entry to the batch summary
        echo -e "${prefix}\t${reads1}\t${reads2}\t${sample_outdir}\tFAILED\t${exit_code}\t${reason}\t${sample_log}" >> "$batch_summary"
        # Increase the failed sample counter by one
        failed_samples=$((failed_samples + 1))

        # If the user requested stopping on the first error...
        if [[ "$stop_on_error" == "true" ]]; then
            # Print a stopping message
            echo "[ERROR] --stop_on_error was used. Stopping batch."
            # Stop the batch using the same exit code returned by FungiFunk.sh
            exit "$exit_code"
        fi
    fi
# End of the while loop
# tail -n +2 skips the header line of the samplesheet
# Feed all samplesheet rows except the header into the while loop
# This avoids using a pipe, which can make counter variables unavailable after the loop in bash
done < <(tail -n +2 "$samplesheet")


########################################################################################################
# Final message
########################################################################################################
echo ""
echo "###################################################################################################"
echo "[INFO] Batch finished."
echo "[INFO] Total samples processed: $total_samples"
echo "[INFO] Successful samples: $successful_samples"
echo "[INFO] Failed samples: $failed_samples"
echo "[INFO] Summary file: $batch_summary"
echo "[INFO] Batch run info: $batch_log"
echo "###################################################################################################"
echo ""

# Check whether at least one sample failed
if [[ "$failed_samples" -gt 0 ]]; then
    # Exit with code 1 if any sample failed
    # This is useful for automated workflows because the batch is marked as not fully successful
    exit 1
# If no samples failed...
else
    # Exit with code 0 to indicate complete success
    exit 0
fi