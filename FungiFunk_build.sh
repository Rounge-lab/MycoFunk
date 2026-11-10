#!/usr/bin/env bash

########################################################################################################
# FungiFunk_build.sh
#
# Build paired custom databases for FungiFunk from a protein FASTA file.
#
# This script creates:
#   1. A DIAMOND database for read mapping.
#   2. An EggNOG-mapper annotation table for functional annotation transfer.
#
# Both outputs are generated from the same input protein FASTA file.
# This is important because FungiFunk uses protein IDs from DIAMOND hits
# to retrieve matching annotations from the EggNOG annotation table.
########################################################################################################

# Treat unset variables as errors
# This helps detect typos and missing variables early
set -u

# Define the path to the location of the script
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

# Define the name and version of the tool
TOOL_NAME="FungiFunk_build"
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
    echo ""
    echo "     _____                  _ _____            _        _           _ _     _     "
    echo "    |  ___|   _ _ __   __ _(_)  ___|   _ _ __ | | __   | |__  _   _(_) | __| |    "
    echo "    | |_ | | | |  _ \ / _  | | |_ | | | |  _ \| |/ /   | '_ \| | | | | |/ _  |    "
    echo "    |  _|| |_| | | | | (_| | |  _|| |_| | | | |   <    | |_) | |_| | | | (_| |    "
    echo "    |_|   \__,_|_| |_|\__, |_|_|   \__,_|_| |_|_|\_\___|_.__/ \__,_|_|_|\__,_|    "
    echo "                      |___/                                                       "
    echo ""
    echo "${TOOL_NAME} v${VERSION}: build custom FungiFunk databases from a protein FASTA file"
    echo ""
    echo "FungiFunk_build creates paired DIAMOND and EggNOG annotation databases"
    echo "from the same input protein FASTA file."
    echo ""
    echo "These paired databases can then be passed to FungiFunk.sh using:"
    echo "  -d <custom_diamond_db.dmnd>"
    echo "  -e <custom_eggnog_annotations.emapper.annotations>"
    echo ""
    echo "##########################"
    echo "         USAGE"
    echo "##########################"
    echo "  FungiFunk_build.sh -i <protein_fasta.faa> -o <outdir> [optional arguments]"
    echo ""
    echo "##########################"
    echo "        ARGUMENTS"
    echo "##########################"
    echo "Required:"
    echo "  -i, --input_fasta          Protein FASTA file used to build the paired databases"
    echo "  -o, --outdir               Output directory for the custom database build"
    echo ""
    echo "Optional:"
    echo "  -p, --db_prefix            Prefix/name for the custom database outputs"
    echo "                             [default: fungifunk_custom_db]"
    echo "  -t, --threads <int>        Number of threads [default: 1]"
    echo "  --eggnog_data_dir          Existing EggNOG-mapper data directory"
    echo "                             [default: download to <outdir>/eggnog_data]"
    echo "  --force                    Overwrite existing output files"
    echo "  --timestamps               Add timestamps to FungiFunk progress messages"
    echo "  --version                  Show FungiFunk_build version"
    echo "  -h, --help                 Show this help message"
    echo ""
    echo "##########################"
    echo "         OUTPUTS"
    echo "##########################"
    echo "  <outdir>/diamond/<db_prefix>.dmnd"
    echo "      DIAMOND database built from the input protein FASTA"
    echo ""
    echo "  <outdir>/eggnog/<db_prefix>.emapper.annotations"
    echo "      EggNOG-mapper annotation table generated from the same input protein FASTA"
    echo ""
    echo "  <outdir>/eggnog_data/"
    echo "      EggNOG-mapper data directory, created only if --eggnog_data_dir is not provided"
    echo ""
    echo "  <outdir>/FungiFunk_build_run_info.log"
    echo "      Build metadata log with version, command, inputs, outputs, and parameters"
    echo ""
    echo "  <outdir>/logs/"
    echo "      Logs from DIAMOND, EggNOG data download, and EggNOG-mapper"
    echo ""
    echo "Notes:"
    echo "  The DIAMOND database and EggNOG annotation table are generated as a paired database set."
    echo "  Both are built from the same protein FASTA so that FungiFunk can correctly match DIAMOND hits to EggNOG annotations."
    echo "  If --eggnog_data_dir is provided, EggNOG data download is skipped automatically."
    echo ""
}


########################################################################################################
# Define default values for arguments
########################################################################################################

input_fasta=""
outdir=""
db_prefix="fungifunk_custom_db"
threads="1"
# Store the EggNOG-mapper data directory.
# If empty, the script will use <outdir>/eggnog_data and download the data there.
eggnog_data_dir=""
user_provided_eggnog_data_dir="false"
force="false"
timestamps="false"


########################################################################################################
# Parse command-line arguments
########################################################################################################

# While there are still command-line arguments left, read and classify them one by one
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input_fasta)
            input_fasta="$2"
            shift 2
            ;;
        -o|--outdir)
            outdir="$2"
            shift 2
            ;;
        -p|--db_prefix)
            db_prefix="$2"
            shift 2
            ;;
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        --eggnog_data_dir)
            eggnog_data_dir="$2"
            user_provided_eggnog_data_dir="true"
            shift 2
            ;;
        --force)
            force="true"
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

if [[ -z "$input_fasta" ]]; then
    echo "ERROR: Missing required argument: -i/--input_fasta" >&2
    Help
    exit 1
fi

if [[ -z "$outdir" ]]; then
    echo "ERROR: Missing required argument: -o/--outdir" >&2
    Help
    exit 1
fi

if [[ ! -f "$input_fasta" ]]; then
    echo "ERROR: Input protein FASTA not found: $input_fasta" >&2
    exit 1
fi

if [[ ! -s "$input_fasta" ]]; then
    echo "ERROR: Input protein FASTA is empty: $input_fasta" >&2
    exit 1
fi

if [[ "$db_prefix" == *"/"* ]]; then
    echo "ERROR: --db_prefix must be a simple name, not a path: $db_prefix" >&2
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


########################################################################################################
# Validate required external tools
########################################################################################################

# Check whether DIAMOND is available in the current environment
if ! command -v diamond >/dev/null 2>&1; then
    echo "ERROR: diamond was not found in PATH." >&2
    echo "Please install DIAMOND or activate the correct Conda environment." >&2
    exit 1
fi

# Check whether EggNOG-mapper is available in the current environment
if ! command -v emapper.py >/dev/null 2>&1; then
    echo "ERROR: emapper.py was not found in PATH." >&2
    echo "Please install eggNOG-mapper or activate the correct Conda environment." >&2
    exit 1
fi

# If the user did not provide an existing EggNOG data directory...
if [[ "$user_provided_eggnog_data_dir" == "false" ]]; then
    # Check whether the EggNOG data download script is available.
    if ! command -v download_eggnog_data.py >/dev/null 2>&1; then
        echo "ERROR: download_eggnog_data.py was not found in PATH." >&2
        echo "This command is needed because --eggnog_data_dir was not provided." >&2
        echo "Please install eggNOG-mapper or provide --eggnog_data_dir." >&2
        exit 1
    fi
fi


########################################################################################################
# Define output paths
########################################################################################################

# Define the output directory for the DIAMOND database
diamond_outdir="${outdir}/diamond"

# Define the output directory for EggNOG-mapper results
eggnog_outdir="${outdir}/eggnog"

# Define the output directory for command logs
logs_outdir="${outdir}/logs"

# If the user did not provide an EggNOG data directory...
if [[ "$user_provided_eggnog_data_dir" == "false" ]]; then
    # Use a default EggNOG data directory inside the build output directory
    eggnog_data_dir="${outdir}/eggnog_data"
fi

# Define the DIAMOND database basename
# DIAMOND will add the .dmnd extension automatically
diamond_db_base="${diamond_outdir}/${db_prefix}"

# Define the expected final DIAMOND database file
diamond_db_file="${diamond_db_base}.dmnd"

# Define the EggNOG-mapper output basename
# emapper.py will add .emapper.annotations and other suffixes
eggnog_output_base="${db_prefix}"

# Define the expected final EggNOG annotation table
eggnog_annotation_file="${eggnog_outdir}/${eggnog_output_base}.emapper.annotations"

# Define the build metadata log path
build_info_log="${outdir}/FungiFunk_build_run_info.log"

# Define the DIAMOND makedb log path
diamond_log="${logs_outdir}/${db_prefix}_diamond_makedb.log"

# Define the EggNOG data download log path
eggnog_download_log="${logs_outdir}/${db_prefix}_eggnog_data_download.log"

# Define the EggNOG-mapper run log path
eggnog_mapper_log="${logs_outdir}/${db_prefix}_eggnog_mapper.log"

# Define a path to the EggNOG-mapper download script that will be used during this run
eggnog_download_script=""

# Define a path to a temporary patched copy of download_eggnog_data.py, created only if the installed script contains an outdated EggNOG URL
patched_eggnog_download_script="${logs_outdir}/download_eggnog_data_patched.py"

# Track whether FungiFunk_build had to create a temporary patched copy of download_eggnog_data.py
eggnog_download_script_was_patched="false"


########################################################################################################
# Create output directories
########################################################################################################

mkdir -p "$diamond_outdir"
mkdir -p "$eggnog_outdir"
mkdir -p "$logs_outdir"

# If the user did not provide an existing EggNOG data directory...
if [[ "$user_provided_eggnog_data_dir" == "false" ]]; then
    # Create the EggNOG data directory where EggNOG-mapper data will be downloaded
    mkdir -p "$eggnog_data_dir"
fi


########################################################################################################
# Check for existing outputs
########################################################################################################

# If the user provided an existing EggNOG data directory...
if [[ "$user_provided_eggnog_data_dir" == "true" ]]; then
    # Check whether that directory exists
    if [[ ! -d "$eggnog_data_dir" ]]; then
        echo "ERROR: EggNOG data directory not found: $eggnog_data_dir" >&2
        exit 1
    fi
fi

# Check whether the DIAMOND database already exists and --force was not used
if [[ -f "$diamond_db_file" && "$force" != "true" ]]; then
    echo "ERROR: DIAMOND database already exists: $diamond_db_file" >&2
    echo "Use --force to overwrite existing outputs." >&2
    exit 1
fi

# Check whether the EggNOG annotation file already exists and --force was not used.
if [[ -f "$eggnog_annotation_file" && "$force" != "true" ]]; then
    echo "ERROR: EggNOG annotation file already exists: $eggnog_annotation_file" >&2
    echo "Use --force to overwrite existing outputs." >&2
    exit 1
fi

# If the user requested overwriting existing outputs...
if [[ "$force" == "true" ]]; then
    # Remove the existing DIAMOND database file if present
    rm -f "$diamond_db_file"
    # Remove existing EggNOG-mapper outputs with the same basename if present
    rm -f "${eggnog_outdir}/${eggnog_output_base}".emapper.*
fi

########################################################################################################
# Write build metadata log
########################################################################################################

cat > "$build_info_log" <<EOF
Tool name: ${TOOL_NAME}
Version: ${VERSION}
Run date: $(date)
Working directory: $(pwd)
Command: ${COMMAND_RUN}

Input protein FASTA: ${input_fasta}
Output directory: ${outdir}
Database prefix: ${db_prefix}

Threads: ${threads}
Force overwrite: ${force}
Timestamps: ${timestamps}

DIAMOND output directory: ${diamond_outdir}
DIAMOND database file: ${diamond_db_file}
DIAMOND log: ${diamond_log}

EggNOG output directory: ${eggnog_outdir}
EggNOG annotation file: ${eggnog_annotation_file}
EggNOG data directory: ${eggnog_data_dir}
EggNOG-mapper log: ${eggnog_mapper_log}

User provided EggNOG data directory: ${user_provided_eggnog_data_dir}
EggNOG data download log: ${eggnog_download_log}
EOF


########################################################################################################
# Print start message
########################################################################################################

echo "###################################################################################################"
echo "     _____                  _ _____            _        _           _ _     _     "
echo "    |  ___|   _ _ __   __ _(_)  ___|   _ _ __ | | __   | |__  _   _(_) | __| |    "
echo "    | |_ | | | |  _ \ / _  | | |_ | | | |  _ \| |/ /   | '_ \| | | | | |/ _  |    "
echo "    |  _|| |_| | | | | (_| | |  _|| |_| | | | |   <    | |_) | |_| | | | (_| |    "
echo "    |_|   \__,_|_| |_|\__, |_|_|   \__,_|_| |_|_|\_\___|_.__/ \__,_|_|_|\__,_|    "
echo "                      |___/                                                       "
echo ""
echo "Running ${TOOL_NAME} v${VERSION}"
echo "Input protein FASTA: $input_fasta"
echo "Output directory: $outdir"
echo "Database prefix: $db_prefix"
echo "DIAMOND database: $diamond_db_file"
echo "EggNOG annotation table: $eggnog_annotation_file"
echo "EggNOG data directory: $eggnog_data_dir"
echo "Build metadata log: $build_info_log"
echo ""


########################################################################################################
# Create additional functions
########################################################################################################

log_msg() {
    if [[ "$timestamps" == "true" ]]; then
        echo "[$(date '+%H:%M:%S')] $*"
    else
        echo "$*"
    fi
}


########################################################################################################
# Download EggNOG-mapper data if needed
########################################################################################################

# If the user did not provide an existing EggNOG data directory...
if [[ "$user_provided_eggnog_data_dir" == "false" ]]; then

    echo "---------------------------------------------------------------------------"
    log_msg "[INFO] Preparing EggNOG-mapper data download..."
    log_msg "[INFO] EggNOG data directory: $eggnog_data_dir"
    log_msg "[INFO] Log file: $eggnog_download_log"

    # Locate the download_eggnog_data.py executable available in the active environment
    eggnog_download_script="$(command -v download_eggnog_data.py)"

    log_msg "[INFO] Checking if download_eggnog_data.py uses the current EggNOG download domain..."

    # Check if the installed download_eggnog_data.py still points to the outdated EggNOG download domain
    if grep -q "eggnogdb.embl.de" "$eggnog_download_script"; then
        # Inform the user that a temporary patched copy will be used instead of modifying the original Conda script
        log_msg "[WARNING] download_eggnog_data.py contains the outdated EggNOG download domain."
        log_msg "[WARNING] Creating a patched temporary copy for this FungiFunk_build run."
        log_msg "[WARNING] The original Conda environment will not be modified."

        # Copy the installed download script to the FungiFunk_build logs directory
        cp "$eggnog_download_script" "$patched_eggnog_download_script"

        # Replace the outdated EggNOG domain only in the temporary copy
        sed -i 's|eggnogdb\.embl\.de|eggnog5.embl.de|g' "$patched_eggnog_download_script"

        # Use the patched temporary copy for the EggNOG data download step
        eggnog_download_script="$patched_eggnog_download_script"

        # Record that the download script was patched for this run
        eggnog_download_script_was_patched="true"
    else
        log_msg "[INFO] download_eggnog_data.py appears to use the current EggNOG download domain."
    fi

    log_msg "[INFO] Downloading EggNOG-mapper data..."
    
    # Run the EggNOG data download script
     python "$eggnog_download_script" \
        -y \
        --data_dir "$eggnog_data_dir" \
        > "$eggnog_download_log" 2>&1

    # Update build_info_log with the script used for the download
    echo "EggNOG download script patched: $eggnog_download_script_was_patched" >> "$build_info_log"
    echo "EggNOG download script used: $eggnog_download_script" >> "$build_info_log"

    # Store the exit code from download_eggnog_data.py
    download_exit_code=$?

    # Check whether the download failed.
    if [[ "$download_exit_code" -ne 0 ]]; then
        log_msg "[ERROR] EggNOG data download failed."
        log_msg "        Exit code: $download_exit_code"
        log_msg "        Log file: $eggnog_download_log"
        echo ""
        log_msg "        Last log lines:"
        echo "        ----------------------------------------------------"
        # Print the last 25 lines of the download log with indentation.
        tail -n 25 "$eggnog_download_log" | sed 's/^/        /'
        echo "        ----------------------------------------------------"
        # Exit using the same code returned by the download command.
        exit "$download_exit_code"

    fi

    # Print a success message
    log_msg "[SUCCESS] EggNOG-mapper data download completed."

else
    # If the user provided an EggNOG data directory...
    echo "---------------------------------------------------------------------------"
    log_msg "[INFO] Using existing EggNOG-mapper data directory."
    log_msg "[INFO] EggNOG data directory: $eggnog_data_dir"
    log_msg "[INFO] EggNOG data download skipped."
fi


########################################################################################################
# Build DIAMOND database
########################################################################################################

echo "---------------------------------------------------------------------------"
log_msg "[INFO] Building DIAMOND database..."
log_msg "[INFO] Output: $diamond_db_file"
log_msg "[INFO] Log file: $diamond_log"

# Build the DIAMOND database from the input protein FASTA
diamond makedb \
    --in "$input_fasta" \
    -d "$diamond_db_base" \
    > "$diamond_log" 2>&1

# Store the exit code from DIAMOND
diamond_exit_code=$?

# Check whether DIAMOND failed
if [[ "$diamond_exit_code" -ne 0 ]]; then

    log_msg "[ERROR] DIAMOND database creation failed."
    log_msg "        Exit code: $diamond_exit_code"
    log_msg "        Log file: $diamond_log"
    echo ""
    log_msg "        Last log lines:"
    echo "        ----------------------------------------------------"
    # Print the last 25 lines of the DIAMOND log with indentation
    tail -n 25 "$diamond_log" | sed 's/^/        /'
    echo "        ----------------------------------------------------"
    # Exit using the same code returned by DIAMOND
    exit "$diamond_exit_code"
fi

# Check whether the expected DIAMOND database file was created
if [[ ! -f "$diamond_db_file" ]]; then
    log_msg "[ERROR] DIAMOND finished, but expected output was not found: $diamond_db_file" >&2
    exit 1
fi

# Print a success message
log_msg "[SUCCESS] DIAMOND database created."



########################################################################################################
# Run EggNOG-mapper
########################################################################################################

echo "---------------------------------------------------------------------------"
log_msg "[INFO] Running EggNOG-mapper..."
log_msg "[INFO] Output directory: $eggnog_outdir"
log_msg "[INFO] Expected annotation table: $eggnog_annotation_file"
log_msg "[INFO] Log file: $eggnog_mapper_log"

# Run EggNOG-mapper on the same protein FASTA
emapper.py \
    -i "$input_fasta" \
    --itype proteins \
    --data_dir "$eggnog_data_dir" \
    -m diamond \
    --output "$eggnog_output_base" \
    --output_dir "$eggnog_outdir" \
    --cpu "$threads" \
    > "$eggnog_mapper_log" 2>&1


# Store the exit code from EggNOG-mapper
eggnog_exit_code=$?

# Check whether EggNOG-mapper failed
if [[ "$eggnog_exit_code" -ne 0 ]]; then
    log_msg "[ERROR] EggNOG-mapper failed."
    log_msg "        Exit code: $eggnog_exit_code"
    log_msg "        Log file: $eggnog_mapper_log"
    echo ""
    log_msg "        Last log lines:"
    echo "        ----------------------------------------------------"
    # Print the last 25 lines of the EggNOG-mapper log with indentation.
    tail -n 25 "$eggnog_mapper_log" | sed 's/^/        /'
    echo "        ----------------------------------------------------"
    # Exit using the same code returned by EggNOG-mapper
    exit "$eggnog_exit_code"
fi

# Check whether the expected EggNOG annotation file was created
if [[ ! -f "$eggnog_annotation_file" ]]; then
    log_msg "[ERROR] EggNOG-mapper finished, but expected annotation file was not found: $eggnog_annotation_file" >&2
    exit 1
fi

# Print a success message
log_msg "[SUCCESS] EggNOG annotation table created."


########################################################################################################
# Final message
########################################################################################################

echo "###################################################################################################"
# Print final success message
log_msg "[SUCCESS] FungiFunk custom database build completed."
echo ""

# Introduce output paths
log_msg "Custom database files:"
log_msg "  DIAMOND database:"
log_msg "    $diamond_db_file"
log_msg "  EggNOG annotation table:"
log_msg "    $eggnog_annotation_file"
echo ""

# Print the exact options to use later for FungiFunk.sh
log_msg "Use these options with FungiFunk.sh:"
log_msg "  -d $diamond_db_file"
log_msg "  -e $eggnog_annotation_file"
echo ""

# Print metadata log path
log_msg "Build metadata log:"
log_msg "  $build_info_log"
echo ""

# Print command log paths
log_msg "Command logs:"
log_msg "  DIAMOND:"
log_msg "    $diamond_log"
if [[ "$user_provided_eggnog_data_dir" == "false" ]]; then
    log_msg "  EggNOG data download:"
    log_msg "    $eggnog_download_log"
fi
log_msg "  EggNOG-mapper:"
log_msg "    $eggnog_mapper_log"
echo "###################################################################################################"
echo ""
# Exit with code 0 to indicate complete success
exit 0