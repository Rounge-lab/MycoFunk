#!/usr/bin/env bash

# Stop the script on command failures, undefined variables or hidden errors in pipes
set -euo pipefail

# Define the path to the location of the script
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"

# Define the name and version of the tool
TOOL_NAME="MycoFunk"
VERSION="0.2.0"

# Preserve the original command-line arguments exactly as supplied
ORIGINAL_ARGS=("$@")

# Save a shell-escaped representation of the command for logging
COMMAND_RUN="$(printf '%q ' "$0" "${ORIGINAL_ARGS[@]}")"
COMMAND_RUN="${COMMAND_RUN% }"


########################################################################################################
# Help
########################################################################################################

Help() {
    echo "                                                      "
    echo " __  __                 _____            _            "
    echo "|  \/  |_   _  ___ ___ |  ___|   _ _ __ | | __        "
    echo "| |\/| | | | |/ __/ _ \| |_ | | | | '_ \| |/ /        "
    echo "| |  | | |_| | (_| (_) |  _|| |_| | | | |   <         "
    echo "|_|  |_|\__, |\___\___/|_|   \__,_|_| |_|_|\_\        "
    echo "        |___/                                         "
    echo "                                                      "
    echo "${TOOL_NAME} v${VERSION}: fungal functional profiling from paired-end reads"
    echo ""
    echo "MycoFunk maps paired-end reads against a fungal protein reference database with DIAMOND,"
    echo "then transfers functional annotation from EggNOG to generate annotation and"
    echo "read-count summary tables."
    echo ""
    echo ""
    echo "##########################"
    echo "         USAGE            "
    echo "##########################"
    echo "  MycoFunk.sh -1 <reads1> -2 <reads2> -p <prefix> -o <outdir> [optional arguments]"
    echo ""
    echo ""
    echo "##########################"
    echo "         ARGUMENTS        "
    echo "##########################"
    echo "Required:"
    echo "  -1, --reads1 <file>             Path to paired-end reads 1 (FASTQ/FASTQ.GZ)"
    echo "  -2, --reads2 <file>             Path to paired-end reads 2 (FASTQ/FASTQ.GZ)"
    echo "  -p, --prefix                    Sample prefix"
    echo "  -o, --outdir <dir>              Output directory"
    echo ""
    echo "Optional:"
    echo "  -d, --diamond_db <dir>          Path to DIAMOND (.dmnd) fungal database file"
    echo "                                  [default = script_dir/MycoFunk_diamond_db/mycofunk_diamond_db.dmnd]"
    echo "  -e, --eggnog_db <dir>           Path to EggNOG (.emapper.annotations) fungal database file"
    echo "                                  [default = script_dir/MycoFunk_eggnog_db/mycofunk_eggnog_db.emapper.annotations]"
    echo "  -t, --threads <int>             Number of threads [default = 1]"
    echo "  --write_full_annotation         Write the full extended annotation table [default = disabled]"
    echo "  --skip_unannotated_hits         Do not write DIAMOND hits without EggNOG annotation"
    echo "  --skip_count_tables             Do not write any read count tables"
    echo "  --skip_protein_read_counts      Do not write read counts per reference protein"
    echo "  --skip_ko_read_counts           Do not write read counts per KEGG ortholog (KO)"
    echo "  --skip_eggnog_ogs_read_counts   Do not write read counts per eggNOG orthologous group"
    echo "  --skip_cazy_read_counts         Do not write read counts per CAZy family"
    echo "  --skip_pfam_read_counts         Do not write read counts per PFAM domain"
    echo "  --compress_outputs              Compress large annotation outputs as .tsv.gz [default = disabled]"
    echo "  --evalue <value>                Maximum E-value to report alignments in DIAMOND [default = 0.001]"
    echo "  --min-score <value>             Minimum bit score to report alignments in DIAMOND"
    echo "                                  (in DIAMOND, --min-score overrides --evalue)"
    echo "  --id <value>                    Minimum sequence identity percentage to report an alignment in DIAMOND [deafult = 88]"
    echo "  --query-cover <value>           Minimum query cover percentage to report an alignment in DIAMOND"
    echo "  --top <value>                   Report hits whose scores are within <value>% of the best hit for each query"
    echo "  --timestamps                    Add timestamps to MycoFunk progress messages"
    echo "  --version                       Show MycoFunk version and exit"
    echo "  -h, --help                      Show this help message and exit"
    echo ""
    echo ""
    echo "##########################"
    echo "        FILTERING         "
    echo "##########################"
    echo "These options are passed directly to DIAMOND blastx and control which alignments"
    echo "are reported before the annotation step."
    echo ""
    echo "  --evalue <value>                Maximum E-value to report alignments"
    echo "  --min-score <value>             Minimum bit score to report alignments"
    echo "                                   Note: in DIAMOND, --min-score overrides --evalue"
    echo "  --id <value>                    Minimum sequence identity percentage to report an alignment"
    echo "  --query-cover <value>           Minimum query cover percentage to report an alignment"
    echo "  --top <value>                   Report hits whose scores are within <value>% of the best hit for each query"
    echo ""
    echo "If these options are not provided, MycoFunk does not add them to the DIAMOND command."
    echo "This allows DIAMOND to use its own default behaviour for these parameters."
    echo ""
    echo ""
    echo "##########################"
    echo "         OUTPUTS          "
    echo "##########################"
    echo "Default:"
    echo "  MycoFunk_<prefix>_run_info.log             Run metadata log with version, inputs, databases, and parameters"
    echo "  <prefix>_input_read_counts.tsv              Input read counts for R1 and R2"
    echo ""
    echo "  diamond/"
    echo "      <prefix>_R1.diamond.out                 DIAMOND blastx output for R1 reads"
    echo "      <prefix>_R1.diamond.log                 DIAMOND log file for R1 reads"
    echo "      <prefix>_R2.diamond.out                 DIAMOND blastx output for R2 reads"
    echo "      <prefix>_R2.diamond.log                 DIAMOND log file for R2 reads"
    echo ""
    echo "  annotation/"
    echo "      <prefix>_annotation.tsv                 Basic annotation table with core columns:"
    echo "                                              query_id, protein_id, diamond_evalue, bitscore,"
    echo "                                              mate, seed_ortholog, eggNOG_OGs, Description,"
    echo "                                              Preferred_name, EC, KEGG_ko, KEGG_Pathway,"
    echo "                                              CAZy, PFAMs"
    echo "      <prefix>_unannotated_hits.tsv           DIAMOND hits with no EggNOG annotation"
    echo "      <prefix>_protein_read_counts.tsv        Read counts per reference protein"
    echo "      <prefix>_ko_read_counts.tsv             Read counts per KEGG ortholog (KO)"
    echo "      <prefix>_eggnog_ogs_read_counts.tsv     Read counts per eggNOG orthologous group"
    echo "      <prefix>_cazy_read_counts.tsv           Read counts per CAZy family"
    echo "      <prefix>_pfam_read_counts.tsv           Read counts per PFAM domain"
    echo ""
    echo "Optional:"
    echo "  annotation/"
    echo "      <prefix>_full_annotation.tsv            Full extended annotation table"
    echo "                                              (only if --write_full_annotation is used)"
    echo ""
    echo ""
    echo "##########################"
    echo "         NOTES            "
    echo "##########################"
    echo "  Count tables are based on read-level evidence (read_count), not fragment-level counts."
    echo "  If --compress_outputs is used, annotation tables are written as .tsv.gz."
    echo "  Read-count summary tables remain uncompressed."
    echo ""
}


########################################################################################################
# Define default values for arguments
########################################################################################################

reads1=""
reads2=""
prefix=""
outdir=""
diamond_db="${script_dir}/MycoFunk_diamond_db/mycofunk_diamond_db.dmnd"
eggnog_db="${script_dir}/MycoFunk_eggnog_db/mycofunk_eggnog_db.emapper.annotations"
threads="1"
write_full_annotation="false"
skip_unannotated_hits="false"
skip_count_tables="false"
skip_protein_read_counts="false"
skip_ko_read_counts="false"
skip_eggnog_ogs_read_counts="false"
skip_cazy_read_counts="false"
skip_pfam_read_counts="false"
compress_outputs="false"
evalue=""
min_score=""
id="88"
query_cover=""
top=""
timestamps="false"


########################################################################################################
# Parse arguments
########################################################################################################

# While there are still command-line arguments left, read and classify them one by one
while [[ $# -gt 0 ]]; do
    case "$1" in
        # If the current argument matches reads1, store its value and move to the next option
        -1|--reads1)
            reads1="$2"
            shift 2
            ;;
        -2|--reads2)
            reads2="$2"
            shift 2
            ;;
        -p|--prefix)
            prefix="$2"
            shift 2
            ;;
        -o|--outdir)
            outdir="$2"
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
        -t|--threads)
            threads="$2"
            shift 2
            ;;
        --write_full_annotation)
            write_full_annotation="true"
            shift 1
            ;;
        --skip_unannotated_hits)
            skip_unannotated_hits="true"
            shift 1
            ;;
        --skip_count_tables)
            skip_count_tables="true"
            shift 1
            ;;
        --skip_protein_read_counts)
            skip_protein_read_counts="true"
            shift 1
            ;;
        --skip_ko_read_counts)
            skip_ko_read_counts="true"
            shift 1
            ;;
        --skip_eggnog_ogs_read_counts)
            skip_eggnog_ogs_read_counts="true"
            shift 1
            ;;
        --skip_cazy_read_counts)
            skip_cazy_read_counts="true"
            shift 1
            ;;
        --skip_pfam_read_counts)
            skip_pfam_read_counts="true"
            shift 1
            ;;
        --compress_outputs)
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

# Stop the script if any required argument is missing, then print the help message
if [[ -z "$reads1" || -z "$reads2" || -z "$prefix" || -z "$outdir" ]]; then
    echo "ERROR: Missing required arguments:  --reads1, --reads2, --prefix and/or --outdir." >&2
    Help
    exit 1
fi
# Check that the first input read file exists before starting the pipeline
if [[ ! -f "$reads1" ]]; then
    echo "ERROR: reads1 file not found: $reads1" >&2
    exit 1
fi

if [[ ! -f "$reads2" ]]; then
    echo "ERROR: reads2 file not found: $reads2" >&2
    exit 1
fi

if [[ ! -f "$diamond_db" ]]; then
    echo "ERROR: DIAMOND database file not found: $diamond_db" >&2
    exit 1
fi

if [[ "${diamond_db}" != *.dmnd ]]; then
    echo "ERROR: DIAMOND database file must have .dmnd extension: $diamond_db" >&2
    exit 1
fi

if [[ ! -f "$eggnog_db" ]]; then
    echo "ERROR: EggNOG database file not found: $eggnog_db" >&2
    exit 1
fi

if [[ "${eggnog_db}" != *.emapper.annotations ]]; then
    echo "ERROR: EggNOG database file must have .emapper.annotations extension: $eggnog_db" >&2
    exit 1
fi
# Check that threads is a positive integer
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

if [[ ! -f "$script_dir/MycoFunk_count_input_reads.py" ]]; then
    echo "ERROR: Input read counting script not found: $script_dir/MycoFunk_count_input_reads.py" >&2
    exit 1
fi

if [[ ! -f "$script_dir/MycoFunk_annotate.py" ]]; then
    echo "ERROR: Annotation script not found: $script_dir/MycoFunk_annotate.py" >&2
    exit 1
fi


########################################################################################################
# Prepare output directories and log file
########################################################################################################

# Create output directories if they do not exist
mkdir -p "$outdir"
mkdir -p "$outdir/diamond"
mkdir -p "$outdir/annotation"

# Define run metadata log file
run_info_log="${outdir}/MycoFunk_${prefix}_run_info.log"

# Write run metadata to log file
cat > "$run_info_log" <<EOF
Tool: ${TOOL_NAME}
Version: ${VERSION}
Run date: $(date '+%Y-%m-%d %H:%M:%S')
Working directory: $(pwd)
Command: ${COMMAND_RUN}

Sample information:
  prefix: ${prefix}
  outdir: ${outdir}
  reads1: ${reads1}
  reads2: ${reads2}

Databases:
  diamond_db: ${diamond_db}
  eggnog_db: ${eggnog_db}

Run options:
  threads: ${threads}
  write_full_annotation: ${write_full_annotation}
  Skip unannotated hits table: ${skip_unannotated_hits}
  Skip all count tables: ${skip_count_tables}
  Skip protein read counts table: ${skip_protein_read_counts}
  Skip KO read counts table: ${skip_ko_read_counts}
  Skip eggNOG OG read counts table: ${skip_eggnog_ogs_read_counts}
  Skip CAZy read counts table: ${skip_cazy_read_counts}
  Skip PFAM read counts table: ${skip_pfam_read_counts}
  compress_outputs: ${compress_outputs}
  timestamps: ${timestamps}

Filtering:
  evalue filter: ${evalue:-DIAMOND default = 0.001}
  min-score filter: ${min_score:-not used}
  identity filter: ${id:-not used}
  query-cover filter: ${query_cover:-not used}
  top filter: ${top:-not used}

EOF


########################################################################################################
# Print run start messages
########################################################################################################

echo "###################################################################################################"
echo " _____                  _ _____            _     "
echo "|  ___|   _ _ __   __ _(_)  ___|   _ _ __ | | __ "
echo "| |_ | | | | '_ \ / _  | | |_ | | | | '_ \| |/ / "
echo "|  _|| |_| | | | | (_| | |  _|| |_| | | | |   <  "
echo "|_|   \__,_|_| |_|\__, |_|_|   \__,_|_| |_|_|\_\ "
echo "                  |___/                          "
echo ""
echo "Running ${TOOL_NAME} v${VERSION}"
echo "Sample: ${prefix}"
echo "Output directory: ${outdir}"
echo "Run information log: ${run_info_log}"
echo ""


########################################################################################################
# Create flags and additional functions
########################################################################################################

# Create a flag for "--evalue" option and check if it is active
evalue_flag=""
if [[ -n "$evalue" ]]; then
    evalue_flag="--evalue $evalue"
fi

# Create a flag for "--min-score" option and check if it is active
min_score_flag=""
if [[ -n "$min_score" ]]; then
    min_score_flag="--min-score $min_score"
fi

# Create a flag for "--id" option and check if it is active
id_flag=""
if [[ -n "$id" ]]; then
    id_flag="--id $id"
fi

# Create a flag for "--query-cover" option and check if it is active
query_cover_flag=""
if [[ -n "$query_cover" ]]; then
    query_cover_flag="--query-cover $query_cover"
fi

# Create a flag for "--top" option and check if it is active
top_flag=""
if [[ -n "$top" ]]; then
    top_flag="--top $top"
fi

log_msg() {
    if [[ "$timestamps" == "true" ]]; then
        echo "[$(date '+%H:%M:%S')] $*"
    else
        echo "$*"
    fi
}


########################################################################################################
# Run count of input reads script
########################################################################################################

echo "####################################### CHECKING INPUT FILES ######################################"

# Count input reads before running DIAMOND
log_msg "Counting input reads..."
log_msg "Writing input read counts table.."

python "$script_dir/MycoFunk_count_input_reads.py" \
    --reads1 "$reads1" \
    --reads2 "$reads2" \
    --prefix "$prefix" \
    --outdir "$outdir"

# Stop if input read counting failed
if [[ $? -ne 0 ]]; then
    echo "ERROR: Input read counting failed" >&2
    exit 1
fi

# Tell the user where the input read counts were written
log_msg "    (Written to: ${outdir}/${prefix}_input_read_counts.tsv)"


########################################################################################################
# Run DIAMOND on R1
########################################################################################################

echo "######################################### RUNNING DIAMOND #########################################"

log_msg "Running DIAMOND blastx on R1..."

# Run DIAMOND blastx on reads1 to align DNA reads against the fungal protein database
diamond blastx \
    $evalue_flag \
    $min_score_flag \
    $id_flag \
    $query_cover_flag \
    $top_flag \
    --outfmt 6 \
    --threads "$threads" \
    -d "$diamond_db" \
    -q "$reads1" \
    -o "$outdir/diamond/${prefix}_R1.diamond.out" \
    2> "$outdir/diamond/${prefix}_R1.diamond.log"


########################################################################################################
# Run DIAMOND on R2
########################################################################################################

log_msg "Running DIAMOND blastx on R2..."

# Run DIAMOND blastx on reads2 to align DNA reads against the fungal protein database
diamond blastx \
    $evalue_flag \
    $min_score_flag \
    $id_flag \
    $query_cover_flag \
    $top_flag \
    --outfmt 6 \
    --threads "$threads" \
    -d "$diamond_db" \
    -q "$reads2" \
    -o "$outdir/diamond/${prefix}_R2.diamond.out" \
    2> "$outdir/diamond/${prefix}_R2.diamond.log"


########################################################################################################
# Run annotation script
########################################################################################################

echo "################################## RUNNING FUNCTIONAL ANNOTATION ##################################"


# Create a flag for "--write_full_annotation" option and check if it is active
full_annotation_flag=""
if [[ "$write_full_annotation" == "true" ]]; then
    full_annotation_flag="--write_full_annotation"
fi

# Create a flag for "--skip_unannotated_hits" option and check if it is active
skip_unannotated_hits_flag=""
if [[ "$skip_unannotated_hits" == "true" ]]; then
    skip_unannotated_hits_flag="--skip_unannotated_hits"
fi

# Create a flag for "--skip_count_tables" option and check if it is active
skip_count_tables_flag=""
if [[ "$skip_count_tables" == "true" ]]; then
    skip_count_tables_flag="--skip_count_tables"
fi

# Create a flag for "--skip_protein_read_counts" option and check if it is active
skip_protein_read_counts_flag=""
if [[ "$skip_protein_read_counts" == "true" ]]; then
    skip_protein_read_counts_flag="--skip_protein_read_counts"
fi

# Create a flag for "--skip_ko_read_counts" option and check if it is active
skip_ko_read_counts_flag=""
if [[ "$skip_ko_read_counts" == "true" ]]; then
    skip_ko_read_counts_flag="--skip_ko_read_counts"
fi

# Create a flag for "--skip_eggnog_ogs_read_counts" option and check if it is active
skip_eggnog_ogs_read_counts_flag=""
if [[ "$skip_eggnog_ogs_read_counts" == "true" ]]; then
    skip_eggnog_ogs_read_counts_flag="--skip_eggnog_ogs_read_counts"
fi

# Create a flag for "--skip_cazy_read_counts" option and check if it is active
skip_cazy_read_counts_flag=""
if [[ "$skip_cazy_read_counts" == "true" ]]; then
    skip_cazy_read_counts_flag="--skip_cazy_read_counts"
fi

# Create a flag for "--skip_pfam_read_counts" option and check if it is active
skip_pfam_read_counts_flag=""
if [[ "$skip_pfam_read_counts" == "true" ]]; then
    skip_pfam_read_counts_flag="--skip_pfam_read_counts"
fi

# Create a flag for "--compress_outputs" option and check if it is active
compress_outputs_flag=""
if [[ "$compress_outputs" == "true" ]]; then
    compress_outputs_flag="--compress_outputs"
fi

# Create a flag for "--timestamps" option and check if it is active
timestamps_flag=""
if [[ "$timestamps" == "true" ]]; then
    timestamps_flag="--timestamps"
fi


# Pass both DIAMOND outputs to the annotation python script to assign functions and summarize pathways
python "$script_dir/MycoFunk_annotate.py" \
    --r1 "$outdir/diamond/${prefix}_R1.diamond.out" \
    --r2 "$outdir/diamond/${prefix}_R2.diamond.out" \
    --eggnog_db "$eggnog_db" \
    --prefix "$prefix" \
    --outdir "$outdir/annotation" \
    $full_annotation_flag \
    $skip_unannotated_hits_flag \
    $skip_count_tables_flag \
    $skip_protein_read_counts_flag \
    $skip_ko_read_counts_flag \
    $skip_eggnog_ogs_read_counts_flag \
    $skip_cazy_read_counts_flag \
    $skip_pfam_read_counts_flag \
    $compress_outputs_flag \
    $timestamps_flag

########################################################################################################
# Print a final message when all steps finish successfully
########################################################################################################

echo "############### ${TOOL_NAME} v${VERSION} completed successfully for sample: $prefix ###############"