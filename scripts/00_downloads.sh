# AUTOMATED DOWNLOAD AND VALIDATION SCRIPT
#
# This script downloads the MicrobeAtlas input files required by the pipeline.
# The official download URLs are embedded in the script, so the repository can
# be cloned and executed without requiring an external links file.
#
# Basic usage:
#   bash 00_downloads.sh
#
# Optional configuration:
#   ATLASPIPE_ROOT=/path/to/project bash 00_downloads.sh
#   MICROBEATLAS_RAW_DATA_DIR=/path/to/raw_data bash 00_downloads.sh


# 0. Paths and configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${ATLASPIPE_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
DATA_DIR="${ATLASPIPE_DATA_DIR:-$PROJECT_ROOT/data/intermediate}"

OUT_DIR="${MICROBEATLAS_RAW_DATA_DIR:-$PROJECT_ROOT/data/raw}"
MANIFEST_FILE="${ATLASPIPE_DOWNLOAD_MANIFEST:-$OUT_DIR/download_manifest.tsv}"

mkdir -p "$OUT_DIR"

downloaded_files=()
final_files=()

verify_downloaded_file() {
    local PATH_TO_FILE="$1"

    if [ ! -s "$PATH_TO_FILE" ]; then
        echo "Downloaded file is empty or missing: $PATH_TO_FILE"
        exit 1
    fi

    if [[ "$PATH_TO_FILE" == *.gz ]]; then
        gzip -t "$PATH_TO_FILE"
    fi
}

download_one() {
    local URL="$1"
    local FINAL_NAME="$2"

    if [[ "$URL" == \#* || -z "$URL" ]]; then
        return 0
    fi

    if [ -z "$FINAL_NAME" ]; then
        echo "No final file name was provided for URL: $URL"
        exit 1
    fi

    local DOWNLOAD_NAME="$FINAL_NAME"
    if [[ "$URL" == *.gz && "$FINAL_NAME" != *.gz ]]; then
        DOWNLOAD_NAME="$FINAL_NAME.gz"
    fi

    local DOWNLOAD_PATH="$OUT_DIR/$DOWNLOAD_NAME"
    local FINAL_PATH="$OUT_DIR/$FINAL_NAME"
    local TMP_PATH="$DOWNLOAD_PATH.part"

    downloaded_files+=("$DOWNLOAD_PATH")
    final_files+=("$FINAL_PATH")

    if [ -s "$DOWNLOAD_PATH" ]; then
        echo "File already downloaded, validating: $DOWNLOAD_NAME"
        verify_downloaded_file "$DOWNLOAD_PATH"
    else
        echo "Downloading file $DOWNLOAD_NAME:"
        wget "$URL" -O "$TMP_PATH" -c
        verify_downloaded_file "$TMP_PATH"
        mv "$TMP_PATH" "$DOWNLOAD_PATH"
    fi
}

download_from_stream() {
    while IFS=$'\t' read -r URL FINAL_NAME; do
        download_one "$URL" "$FINAL_NAME"
    done
}

# 1. Download MicrobeAtlas files

download_from_stream <<'MICROBEATLAS_DOWNLOADS'
# URL	final_name
https://microbeatlas.org/downloads/paper/community_cluster_metadata.tsv	community_cluster_metadata.tsv
https://microbeatlas.org/downloads/paper/community_cluster_mapping.tsv	community_cluster_mapping.tsv
https://microbeatlas.org/downloads/paper/metag_minfilter/samples-otus.97.metag.minfilter.rod2025companion.biom.gz	samples-otus.97.metag.minfilter.rod2025companion.biom
https://microbeatlas.org/downloads/otus/otus.97.allinfo	otus.97.allinfo
MICROBEATLAS_DOWNLOADS

if [ $? -ne 0 ]; then
    echo "Error: MicrobeAtlas files could not be downloaded."
    exit 1
fi

echo "Download completed."

# 2. Decompress downloaded files
# 2.1. Pre-checks
if [ ! -d "$OUT_DIR" ]; then
    echo "Directory does not exist: $OUT_DIR"
    exit 1
fi

# 2.2. Decompress only the .gz files managed by this script
echo "Decompressing .gz files downloaded by this script"
for DOWNLOAD_PATH in "${downloaded_files[@]}"; do
    if [[ "$DOWNLOAD_PATH" == *.gz ]]; then
        gzip -dkf "$DOWNLOAD_PATH"
    fi
done

# 3. Validate final files

{
    printf "parameter\tvalue\n"
    printf "output_dir\t%s\n" "$OUT_DIR"
} > "$MANIFEST_FILE"

for FINAL_PATH in "${final_files[@]}"; do

    if [ ! -s "$FINAL_PATH" ]; then
        echo "Error: final file is missing or empty:"
        echo "  $FINAL_PATH"
        exit 1
    fi

    printf "final_file\t%s\n" "$FINAL_PATH" >> "$MANIFEST_FILE"
done

echo "File validation completed."
echo "Verified final files:"
cat "$MANIFEST_FILE"
