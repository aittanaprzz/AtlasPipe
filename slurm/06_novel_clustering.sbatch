#!/bin/bash
  
#SBATCH --job-name=novel_clustering
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=80G
#SBATCH --time=24:00:00
#SBATCH --qos=medium

set -euo pipefail

source "${ATLASPIPE_ROOT:-${SLURM_SUBMIT_DIR:?}}/slurm/common.sh"
atlaspipe_init_slurm_step "06" "novel_clustering"

echo "Iniciando clustering de los amplicones no mapeados con VSEARCH en Garnatxa..."

ml anaconda
conda activate MicrobeAtlas

SCRIPTS_DIR="$PROJECT_ROOT/scripts"
DATA_DIR="${ATLASPIPE_DATA_DIR:-$PROJECT_ROOT/data/intermediate}"
ENV_LABEL="${MICROBEATLAS_ENV_LABEL:-aquatic}"
SAMPLES_PER_CC="${MICROBEATLAS_SAMPLES_PER_CC:-1}"

selection_label_from_samples_per_cc() {
    case "$1" in
        all|todos)
            printf "allsamplesperCC"
            ;;
        1)
            printf "1sampelperCC"
            ;;
        ''|*[!0-9]*)
            echo "MICROBEATLAS_SAMPLES_PER_CC debe ser un entero positivo o 'all': $1" >&2
            exit 1
            ;;
        *)
            if [ "$1" -lt 1 ]; then
                echo "MICROBEATLAS_SAMPLES_PER_CC debe ser un entero positivo o 'all': $1" >&2
                exit 1
            fi
            printf "%ssamplesperCC" "$1"
            ;;
    esac
}
SELECTION_LABEL="${MICROBEATLAS_SELECTION_LABEL:-$(selection_label_from_samples_per_cc "$SAMPLES_PER_CC")}"
if [ "$SELECTION_LABEL" = "1sampelperCC" ]; then
    VSEARCH_SUFFIX=""
else
    VSEARCH_SUFFIX="_${ENV_LABEL}_${SELECTION_LABEL}"
fi

IN_DIR="${ATLASPIPE_OWN_AMPLICON_MAPPING_DIR:-$DATA_DIR/05_amplicon_mapping}"
AMPLICON_FILE="${ATLASPIPE_OWN_AMPLICON_FASTA:-$IN_DIR/amplicons.fasta}"
VSEARCH_FILE="${ATLASPIPE_OWN_AMPLICON_OTU97_MAP:-$IN_DIR/OTUs97_vsearch${VSEARCH_SUFFIX}.blast6}"
NO_HITS_FILE="${ATLASPIPE_OWN_NO_HITS_FASTA:-$IN_DIR/amplicons_no_hits${VSEARCH_SUFFIX}.fasta}"
UC_FILE="${ATLASPIPE_OWN_NO_HITS_UC:-$IN_DIR/amplicon_clusters97${VSEARCH_SUFFIX}.uc}"
CENTROIDS_FILE="${ATLASPIPE_OWN_NO_HITS_CENTROIDS:-$IN_DIR/cluster_centroids97${VSEARCH_SUFFIX}.fasta}"
FEATURE_MAP_FILE="${ATLASPIPE_OWN_NO_HITS_FEATURE_MAP:-$IN_DIR/amplicon_cluster_map${VSEARCH_SUFFIX}.tsv}"
CLUSTER_SIZES_FILE="${ATLASPIPE_OWN_NO_HITS_CLUSTER_SIZES:-$IN_DIR/cluster_sizes${VSEARCH_SUFFIX}.tsv}"
OUT_BLAST6OUT="${ATLASPIPE_OWN_NO_HITS_BLAST6:-$IN_DIR/amplicon_clusters97${VSEARCH_SUFFIX}.blast6}"
SUMMARY_FILE="${ATLASPIPE_OWN_NO_HITS_CLUSTER_SUMMARY:-$IN_DIR/cluster_summary97${VSEARCH_SUFFIX}.tsv}"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

cd "$IN_DIR"

if [ ! -f "$AMPLICON_FILE" ]; then
    echo "No existe el FASTA de amplicones de entrada: $AMPLICON_FILE"
    exit 1
fi

if [ ! -f "$VSEARCH_FILE" ]; then
    echo "No existe la salida de mapeo del paso 05: $VSEARCH_FILE"
    echo "Error en el paso 05 de mapeo de amplicones a OTUs 97%"
    exit 1
fi

awk '
  NR == FNR {
    if (NF >= 2 && $2 != "*") {
      hit[$1] = 1
    }
    next
  }
  /^>/ {
    id = substr($0, 2)
    sub(/[ \t].*/, "", id)
    keep = !(id in hit)
  }
  keep {
    print
  }
' "$VSEARCH_FILE" "$AMPLICON_FILE" > "$NO_HITS_FILE"

NO_HITS_COUNT=$(awk 'BEGIN { n = 0 } /^>/ { n++ } END { print n }' "$NO_HITS_FILE")
echo "Amplicones no mapeados extraidos: $NO_HITS_COUNT"

if [ "$NO_HITS_COUNT" -eq 0 ]; then
    echo "No hay amplicones no mapeados para clusterizar."
    conda activate Renv
    Rscript "$SCRIPTS_DIR/06_novel_clustering_parse.R" \
      --uc "$UC_FILE" \
      --nohit-fasta "$NO_HITS_FILE" \
      --feature-map "$FEATURE_MAP_FILE" \
      --cluster-sizes "$CLUSTER_SIZES_FILE" \
      --summary "$SUMMARY_FILE"
    exit 0
fi

rm -f "$UC_FILE" "$CENTROIDS_FILE" "$OUT_BLAST6OUT"

vsearch --cluster_fast "$NO_HITS_FILE" \
  --id 0.97 \
  --blast6out "$OUT_BLAST6OUT" \
  --centroids "$CENTROIDS_FILE" \
  --uc "$UC_FILE" \
  --strand both \
  --threads "$THREADS"

conda activate Renv
Rscript "$SCRIPTS_DIR/06_novel_clustering_parse.R" \
  --uc "$UC_FILE" \
  --nohit-fasta "$NO_HITS_FILE" \
  --feature-map "$FEATURE_MAP_FILE" \
  --cluster-sizes "$CLUSTER_SIZES_FILE" \
  --summary "$SUMMARY_FILE"

UC_SEEDS=$(awk 'BEGIN { n = 0 } $1 == "S" { n++ } END { print n }' "$UC_FILE")
UC_HITS=$(awk 'BEGIN { n = 0 } $1 == "H" { n++ } END { print n }' "$UC_FILE")

echo "Clusters centroides: $UC_SEEDS"
echo "Amplicones asignados a centroides: $UC_HITS"
cat "$SUMMARY_FILE"

echo "Clustering con VSEARCH finalizado:"
echo "  BLAST6: $OUT_BLAST6OUT"
echo "  Centroides: $CENTROIDS_FILE"
echo "  Clusters UC: $UC_FILE"
echo "  Mapa amplicon no-hit -> cluster: $FEATURE_MAP_FILE"
echo "  Tamanos de clusters: $CLUSTER_SIZES_FILE"
echo "  Resumen: $SUMMARY_FILE"
