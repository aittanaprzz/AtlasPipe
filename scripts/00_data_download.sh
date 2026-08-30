#!/bin/bash

# SCRIPT DE DESCARGA AUTOMÁTICA Y DESCOMPRESIÓN
# Este script descarga los archivos de MicrobeAtlas necesarios para ejecutar el
# pipeline. Los enlaces oficiales están integrados aquí para que el repositorio
# pueda descargarse de GitHub y ejecutarse sin un archivo auxiliar de links.
#
# Si un usuario avanzado quiere usar una lista distinta, puede pasar un TSV con:
#   MICROBEATLAS_LINKS_FILE=/ruta/links.tsv bash 00_data_download.sh
# El TSV debe tener de 2 a 4 columnas separadas por tabulador:
#   URL  archivo_descargado  archivo_final  sha256_esperado
# Las columnas 3 y 4 son opcionales. Si archivo_final no se indica y el archivo
# descargado termina en .gz, se espera el mismo nombre sin .gz tras descomprimir.
# Si sha256_esperado se indica, se compara contra el archivo descargado (huella digital del doc).
# Esa opción es opcional; el pipeline no depende de ella.

# CHECKSUM guarda el valor SHA-256 esperado, comprobaciones de la integridad de docs


# 0. Definición de variables y rutas
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${MICROBEATLAS_ROOT:-$(cd "$SCRIPT_DIR/.." && pwd)}"
DATA_DIR="${MICROBEATLAS_DATA_DIR:-$PROJECT_ROOT/data/intermediate}"

FILES_FILE="${MICROBEATLAS_LINKS_FILE:-}"
OUT_DIR="${MICROBEATLAS_RAW_DATA_DIR:-$PROJECT_ROOT/data/raw/microbeatlas}"
MANIFEST_FILE="${MICROBEATLAS_DOWNLOAD_MANIFEST:-$OUT_DIR/download_manifest.tsv}"
SHA256_FILE="${MICROBEATLAS_DOWNLOAD_SHA256:-$OUT_DIR/download_sha256.tsv}"

mkdir -p "$OUT_DIR"

downloaded_files=()     #comprobaciones de descarga
final_files=()

is_missing_checksum() {
    local CHECKSUM="$1"
    [ -z "$CHECKSUM" ] || [ "$CHECKSUM" = "-" ] || [ "$CHECKSUM" = "NA" ]
}

sha256_of() {
    sha256sum "$1" | awk '{ print $1 }'
}

verify_downloaded_file() {
    local PATH_TO_FILE="$1"
    local EXPECTED_SHA256="$2"

    if [ ! -s "$PATH_TO_FILE" ]; then
        echo "Archivo descargado vacío o inexistente: $PATH_TO_FILE"
        exit 1
    fi

    if [[ "$PATH_TO_FILE" == *.gz ]]; then
        gzip -t "$PATH_TO_FILE"
    fi

    local OBSERVED_SHA256
    OBSERVED_SHA256="$(sha256_of "$PATH_TO_FILE")"

    if ! is_missing_checksum "$EXPECTED_SHA256"; then
        if [ "$OBSERVED_SHA256" != "$EXPECTED_SHA256" ]; then
            echo "Checksum SHA256 incorrecto para $PATH_TO_FILE"
            echo "  esperado:  $EXPECTED_SHA256"
            echo "  observado: $OBSERVED_SHA256"
            exit 1
        fi
    fi
}

download_one() {
    local URL="$1"
    local FILE_NAME="$2"
    local FINAL_NAME="${3:-}"
    local EXPECTED_SHA256="${4:-}"

    if [[ "$URL" == \#* || -z "$URL" ]]; then
        return 0
    fi

    if [ -z "$FILE_NAME" ]; then
        echo "No se indicó nombre de archivo para la URL: $URL"
        exit 1
    fi

    if [ -z "$FINAL_NAME" ]; then
        if [[ "$FILE_NAME" == *.gz ]]; then
            FINAL_NAME="${FILE_NAME%.gz}"
        else
            FINAL_NAME="$FILE_NAME"
        fi
    fi

    local DOWNLOAD_PATH="$OUT_DIR/$FILE_NAME"
    local FINAL_PATH="$OUT_DIR/$FINAL_NAME"
    local TMP_PATH="$DOWNLOAD_PATH.part"

    downloaded_files+=("$DOWNLOAD_PATH")       #comprobaciones de descarga
    final_files+=("$FINAL_PATH")

    if [ -s "$DOWNLOAD_PATH" ]; then
        echo "Archivo ya descargado, verificando: $FILE_NAME"
        verify_downloaded_file "$DOWNLOAD_PATH" "$EXPECTED_SHA256"
    else
        echo "Descargando archivo $FILE_NAME:"
        wget "$URL" -O "$TMP_PATH" -c
        verify_downloaded_file "$TMP_PATH" "$EXPECTED_SHA256"
        mv "$TMP_PATH" "$DOWNLOAD_PATH"
    fi
}

download_from_stream() {
    while IFS=$'\t' read -r URL FILE_NAME FINAL_NAME EXPECTED_SHA256; do
        download_one "$URL" "$FILE_NAME" "$FINAL_NAME" "$EXPECTED_SHA256"
    done
}

# 1. Descargar archivos. Por defecto se usa la lista integrada. Si se define
# MICROBEATLAS_LINKS_FILE, se usa esa lista personalizada.
if [ -n "$FILES_FILE" ]; then
    if [ ! -f "$FILES_FILE" ]; then
        echo "No existe MICROBEATLAS_LINKS_FILE: $FILES_FILE"
        exit 1
    fi

    download_from_stream < "$FILES_FILE"
else
    download_from_stream <<'MICROBEATLAS_DOWNLOADS'
# URL	download_name	final_name	sha256
https://microbeatlas.org/downloads/paper/community_cluster_metadata.tsv	community_cluster_metadata.tsv	community_cluster_metadata.tsv	-
https://microbeatlas.org/downloads/paper/community_cluster_mapping.tsv	community_cluster_mapping.tsv	community_cluster_mapping.tsv	-
https://microbeatlas.org/downloads/otutables/samples-otus.97.mapped.biom.gz	samples-otus.97.mapped.biom.gz	samples-otus.97.mapped.biom	-
https://microbeatlas.org/downloads/otus/otus.97.allinfo	otus.97.allinfo	otus.97.allinfo	-
https://microbeatlas.org/downloads/mapref/otus.info.tsv	otus.info.tsv	otus.info.tsv	-
https://microbeatlas.org/downloads/samples/samples.env.info.tsv	samples.env.info.tsv	samples.env.info.tsv	-
MICROBEATLAS_DOWNLOADS
fi

echo "Descarga completada."

# 4. Descompresión de archivos
# 4.1. comprobaciones previas
if [ ! -d "$OUT_DIR" ]; then
    echo "No existe el directorio $OUT_DIR" 
    exit 1
fi

# 4.2. Descomprimir solo los .gz gestionados por este script
echo "Descomprimiendo archivos .gz descargados por este script"
for DOWNLOAD_PATH in "${downloaded_files[@]}"; do
    if [[ "$DOWNLOAD_PATH" == *.gz ]]; then
        gzip -dkf "$DOWNLOAD_PATH"
    fi
done

# 5. Comprobar archivos finales esperados y guardar checksums observados
{
    printf "file\tsha256\n"
    for DOWNLOAD_PATH in "${downloaded_files[@]}"; do
        if [ -f "$DOWNLOAD_PATH" ]; then
            printf "%s\t%s\n" "$DOWNLOAD_PATH" "$(sha256_of "$DOWNLOAD_PATH")"
        fi
    done
    for FINAL_PATH in "${final_files[@]}"; do
        if [ -f "$FINAL_PATH" ]; then
            printf "%s\t%s\n" "$FINAL_PATH" "$(sha256_of "$FINAL_PATH")"
        fi
    done
} > "$SHA256_FILE"

{
    printf "parameter\tvalue\n"
    printf "output_dir\t%s\n" "$OUT_DIR"
    printf "sha256_file\t%s\n" "$SHA256_FILE"
    for FINAL_PATH in "${final_files[@]}"; do
        if [ ! -s "$FINAL_PATH" ]; then
            echo "Falta archivo final esperado o está vacío: $FINAL_PATH"
            exit 1
        fi
        printf "final_file\t%s\n" "$FINAL_PATH"
    done
} > "$MANIFEST_FILE"

echo "Archivos finales verificados:"
cat "$MANIFEST_FILE"
echo "Checksums SHA256 guardados en: $SHA256_FILE"
