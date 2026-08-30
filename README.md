# Pipeline de contextualización de ASVs con MicrobeAtlas

Pipeline bioinformática para contextualizar ASVs propios frente a OTUs al 97% de MicrobeAtlas. La pipeline combina datos públicos de MicrobeAtlas con una tabla ASV propia, asigna ASVs a OTUs de referencia mediante VSEARCH, agrupa los ASVs no mapeados y genera una tabla final de abundancias contextualizadas.

## Modelo de ejecución

Esta pipeline está pensada para ejecutarse en un servidor HPC con SLURM.

Los archivos `.sbatch` forman parte oficial de la pipeline. Algunos pasos también pueden ejecutarse en una sesión interactiva si el tamaño de los datos lo permite, pero los pasos pesados deben enviarse al clúster con `sbatch`.

La ejecución es secuencial: cada paso debe terminar correctamente antes de lanzar el siguiente, salvo que se usen dependencias de SLURM.

## Estructura del proyecto

```text
MicrobeAtlas-pipeline/
|-- README.md
|-- environment.yml
|-- config/
|   |-- example_env.sh
|   |-- env_aquatic_1sample.sh
|   `-- env.sh                  # local, ignorado por Git
|-- scripts/
|   |-- 00_data_download.sh
|   |-- 01_prepare_own_asv_inputs.R
|   |-- 02_filter_samples_1sampelperCC.R
|   |-- 04_multifasta.R
|   |-- 06_vsearch_clustering_parse.R
|   `-- 07_tabla_abundancias_profesor.R
|-- slurm/
|   |-- 00_data_download.sbatch
|   |-- 01_prepare_own_asv_inputs.sbatch
|   |-- 02_filter_samples_1sampelperCC.sbatch
|   |-- 03_extract_biome.sbatch
|   |-- 04_multifasta.sbatch
|   |-- 05_vsearch.sbatch
|   |-- 06_vsearch_clustering.sbatch
|   |-- 07_tabla_abundancias_profesor.sbatch
|   `-- common.sh
|-- docs/
|-- data/
|   |-- input/
|   |-- raw/
|   |   `-- microbeatlas/
|   `-- intermediate/
|-- results/
`-- log/
```

## Requisitos

- Conda o Mamba
- SLURM para la ejecución en HPC
- R
- VSEARCH
- Herramientas de formato BIOM
- `wget`, `gzip`, `awk`

## Entorno

Crea el entorno de Conda:

```bash
conda env create -f environment.yml
conda activate microbeatlas-pipeline
```

El archivo `environment.yml` documenta un entorno de referencia con el software principal requerido por la pipeline. Es útil para la reproducibilidad y para instalaciones locales o de prueba.

En el clúster Garnatxa, los archivos `.sbatch` oficiales activan actualmente varios entornos de Conda existentes en lugar de un único entorno unificado:

| Entorno | Uso |
|---|---|
| `Renv` | Scripts de R y dependencias de paquetes de R |
| `Biom` | Extracción/conversión de tablas BIOM |
| `MicrobeAtlas` | Pasos basados en VSEARCH |

Esta separación es intencional y específica de Garnatxa. Las herramientas se instalaron en entornos separados debido a conflictos de Python/versiones entre dependencias. Si la pipeline se despliega en otra máquina, crea nombres de entorno equivalentes o edita las líneas `conda activate ...` en `slurm/*.sbatch` para que coincidan con la instalación local.

## Configuración

El repositorio incluye `config/example_env.sh` como plantilla versionada. No lo edites para ejecuciones rutinarias. En su lugar, crea un archivo de configuración local:

```bash
cd /storage/microbiome/home_members/pesai/MicrobeAtlas-pipeline
cp config/example_env.sh config/env.sh
nano config/env.sh
```

Para ejecuciones simples de benchmark/prueba, mantén una configuración con nombre por cada prueba y copia la que quieras usar en `config/env.sh`. Por ejemplo:

```bash
cp config/env_aquatic_1sample.sh config/env.sh
```

Después, carga la configuración local antes de ejecutar pasos interactivos:

```bash
source config/env.sh
```

El archivo local `config/env.sh` está ignorado por Git y debe contener la ruta del proyecto específica de la máquina, la tabla ASV de entrada, la etiqueta del entorno y el número de muestras por clúster de comunidad. El script auxiliar `submit_pipeline_slurm.sh` también carga `config/env.sh` antes de enviar trabajos.

Cada banco de pruebas debe usar un `MICROBEATLAS_RUN_LABEL` corto. Los datos crudos descargados de MicrobeAtlas se comparten, pero las salidas intermedias y finales se escriben en directorios específicos de cada ejecución:

```text
data/raw/microbeatlas/
data/intermediate/<run_label>/
results/<run_label>/
log/<run_label>/
```

Para la configuración incluida `aquatic_1sample`, las salidas se escriben en:

```text
data/intermediate/aquatic_1sample/
results/aquatic_1sample/
log/aquatic_1sample/
```

Variables principales:

| Variable | Significado |
|---|---|
| `MICROBEATLAS_ROOT` | Directorio raíz del proyecto |
| `MICROBEATLAS_RUN_LABEL` | Nombre corto para el banco de pruebas/ejecución actual |
| `MICROBEATLAS_INPUT_DIR` | Directorio para archivos de entrada proporcionados por el usuario |
| `MICROBEATLAS_RAW_DATA_DIR` | Directorio para archivos descargados de MicrobeAtlas |
| `MICROBEATLAS_INTERMEDIATE_DIR` | Directorio para archivos intermedios de la pipeline |
| `MICROBEATLAS_RESULTS_DIR` | Directorio para archivos de resultados finales |
| `MICROBEATLAS_LOG_DIR` | Directorio para stdout/stderr de SLURM de la ejecución actual |
| `MICROBEATLAS_DATA_DIR` | Alias retrocompatible de `MICROBEATLAS_INTERMEDIATE_DIR` |
| `MICROBEATLAS_OWN_ASV_TABLE_RAW` | Tabla de entrada con ASVs del usuario |
| `MICROBEATLAS_ENV_LABEL` | Etiqueta usada en los nombres de salida |
| `MICROBEATLAS_SAMPLES_PER_CC` | Número de muestras seleccionadas por clúster de comunidad |

Los pasos enviados con `submit_pipeline_slurm.sh` o con `sbatch slurm/<paso>.sbatch`
guardan stdout y stderr dentro de `MICROBEATLAS_LOG_DIR`. Por defecto, cada
banco de pruebas usa una carpeta independiente bajo `log/<run_label>/`, con un
par de archivos por paso y job:

```text
log/<run_label>/02_filter_samples_<job_id>.out
log/<run_label>/02_filter_samples_<job_id>.err
```

## Tabla ASV de entrada

Coloca la tabla ASV del usuario en:

```text
data/input/final_table_filtered.tsv
```

La tabla ASV del usuario debe contener al menos:

| Columna | Descripción |
|---|---|
| `asv_code` | Identificador del ASV. Formato esperado: `asv_*` |
| `asv` | Secuencia nucleotídica del ASV |

Formato exacto esperado:

- La tabla debe tener una fila de cabecera.
- Cada fila debe representar un ASV.
- La primera columna de identificadores debe llamarse `asv_code`.
- La columna de secuencia debe llamarse `asv`.
- Los identificadores de `asv_code` deben empezar por `asv_`, por ejemplo `asv_0001`.
- Las secuencias de `asv` deben contener solo caracteres nucleotídicos válidos (`A`, `C`, `G`, `T`, `U`, `N` o `-`).
- Cada columna de muestra debe contener abundancias numéricas.
- No debe haber ASVs duplicados ni nombres de muestra duplicados.

Ejemplo mínimo:

```text
asv_code	asv	muestra_1	muestra_2	muestra_3
asv_0001	ACGTACGTACGT	10	0	3
asv_0002	TGCAACGTACGT	0	5	1
```

Todas las columnas que no sean metadatos se interpretan como columnas de abundancia de muestras. Las columnas de metadatos/taxonomía como `bim1`, `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus` y `Species` se ignoran para la extracción de abundancias. La matriz debe estar orientada con ASVs en filas y muestras en columnas; las tablas con ASVs en columnas y muestras en filas deben transponerse antes de usarlas como entrada.

## Contrato de la pipeline

| Paso | Script/SBATCH | Entrada principal | Salida principal | Entorno |
|---|---|---|---|---|
| 00 | `scripts/00_data_download.sh`, `slurm/00_data_download.sbatch` | URLs de descarga de MicrobeAtlas incluidas en el script | `data/raw/microbeatlas/` | herramientas de shell |
| 01 | `scripts/01_prepare_own_asv_inputs.R`, `slurm/01_prepare_own_asv_inputs.sbatch` | `data/input/final_table_filtered.tsv` | `$MICROBEATLAS_INTERMEDIATE_DIR/01_own_asv_inputs/ASV_abundances.tsv`; `$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/ASV.fasta` | `Renv` |
| 02 | `scripts/02_write_biom_sample_ids.py`, `scripts/02_filter_samples_1sampelperCC.R`, `slurm/02_filter_samples_1sampelperCC.sbatch` | `data/raw/microbeatlas/community_cluster_metadata.tsv`; `data/raw/microbeatlas/community_cluster_mapping.tsv`; `data/raw/microbeatlas/samples-otus.97.mapped.biom` | `$MICROBEATLAS_INTERMEDIATE_DIR/02_filters/biom_sample_ids.txt`; `$MICROBEATLAS_INTERMEDIATE_DIR/02_filters_<env>_<selection>/` | `python3` con `h5py`, después `Renv` |
| 03 | `slurm/03_extract_biome.sbatch` | `data/raw/microbeatlas/samples-otus.97.mapped.biom`; `$MICROBEATLAS_INTERMEDIATE_DIR/02_filters_<env>_<selection>/community_cluster_selection_<selection>_forbiomsubset.txt` | `$MICROBEATLAS_INTERMEDIATE_DIR/03_samples_subsets/` | `Biom` |
| 04 | `scripts/04_multifasta.R`, `slurm/04_multifasta.sbatch` | `$MICROBEATLAS_INTERMEDIATE_DIR/03_samples_subsets/samples-otus.97_selection_<env>_<selection>.tsv`; `data/raw/microbeatlas/otus.97.allinfo` | `$MICROBEATLAS_INTERMEDIATE_DIR/04_multifasta/` | `Renv` |
| 05 | `slurm/05_vsearch.sbatch` | `$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/ASV.fasta`; `$MICROBEATLAS_INTERMEDIATE_DIR/04_multifasta/MicrobeAtlas_selection_<env>_<selection>_otus.97.fasta` | `$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/OTUs97_vsearch*.blast6` | `MicrobeAtlas` |
| 06 | `slurm/06_vsearch_clustering.sbatch`, `scripts/06_vsearch_clustering_parse.R` | `$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/ASV.fasta`; `$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/OTUs97_vsearch*.blast6` | `$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/ASVs_no_hits*`; `$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/ASVs_no_hits_feature_map*.tsv` | `MicrobeAtlas`, después `Renv` |
| 07 | `scripts/07_tabla_abundancias_profesor.R`, `slurm/07_tabla_abundancias_profesor.sbatch` | Salidas de los pasos 01, 03, 05 y 06 | `$MICROBEATLAS_RESULTS_DIR/07_tabla_abundancias_profesor/tabla_abundancias_features_contextualizadas.tsv` | `Renv` |

## Ejecución paso a paso

Ejecuta los pasos en orden.

### 00. Descargar los datos de MicrobeAtlas

Recomendado: SLURM.

```bash
sbatch slurm/00_data_download.sbatch
```

Ejecución interactiva alternativa:

```bash
bash scripts/00_data_download.sh
```

Salidas principales:

```text
data/raw/microbeatlas/
```

### 01. Preparar las entradas ASV del usuario

Recomendado: interactivo o SLURM.

```bash
Rscript scripts/01_prepare_own_asv_inputs.R
```

O bien:

```bash
sbatch slurm/01_prepare_own_asv_inputs.sbatch
```

Salidas principales:

```text
$MICROBEATLAS_INTERMEDIATE_DIR/01_own_asv_inputs/ASV_abundances.tsv
$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/ASV.fasta
```

### 02. Filtrar muestras de MicrobeAtlas

Recomendado: interactivo o SLURM.

Listar entornos disponibles:

```bash
Rscript scripts/02_filter_samples_1sampelperCC.R --list-envs
```

Ejecutar el filtrado:

```bash
python3 scripts/02_write_biom_sample_ids.py \
  --biom "$MICROBEATLAS_RAW_DATA_DIR/samples-otus.97.mapped.biom" \
  --output "$MICROBEATLAS_INTERMEDIATE_DIR/02_filters/biom_sample_ids.txt"

Rscript scripts/02_filter_samples_1sampelperCC.R \
  --env-prefix "$MICROBEATLAS_ENV_LABEL" \
  --samples-per-cc "$MICROBEATLAS_SAMPLES_PER_CC"
```

O bien con SLURM:

```bash
sbatch slurm/02_filter_samples_1sampelperCC.sbatch
```

Salidas principales:

```text
$MICROBEATLAS_INTERMEDIATE_DIR/02_filters/biom_sample_ids.txt
$MICROBEATLAS_INTERMEDIATE_DIR/02_filters_<env>_<selection>/
```

El envoltorio SLURM del paso 02 genera automáticamente `biom_sample_ids.txt` a partir de la tabla BIOM antes de seleccionar muestras. Esto evita muestrear registros que existen en los metadatos de MicrobeAtlas pero están ausentes en la tabla BIOM de abundancias usada por el paso 03.

### 03. Extraer el subconjunto BIOM

Recomendado: SLURM.

```bash
sbatch slurm/03_extract_biome.sbatch
```

Salidas principales:

```text
$MICROBEATLAS_INTERMEDIATE_DIR/03_samples_subsets/
```

### 04. Generar el multifasta de OTUs de MicrobeAtlas

Recomendado: interactivo o SLURM si el tamaño de los datos es alto.

```bash
Rscript scripts/04_multifasta.R
```

O bien:

```bash
sbatch slurm/04_multifasta.sbatch
```

Salidas principales:

```text
$MICROBEATLAS_INTERMEDIATE_DIR/04_multifasta/
```

### 05. Mapear ASVs contra OTUs de MicrobeAtlas con VSEARCH

Recomendado: SLURM.

```bash
sbatch slurm/05_vsearch.sbatch
```

Salidas principales:

```text
$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/OTUs97_vsearch*.blast6
```

### 06. Agrupar ASVs no mapeados

Recomendado: SLURM.

```bash
sbatch slurm/06_vsearch_clustering.sbatch
```

Este paso de SLURM extrae los ASVs sin aciertos contra OTUs de MicrobeAtlas, agrupa esos ASVs no mapeados con VSEARCH y después ejecuta `scripts/06_vsearch_clustering_parse.R` para construir el mapa de características sin acierto y el resumen de tamaños de clúster.

Salidas principales:

```text
$MICROBEATLAS_INTERMEDIATE_DIR/05_vsearch/ASVs_no_hits*
```

### 07. Construir la tabla final de abundancias contextualizadas

Recomendado: interactivo o SLURM.

```bash
Rscript scripts/07_tabla_abundancias_profesor.R
```

O bien:

```bash
sbatch slurm/07_tabla_abundancias_profesor.sbatch
```

Salidas principales:

```text
$MICROBEATLAS_RESULTS_DIR/07_tabla_abundancias_profesor/
```

La tabla final principal es:

```text
$MICROBEATLAS_RESULTS_DIR/07_tabla_abundancias_profesor/tabla_abundancias_features_contextualizadas.tsv
```

## Envío opcional a SLURM con dependencias

El script auxiliar `submit_pipeline_slurm.sh` envía los pasos disponibles de SLURM usando dependencias `afterok`, de modo que cada trabajo empieza solo si el anterior termina correctamente.

```bash
bash submit_pipeline_slurm.sh
```

Revisa los trabajos enviados con:

```bash
squeue -u "$USER"
```

## Notas importantes

No lances todos los trabajos de SLURM a la vez salvo que las dependencias estén configuradas.

La pipeline escribe archivos intermedios grandes en `data/`. Estos archivos están ignorados por Git.

Esta versión se distribuye intencionalmente como una pipeline secuencial documentada mediante scripts. Una versión futura puede migrar los mismos pasos a un gestor de flujos de trabajo como Snakemake o Nextflow.
