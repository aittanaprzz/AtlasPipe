# AtlasPipe

Pipeline bioinformática para contextualizar amplicones propios frente a OTUs al 97% de MicrobeAtlas. AtlasPipe combina datos públicos de MicrobeAtlas con una tabla de amplicones propia, asigna amplicones a OTUs de referencia mediante VSEARCH, agrupa los amplicones no mapeados y genera una tabla final de abundancias contextualizadas.

## Modelo de ejecución

Esta pipeline está pensada para ejecutarse en un servidor HPC con SLURM.

Los archivos `.sbatch` forman parte oficial de la pipeline. Algunos pasos también pueden ejecutarse en una sesión interactiva si el tamaño de los datos lo permite, pero los pasos pesados deben enviarse al clúster con `sbatch`.

La ejecución es secuencial: cada paso debe terminar correctamente antes de lanzar el siguiente, salvo que se usen dependencias de SLURM.

## Estructura del proyecto

```text
AtlasPipe/
|-- README.md
|-- environment.yml
|-- config/
|   |-- example_env.sh
|   |-- env_aquatic_1sample.sh
|   `-- env.sh                  # local, ignorado por Git
|-- scripts/
|   |-- 00_downloads.sh
|   |-- 01_prepare_inputs.R
|   |-- 02_sample_selection.R
|   |-- 03_biom_subset.sh
|   |-- 04_reference_fasta.R
|   |-- 05_amplicon_mapping.sh
|   |-- 06_novel_clustering_parse.R
|   `-- 07_final_table.R
|-- slurm/
|   |-- 00_downloads.sbatch
|   |-- 01_prepare_inputs.sbatch
|   |-- 02_sample_selection.sbatch
|   |-- 03_biom_subset.sbatch
|   |-- 04_reference_fasta.sbatch
|   |-- 05_amplicon_mapping.sbatch
|   |-- 06_novel_clustering.sbatch
|   |-- 07_final_table.sbatch
|   `-- common.sh
|-- docs/
|-- data/
|   |-- input/
|   |-- raw/
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
conda activate atlaspipe
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
cd /storage/microbiome/home_members/pesai/AtlasPipe
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

El archivo local `config/env.sh` está ignorado por Git y debe contener la ruta del proyecto específica de la máquina, la tabla de amplicones de entrada, la etiqueta del entorno y el número de muestras por clúster de comunidad. El script auxiliar `submit_pipeline_slurm.sh` también carga `config/env.sh` antes de enviar trabajos.

Cada banco de pruebas debe usar un `ATLASPIPE_RUN_LABEL` corto. Los datos crudos descargados de MicrobeAtlas se comparten, pero las salidas intermedias y finales se escriben en directorios específicos de cada ejecución:

```text
data/raw/
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
| `ATLASPIPE_ROOT` | Directorio raíz del proyecto |
| `ATLASPIPE_RUN_LABEL` | Nombre corto para el banco de pruebas/ejecución actual |
| `ATLASPIPE_INPUT_DIR` | Directorio para archivos de entrada proporcionados por el usuario |
| `MICROBEATLAS_RAW_DATA_DIR` | Directorio para archivos descargados de MicrobeAtlas |
| `ATLASPIPE_INTERMEDIATE_DIR` | Directorio para archivos intermedios de la pipeline |
| `ATLASPIPE_RESULTS_DIR` | Directorio para archivos de resultados finales |
| `ATLASPIPE_LOG_DIR` | Directorio para stdout/stderr de SLURM de la ejecución actual |
| `ATLASPIPE_DATA_DIR` | Alias retrocompatible de `ATLASPIPE_INTERMEDIATE_DIR` |
| `MICROBEATLAS_BIOM_TABLE` | Tabla BIOM de abundancias de MicrobeAtlas usada para seleccionar muestras presentes y extraer el subconjunto del paso 03 |
| `MICROBEATLAS_REFERENCE_FASTA_DIR` | Directorio del FASTA de referencia generado en el paso 04 |
| `ATLASPIPE_OWN_AMPLICON_MAPPING_DIR` | Directorio de FASTA de amplicones, mapeos contra OTU97 y clustering de no-hit |
| `ATLASPIPE_OWN_AMPLICON_TABLE_RAW` | Tabla de entrada con amplicones del usuario |
| `MICROBEATLAS_ENV_LABEL` | Etiqueta usada en los nombres de salida |
| `MICROBEATLAS_ENV_IDS` | IDs numericos de ambientes a seleccionar en el paso 02; por ejemplo `30` selecciona `aquatic` y sus subambientes, y `30,45` selecciona ambos arboles completos |
| `MICROBEATLAS_EXCLUDE_ENV_IDS` | IDs numericos opcionales de ambientes a excluir despues de la seleccion; tambien excluyen sus subambientes |
| `MICROBEATLAS_SAMPLES_PER_CC` | Número de muestras seleccionadas por clúster de comunidad |

Por defecto, `MICROBEATLAS_BIOM_TABLE` apunta a:

```text
data/raw/samples-otus.97.metag.minfilter.rod2025companion.biom
```

Este BIOM `rod2025companion` se usa como valor por defecto porque conserva muestras para todos los community clusters `AMPLICON + robust` representados en `community_cluster_mapping.tsv`. Esto permite que el paso 02 seleccione muestras por community cluster sin perder clústeres que sí existen en el mapping.

Los pasos enviados con `submit_pipeline_slurm.sh` o con `sbatch slurm/<paso>.sbatch`
guardan stdout y stderr dentro de `ATLASPIPE_LOG_DIR`. Por defecto, cada
banco de pruebas usa una carpeta independiente bajo `log/<run_label>/`, con un
par de archivos por paso y job:

```text
log/<run_label>/02_sample_selection_<job_id>.out
log/<run_label>/02_sample_selection_<job_id>.err
```

## Tabla de amplicones de entrada

Coloca la tabla de amplicones del usuario en:

```text
data/input/final_table_filtered.tsv
```

La tabla de amplicones del usuario debe contener al menos tres columnas. La pipeline interpreta las dos primeras por posición, no por nombre:

| Posición | Descripción |
|---|---|
| Columna 1 | Identificador original del amplicón proporcionado por el usuario |
| Columna 2 | Secuencia nucleotídica del amplicón |
| Columnas restantes | Abundancias de muestras, salvo columnas de metadatos/taxonomía reconocidas |

Formato exacto esperado:

- La tabla debe tener una fila de cabecera.
- Cada fila debe representar un amplicón.
- Los nombres de las dos primeras columnas no importan: la primera se toma como ID original y la segunda como secuencia.
- Los identificadores originales pueden tener cualquier formato no vacío, siempre que no estén duplicados.
- En el paso 01, la pipeline renombra internamente cada identificador propio con el prefijo `ext_` para evitar colisiones con IDs de MicrobeAtlas. La correspondencia se guarda en `own_amplicon_id_map.tsv`.
- Las secuencias deben contener solo caracteres nucleotídicos válidos (`A`, `C`, `G`, `T`, `U`, `N` o `-`).
- Cada columna de muestra debe contener abundancias numéricas.
- No debe haber amplicones duplicados ni nombres de muestra duplicados.

Ejemplo mínimo:

```text
id	sequence	muestra_1	muestra_2	muestra_3
ASV001	ACGTACGTACGT	10	0	3
OTU_45	TGCAACGTACGT	0	5	1
```

Todas las columnas que no sean metadatos se interpretan como columnas de abundancia de muestras. Las columnas de metadatos/taxonomía como `bim1`, `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus` y `Species` se ignoran para la extracción de abundancias. La matriz debe estar orientada con amplicones en filas y muestras en columnas; las tablas con amplicones en columnas y muestras en filas deben transponerse antes de usarlas como entrada.

## Contrato de la pipeline

| Paso | Script/SBATCH | Entrada principal | Salida principal | Entorno |
|---|---|---|---|---|
| 00 | `scripts/00_downloads.sh`, `slurm/00_downloads.sbatch` | URLs de descarga de MicrobeAtlas incluidas en el script | `data/raw/` | herramientas de shell |
| 01 | `scripts/01_prepare_inputs.R`, `slurm/01_prepare_inputs.sbatch` | `data/input/final_table_filtered.tsv` | `$ATLASPIPE_INTERMEDIATE_DIR/01_own_amplicon_inputs/amplicon_abundances.tsv`; `$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/amplicons.fasta` | `Renv` |
| 02 | `scripts/02_write_biom_sample_ids.py`, `scripts/02_sample_selection.R`, `slurm/02_sample_selection.sbatch` | `data/raw/community_cluster_metadata.tsv`; `data/raw/community_cluster_mapping.tsv`; `data/raw/samples-otus.97.metag.minfilter.rod2025companion.biom` | `$ATLASPIPE_INTERMEDIATE_DIR/02_filters/biom_sample_ids.txt`; `$ATLASPIPE_INTERMEDIATE_DIR/02_filters_<env>_<selection>/` | `python3` con `h5py`, después `Renv` |
| 03 | `scripts/03_biom_subset.sh`, `slurm/03_biom_subset.sbatch` | `data/raw/samples-otus.97.metag.minfilter.rod2025companion.biom`; `$ATLASPIPE_INTERMEDIATE_DIR/02_filters_<env>_<selection>/community_cluster_selection_<selection>_forbiomsubset.txt` | `$ATLASPIPE_INTERMEDIATE_DIR/03_samples_subsets/` | `Biom` |
| 04 | `scripts/04_reference_fasta.R`, `slurm/04_reference_fasta.sbatch` | `$ATLASPIPE_INTERMEDIATE_DIR/03_samples_subsets/samples-otus.97_selection_<env>_<selection>.tsv`; `data/raw/otus.97.allinfo` | `$ATLASPIPE_INTERMEDIATE_DIR/04_reference_fasta/` | `Renv` |
| 05 | `scripts/05_amplicon_mapping.sh`, `slurm/05_amplicon_mapping.sbatch` | `$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/amplicons.fasta`; `$ATLASPIPE_INTERMEDIATE_DIR/04_reference_fasta/MicrobeAtlas_selection_<env>_<selection>_otus.97.fasta` | `$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/OTUs97_vsearch*.blast6` | `MicrobeAtlas` |
| 06 | `slurm/06_novel_clustering.sbatch`, `scripts/06_novel_clustering_parse.R` | `$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/amplicons.fasta`; `$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/OTUs97_vsearch*.blast6` | `$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/amplicons_no_hits*`; `$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/amplicon_cluster_map*.tsv` | `MicrobeAtlas`, después `Renv` |
| 07 | `scripts/07_final_table.R`, `slurm/07_final_table.sbatch` | Salidas de los pasos 01, 03, 05 y 06 | `$ATLASPIPE_RESULTS_DIR/07_final_table/contextualized_abundance_table.tsv`; `$ATLASPIPE_RESULTS_DIR/07_final_table/amplicon_abundance_table.tsv` | `Renv` |

Notas importantes:

- La seleccion de ambientes en el paso 02 se hace solo por numero. Primero lista el catalogo con `Rscript scripts/02_sample_selection.R --list-envs` y despues usa los IDs de la columna `idx`.
- Seleccionar un ID padre incluye automaticamente sus subambientes. Por ejemplo, si el ID `30` corresponde a `aquatic`, se incluyen `aquatic` y todos los ambientes que empiezan por `aquatic;`.
- La exclusion tambien se hace por numero y se aplica despues de la inclusion. Por ejemplo, `--env-id 30 --exclude-env-id 32` selecciona todo el arbol de `30` salvo el arbol de `32`.
- Para seleccionar varios ambientes completos, separa los IDs con comas o espacios: `--env-id 30,45` o `--env-id 30 45`. Tambien puedes excluir varios: `--exclude-env-id 32,33`.
- El paso 02 escribe `sample_metadata_<selection>.tsv` y `selected_environment_counts_<selection>.tsv` para auditar qué ambientes exactos entraron.
- El paso 07 conserva la tabla final en cuentas crudas y añade `sample_metadata.tsv` y `sample_depths.tsv` para trazabilidad. La normalización queda fuera de la pipeline y debe hacerla cada análisis según su objetivo.

## Ejecución paso a paso

Ejecuta los pasos en orden.

### 00. Descargar los datos de MicrobeAtlas

Recomendado: SLURM.

```bash
sbatch slurm/00_downloads.sbatch
```

Ejecución interactiva alternativa:

```bash
bash scripts/00_downloads.sh
```

Salidas principales:

```text
data/raw/
```

### 01. Preparar las entradas de amplicones del usuario

Recomendado: interactivo o SLURM.

```bash
Rscript scripts/01_prepare_inputs.R
```

O bien:

```bash
sbatch slurm/01_prepare_inputs.sbatch
```

Salidas principales:

```text
$ATLASPIPE_INTERMEDIATE_DIR/01_own_amplicon_inputs/amplicon_abundances.tsv
$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/amplicons.fasta
```

### 02. Filtrar muestras de MicrobeAtlas

Recomendado: interactivo o SLURM.

Listar entornos disponibles:

```bash
Rscript scripts/02_sample_selection.R --list-envs
```

Ejecutar el filtrado:

```bash
python3 scripts/02_write_biom_sample_ids.py \
  --biom "$MICROBEATLAS_BIOM_TABLE" \
  --output "$ATLASPIPE_INTERMEDIATE_DIR/02_filters/biom_sample_ids.txt"

Rscript scripts/02_sample_selection.R \
  --env-id "$MICROBEATLAS_ENV_IDS" \
  --exclude-env-id "${MICROBEATLAS_EXCLUDE_ENV_IDS:-}" \
  --samples-per-cc "$MICROBEATLAS_SAMPLES_PER_CC"
```

O bien con SLURM:

```bash
sbatch slurm/02_sample_selection.sbatch
```

Salidas principales:

```text
$ATLASPIPE_INTERMEDIATE_DIR/02_filters/biom_sample_ids.txt
$ATLASPIPE_INTERMEDIATE_DIR/02_filters_<env>_<selection>/
```

El envoltorio SLURM del paso 02 genera automáticamente `biom_sample_ids.txt` a partir de la tabla BIOM antes de seleccionar muestras. Esto evita muestrear registros que existen en los metadatos de MicrobeAtlas pero están ausentes en la tabla BIOM de abundancias usada por el paso 03.

### 03. Extraer el subconjunto BIOM

Recomendado: SLURM.

```bash
bash scripts/03_biom_subset.sh
```

O bien:

```bash
sbatch slurm/03_biom_subset.sbatch
```

Salidas principales:

```text
$ATLASPIPE_INTERMEDIATE_DIR/03_samples_subsets/
```

### 04. Generar el FASTA de referencia de OTUs de MicrobeAtlas

Recomendado: interactivo o SLURM si el tamaño de los datos es alto.

```bash
Rscript scripts/04_reference_fasta.R
```

O bien:

```bash
sbatch slurm/04_reference_fasta.sbatch
```

Salidas principales:

```text
$ATLASPIPE_INTERMEDIATE_DIR/04_reference_fasta/
```

### 05. Mapear amplicones contra OTUs de MicrobeAtlas con VSEARCH

Recomendado: SLURM.

```bash
bash scripts/05_amplicon_mapping.sh
```

O bien:

```bash
sbatch slurm/05_amplicon_mapping.sbatch
```

Salidas principales:

```text
$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/OTUs97_vsearch*.blast6
```

### 06. Agrupar amplicones no mapeados

Recomendado: SLURM.

```bash
sbatch slurm/06_novel_clustering.sbatch
```

Este paso de SLURM extrae los amplicones sin aciertos contra OTUs de MicrobeAtlas, agrupa esos amplicones no mapeados con VSEARCH y después ejecuta `scripts/06_novel_clustering_parse.R` para construir el mapa de características sin acierto y el resumen de tamaños de clúster.

Salidas principales:

```text
$ATLASPIPE_INTERMEDIATE_DIR/05_amplicon_mapping/amplicons_no_hits*
```

### 07. Construir la tabla final de abundancias contextualizadas

Recomendado: interactivo o SLURM.

```bash
Rscript scripts/07_final_table.R
```

O bien:

```bash
sbatch slurm/07_final_table.sbatch
```

Salidas principales:

```text
$ATLASPIPE_RESULTS_DIR/07_final_table/
```

La tabla contextualizada principal es:

```text
$ATLASPIPE_RESULTS_DIR/07_final_table/contextualized_abundance_table.tsv
```

La tabla de amplicones preservada mantiene los amplicones originales del usuario como filas y añade su asignación contextual:

```text
$ATLASPIPE_RESULTS_DIR/07_final_table/amplicon_abundance_table.tsv
$ATLASPIPE_RESULTS_DIR/07_final_table/amplicon_abundance_table_long.tsv
$ATLASPIPE_RESULTS_DIR/07_final_table/amplicon_context_map.tsv
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
