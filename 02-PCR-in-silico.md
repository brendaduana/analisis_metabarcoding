# PCR in silico

Este documento corresponde al **Paso 6** del pipeline, es la continuación del [`01-setup-database.md`](01-setup-database.md), que documenta los pasos de preparación del equipo y de la base de datos.

Script asociado: [`02_pcr_in_silico.R`](02_pcr_in_silico.R)

**Qué hace:** Simula en la computadora una PCR con cada set de primers para saber cuántas especies y cuáles podría amplificar. No se usa laboratorio ni DNA real: el programa busca, en bases de datos de secuencias, dónde "se pegarían" los primers.

**Conceptos que se usan en este paso**

| Concepto| Significado |
|---|---|
| Primer | Secuencia corta de DNA (18–27 letras) que se pega al inicio (*forward*) o al final (*reverse*) de la región que se quiere copiar |
| Set de primers | El par forward + reverse |
| Amplicón | El fragmento de DNA que queda entre los dos primers, el que se copiaría en una PCR |
| PCR *in silico* | Buscar en secuencias de referencia dónde se pegan los dos primers, con la orientación y distancia correcta |
| Accesión | El identificador de una secuencia en NCBI, por ejemplo `NC_012920.1` (genoma mitocondrial humano) |
| *Seeds* (semillas) | Las secuencias que contienen los dos primers. rCRUX las usa después en `blast_seeds()` para construir la base de referencia |

**Qué pide rCRUX.** El [diagrama general del repositorio](https://github.com/CalCOFI/rCRUX/tree/main/flowcharts) muestra que la *PCR in silico* es el primer paso del flujo, y hay dos opciones:

1. `get_seeds_local()`: busca los primers con BLAST en una base de datos guardada en la computadora. Es más lenta, pero, según el diagrama, *"not subject to NCBI's memory or taxonomy limitations"*.
2. `get_seeds_remote()`: manda los primers a Primer-BLAST, la herramienta web de NCBI, que busca en su base `nt`.

Aquí se usan **las dos para cada set** y se comparan. 

La local busca en `mito` (genomas mitocondriales completos, Paso 5) y la remota en `nt` (secuencias públicas de NCBI).

**Sets de primers.** Son los primers adquiridos para el proyecto, con las secuencias tal como aparecen en los artículos originales:

| Set | Gen | Forward | Reverse | Tamaño reportado | Referencia |
|---|---|---|---|---|---|
| 12SV5 | 12S rRNA | `TAGAACAGGCTCCTCTAG` | `TTAGATACCCCACTATGC` | 73–110 bp | Riaz et al. (2011) |
| MiMammal-U | 12S rRNA | `GGGTTGGTAAATTTCGTGCCAGC` | `CATAGTGGGGTATCTAATCCCAGTTTG` | aprox 171 bp | Ushio et al. (2017) |

---

## Antes de empezar

Necesitas tener terminados los Pasos 1 a 5 del `01-setup-database.md`:

| Qué | Dónde debe estar |
|---|---|
| rCRUX instalado | dentro de R |
| BLAST+ instalado | en la computadora |
| Base de taxonomía | `C:/rcrux/taxonomizr_db/accessionTaxa.sql` |
| Base BLAST `mito` | `C:/rcrux/blast_db/mito.*` (11 archivos) |

---

## Paso 6. PCR *in silico*

### 6.1 Cargar rCRUX, definir las rutas y crear la carpeta de resultados
 
Se carga rCRUX, guarda las rutas de las dos bases de datos y se crea la carpeta de resultados del paso, con una subcarpeta por tipo de resultado:
 
```r
library(rCRUX)
 
accession_taxa_sql_path <- "C:/rcrux/taxonomizr_db/accessionTaxa.sql"
blast_db_path <- "C:/rcrux/blast_db/mito"
```

**Se deben crear las carpetasantes** rCRUX crea carpetas con `dir.create()` sin `recursive = TRUE`, así que solo puede crear un nivel a la vez:

- `get_seeds_local()` crea la carpeta de la corrida (`12SV5_local`). Por eso la carpeta que la contiene (`01_rcrux`) debe existir antes; si no, la corrida falla al guardar.

> `dir.create(output_directory_path)` — `R/get_seeds_local.R`

- `get_seeds_remote()` solo crea la subcarpeta `get_seeds_local` dentro de la carpeta de la corrida. Por eso, en la vía remota, también hay que crear antes la carpeta de la corrida (`12SV5_remote_params`; paso 6.3).

> `out <- file.path(output_directory_path, "get_seeds_local")`
> `dir.create(out, showWarnings = FALSE)` — `R/get_seeds_remote.R`

```r
dir.create("C:/rcrux/outputs_pcr_insilico")

# Subcarpetas para ordenar los resultados 
dir.create("C:/rcrux/outputs_pcr_insilico/01_rcrux") # 6.2 y 6.3
dir.create("C:/rcrux/outputs_pcr_insilico/02_mamiferos") # 6.5
dir.create("C:/rcrux/outputs_pcr_insilico/03_graficas") # 6.6
dir.create("C:/rcrux/outputs_pcr_insilico/04_venn") # 6.7
dir.create("C:/rcrux/outputs_pcr_insilico/05_cruce_sitio") #6.8
```

- `accession_taxa_sql_path`: base de taxonomía (Paso 4). Con ella rCRUX traduce cada accesión a su especie, género, familia, etc.
- `blast_db_path`: base `mito` (Paso 5). Se escribe sin extensión (`mito`, no `mito.nsq`), porque así la pide BLAST.

### 6.2 Vía local: `get_seeds_local()`

Para cada set, rCRUX:

1. Busca con BLAST el primer forward y el reverse en las 18,034 secuencias de `mito`.
2. Se queda con las secuencias donde aparecen los dos, en la orientación correcta para formar un amplicón.
3. Filtra por tamaño del amplicón y por número de diferencias (*mismatches*) entre el primer y la secuencia.
4. Busca la taxonomía de cada secuencia en `accessionTaxa.sql`.
5. Guarda los resultados en archivos `.csv`.

**Se corre primero 12SV5 y, cuando termina, MiMammal.**

```r
# 12SV5 (Riaz et al., 2011)
get_seeds_local(forward_primer_seq = "TAGAACAGGCTCCTCTAG",
                reverse_primer_seq = "TTAGATACCCCACTATGC",
                metabarcode_name = "12SV5",
                output_directory_path = "C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_local",
                accession_taxa_sql_path = accession_taxa_sql_path,
                blast_db_path = blast_db_path,
                max_to_blast = 1)
 
# MiMammal (Ushio et al., 2017)
get_seeds_local(forward_primer_seq = "GGGTTGGTAAATTTCGTGCCAGC",
                reverse_primer_seq = "CATAGTGGGGTATCTAATCCCAGTTTG",
                metabarcode_name = "MiMammal",
                output_directory_path = "C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_local",
                accession_taxa_sql_path = accession_taxa_sql_path,
                blast_db_path = blast_db_path,
                max_to_blast = 1)
```

¿Qué significa el código?

| Argumento | Requerimiento |
|---|---|
| `forward_primer_seq` | secuencia del forward |
| `reverse_primer_seq` | secuencia del reverse |
| `metabarcode_name` | nombre que rCRUX pone `"12SV5"` o `"MiMammal"` |
| `output_directory_path` | `".../12SV5_local"` rCRUX crea dentro una subcarpeta `get_seeds_local` |
| `accession_taxa_sql_path` | ruta del paso 6.1, Base de taxonomía |
| `blast_db_path` | ruta del paso 6.1, Base `mito` |
| `max_to_blast` | `1` Cuántos primers manda a BLAST a la vez |

*Todo lo que no aparece en el comando queda con el valor por defecto de rCRUX.*

Solo se cambia `max_to_blast`, de 2 a 1, porque el README lo recomienda para computadoras con 16 GB de RAM: *"blasting one at a time is more feasable for personal computers with 16 GB RAM"*.

El ejemplo del README también cambia `evalue` a 300, para que la búsqueda tarde menos con la base `nt` completa. Con `mito`, que es mucho más pequeña, no hace falta.

*Resultado:* esta es una salida real de 12SV5; la de MiMammal tiene la misma estructura.

```text
Output directory: C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_local/get_seeds_local
No previous primer blast files found. Starting pipeline from beginning.
Examining primers for degenerate bases.
  1 forward primer(s) will be blasted.
  1 reverse primer(s) will be blasted.
Reads will be blasted in subsets of up to 1 read(s). To change this, modify max_to_blast.
Calling blastn for primers. This may take a long time.
blastn -db C:/rcrux/blast_db/mito -task blastn-short -query ... -evalue 3e+07 -num_alignments 10000000 -qcov_hsp_perc 90 -perc_identity 50 -reward 2 -word_size 7 -num_threads 1
Calling blastn for primers. This may take a long time.
blastn -db C:/rcrux/blast_db/mito -task blastn-short -query ... (misma línea)
Blasting complete.
Wrangling results.
Joining with `by = join_by(saccver)`
Attaching taxonomies.
Done.
# A tibble: 7,341 × 37
```

- `1 forward primer(s) will be blasted`: los primers no tienen bases degeneradas (letras como R o Y, [ver más](https://www.bioinformatics.org/sms/iupac.html)), así que se busca una sola versión de cada uno.
- `Calling blastn` aparece **dos veces**: con `max_to_blast = 1`, rCRUX busca primero el forward y luego el reverse.
- La línea que empieza con `blastn` muestra exactamente cómo se hizo la búsqueda, con los valores por defecto de rCRUX.
- `# A tibble: 7,341 × 37`: la tabla final tiene **7,341 filas** (una por accesión; con MiMammal son **3,270**) y 37 columnas.

La corrida funcionó si aparece `Done.` y en `C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_local/get_seeds_local/` quedan cuatro archivos:

`12SV5_unfiltered_get_seeds_local_output.csv` - Pares de primers que forman un amplicón en la misma accesión, con la orientación correcta.

`12SV5_filtered_get_seeds_local_output_with_taxonomy.csv` - Lo anterior, filtrado por tamaño y mismatches, con la taxonomía de cada accesión. **Es el archivo principal y el que usará `blast_seeds()`.**

`12SV5_filtered_get_seeds_local_unique_taxonomic_rank_counts.csv` - Número de filos, clases, órdenes, familias, géneros y especies distintos (se usa en 6.4).

`12SV5_primers_selected_for_blastn.fasta` - Los primers que se buscaron.

Con **MiMammal**, los mismos archivos aparecen en `MiMammal_local/get_seeds_local/`, empezando con `MiMammal_`.

### 6.3 Vía remota: `get_seeds_remote()`

Manda los primers a Primer-BLAST, que busca en `nt`, una base mucho más grande que `mito`. Por eso **necesita internet** durante toda la corrida. rCRUX espera la respuesta de NCBI, la lee, le agrega la taxonomía con `accessionTaxa.sql` y guarda los resultados.

Tiene tres partes:

(a) aplicar las correcciones al código de rCRUX [Ver aqui `get_seeds_remote()`](01-setup-database.md#vía-remota-con-get_seeds_remote). Las correcciones solo duran mientras está abierta la sesión de R, así que se corren **en cada sesión nueva**, después de `library(rCRUX)` y antes de `get_seeds_remote()`,

(b) correr 12SV5,

```r
# 12SV5 (Riaz et al., 2011)
dir.create("C:/rcrux/outputs_pcr_insilico/12SV5_remote_params")
get_seeds_remote(forward_primer_seq = "TAGAACAGGCTCCTCTAG",
                 reverse_primer_seq = "TTAGATACCCCACTATGC",
                 output_directory_path = "C:/rcrux/outputs_pcr_insilico/12SV5_remote_params",
                 metabarcode_name = "12SV5",
                 accession_taxa_sql_path = accession_taxa_sql_path,
                 HITSIZE = '1000000',
                 evalue = '100000',
                 word_size = '6',
                 MAX_TARGET_PER_TEMPLATE = '5',
                 NUM_TARGETS_WITH_PRIMERS = '500000',
                 minimum_length = 50,
                 MAX_TARGET_SIZE = 200,
                 organism = c("1476529", "7776"),
                 return_table = FALSE)
warnings()
```

(c) correr MiMammal

```r
# MiMammal (Ushio et al., 2017)
dir.create("C:/rcrux/outputs_pcr_insilico/MiMammal_remote_params")
get_seeds_remote(forward_primer_seq = "GGGTTGGTAAATTTCGTGCCAGC",
                 reverse_primer_seq = "CATAGTGGGGTATCTAATCCCAGTTTG",
                 output_directory_path = "C:/rcrux/outputs_pcr_insilico/MiMammal_remote_params",
                 metabarcode_name = "MiMammal",
                 accession_taxa_sql_path = accession_taxa_sql_path,
                 HITSIZE = '1000000',
                 evalue = '100000',
                 word_size = '6',
                 MAX_TARGET_PER_TEMPLATE = '5',
                 NUM_TARGETS_WITH_PRIMERS = '500000',
                 minimum_length = 170,
                 MAX_TARGET_SIZE = 250,
                 organism = c("1476529", "7776"),
                 return_table = FALSE)
warnings()
```

**Qué significa el código.** Cada argumento de `get_seeds_remote()` y lo que recibe en este paso. Los valores son los del ejemplo del README, salvo los tamaños del amplicón de MiMammal, que se explican abajo:

| Argumento | Requerimiento |
|---|---|
| `forward_primer_seq`, `reverse_primer_seq`, `metabarcode_name`, `accession_taxa_sql_path` | igual que en el paso 6.2 |
| `output_directory_path` | `".../01_rcrux/<set>_remote_params"` rCRUX guarda los archivos en una subcarpeta `get_seeds_local`, aunque sea la vía remota paso 6.1 |
| `organism` | `c("1476529", "7776")` grupos de NCBI donde se busca: `1476529` vertebrados sin mandíbula (lampreas y mixinos) y `7776` vertebrados con mandíbula (peces, anfibios, reptiles, aves y mamíferos); rCRUX hace **una consulta por grupo** |
| `return_table` | `FALSE` no imprime la tabla en la consola (es muy grande), se guarda en archivo |
| `NUM_TARGETS_WITH_PRIMERS` | `'500000'` número máximo de secuencias que devuelve NCBI |
| `HITSIZE` | `'1000000'` número de resultados de BLAST que revisa NCBI antes de buscar los pares de primers |
| `evalue` | `'100000'` qué tan permisiva es la búsqueda; un valor más alto deja pasar más resultados |
| `word_size` | `'6'` coincidencia mínima, en letras, para empezar a alinear; más chico encuentra más |
| `MAX_TARGET_PER_TEMPLATE` | `'5'` número máximo de amplicones por secuencia (por defecto 100) |
| `MAX_TARGET_SIZE` | `200` (12SV5) / `250` (MiMammal) tamaño máximo del amplicón que reporta NCBI |
| `minimum_length` | `50` (12SV5) / `170` (MiMammal) filtro de rCRUX: tamaño mínimo del amplicón |

Con los valores por defecto, Primer-BLAST devuelve como máximo 1,000 secuencias por consulta, así que el número de especies quedaría limitado por NCBI y no por el primer.

La documentación de la función indica que `NUM_TARGETS_WITH_PRIMERS` se pone alto para obtener el mayor número de resultados y que su valor por defecto en NCBI es 1000:
 
> ```
> NUM_TARGETS_WITH_PRIMERS a primer BLAST search parameter set high to
>        maximize the number of observations returned.
>        The default is NCBI NUM_TARGETS_WITH_PRIMERS = 1000
> ```
> — `R/get_seeds_remote.R`, documentación de la función
 
El README lo advierte (*"When using default parameters only 1047 hits are returned"*) y da un ejemplo con parámetros modificados para 12SV5 que, según los autores, devuelve *"approximately 111500 blast seed returns"*. Aquí se usan esos mismos parámetros.
 
**Por qué MiMammal lleva otros tamaños.** El ejemplo del README es para 12SV5, cuyo amplicón mide unos 135 bp con los primers. El de MiMammal mide unos 220 bp con los primers (171 bp de inserto, según Ushio et al. 2017, más los dos primers), así que con `MAX_TARGET_SIZE = 200` NCBI no reportaría ningún amplicón.
 
Salida de 12SV5 (se omiten las líneas repetidas de *refreshing):

```text
Output directory: C:/rcrux/outputs_pcr_insilico/12SV5_remote_params/get_seeds_local
BLASTing 1 primer combinations
Submitting Primer-BLAST query
BLAST alignment processing, refreshing in 20 seconds...
BLAST alignment completed in 22.1111109256744 seconds
BLASTing 1 primer combinations
Submitting Primer-BLAST query
BLAST alignment processing, refreshing in 20 seconds...   (8 veces)
BLAST alignment completed in 3.42817306518555 seconds
Hubo 40 avisos (use warnings() para verlos)
```

Salida real de MiMammal:

```
Output directory: C:/rcrux/outputs_pcr_insilico/MiMammal_remote_params/get_seeds_local
BLASTing 1 primer combinations
Submitting Primer-BLAST query
BLAST alignment processing, refreshing in 20 seconds...   (2 veces)
BLAST alignment completed in 43.8542640209198 seconds
BLASTing 1 primer combinations
Submitting Primer-BLAST query
BLAST alignment processing, refreshing in 20 seconds...   (9 veces)
BLAST alignment completed in 3.41569996674856 seconds
Hubo 40 avisos (use warnings() para verlos)
```

Las corridas remotas se hicieron con la carpeta de la corrida directamente en `outputs_pcr_insilico`, como muestra la línea `Output directory`; después se movieron a `01_rcrux`, que es la ruta que usa el código.

Como se leeria la salida:

- `BLASTing 1 primer combinations` marca el inicio de una consulta a NCBI. Aparece **dos veces**: la primera para el grupo `1476529` y la segunda para `7776`.
- `refreshing in 20 seconds...` significa que NCBI sigue trabajando y rCRUX volverá a preguntar en 20 segundos. Mientras tanto, la consola no cambia: es normal.
- `BLAST alignment completed` indica que esa consulta terminó. **Deben aparecer dos.** El número que lo acompaña no siempre está en segundos: va sin unidades, y R lo expresa en segundos si la consulta duró menos de un minuto y en minutos si duró más, aunque el código de rCRUX (`R/modifiedPrimerTreeFunctions.R`) escribe siempre *seconds*. Por ejemplo, *"43.85 seconds"* (2 esperas de 20 s) son segundos y *"3.43 seconds"* (8 esperas) son minutos.
- `Hubo 40 avisos`: en `warnings()` deben aparecer solo avisos `In file.remove(tmp): … 'Permission denied'`, que son inofensivos: son archivos temporales de taxonomizr que Windows no deja borrar. El código de rCRUX (`R/get_seeds_local.R`) lo explica así: *"issue with taxonomizr remove.file that leads to warnings that can be omitted"*. Si aparece `Timeout was reached`, NCBI cortó la conexión y hay que repetir el comando.

*Resultado:* en las dos corridas aparecieron dos `BLAST alignment completed`, en `warnings()` no apareció `Timeout was reached`, y en `C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_remote_params/get_seeds_local/` quedaron tres archivos, con `_get_seeds_remote_` en el nombre.

Cada archivo contendrá informacion descriptiva de cada corrida:

- `12SV5_unfiltered_get_seeds_remote_output.csv`: todo lo que devolvió NCBI
- `12SV5_filtered_get_seeds_remote_output_with_taxonomy.csv`: **archivo principal** lo anterior, filtrado por tamaño y mismatches, con taxonomía.
- `12SV5_filtered_get_seeds_remote_unique_taxonomic_rank_counts.csv`: contiene número de filos, clases, órdenes, familias, géneros y especies distintos.

A diferencia de la vía local, la remota no deja una sola fila por accesión:

* 12SV5 dio 151,700 filas de 151,593 accesiones
* MiMammal 105,385 filas de 105,292 accesiones.

### 6.4 Conteos de rCRUX por corrida

Al terminar cada corrida, rCRUX esceibe un archivo con el numero de filos, clases, órdenes, familias, géneros y especies distintos. Aqui se leen los cuatro archivos y se juntan en una sola tabla para visualizar las 4 corridas juntas. Según el README, si en un nivel hay registros sin nombre (NA), *"they will be counted as a single unique rank"*: cuentan como uno más.

```r
c12_local  <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_local/get_seeds_local/12SV5_filtered_get_seeds_local_unique_taxonomic_rank_counts.csv")
c12_remota <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_remote_params/get_seeds_local/12SV5_filtered_get_seeds_remote_unique_taxonomic_rank_counts.csv")
cmm_local  <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_local/get_seeds_local/MiMammal_filtered_get_seeds_local_unique_taxonomic_rank_counts.csv")
cmm_remota <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_remote_params/get_seeds_local/MiMammal_filtered_get_seeds_remote_unique_taxonomic_rank_counts.csv")
 
conteos <- rbind(c12_local, c12_remota, cmm_local, cmm_remota)
rownames(conteos) <- c("12SV5 local", "12SV5 remota", "MiMammal local", "MiMammal remota")
conteos
```

*Resultados:* la columna `species` es el número de especies que amplifica cada set en cada vía. 12SV5 amplifica más especies que MiMammal, y la vía remota más que la local.

| | superkingdom | phylum | class | order | family | genus | species |
|---|---|---|---|---|---|---|---|
| 12SV5 local | 1 | 5 | 12 | 152 | 810 | 3,285 | 7,237 |
| 12SV5 remota | 1 | 2 | 9 | 158 | 986 | 6,304 | 22,224 |
| MiMammal local | 1 | 7 | 19 | 127 | 476 | 1,575 | 3,192 |
| MiMammal remota | 1 | 2 | 7 | 136 | 622 | 3,120 | 8,504 |

### 6.5 Gráfica de los resultados
 
Los conteos de 6.4 dicen cuantas especies amplifica cada set, pero no de que grupo son. Esta seccion resume las cuatro corridas en una grafica de barras aplicadas, que muestra cuantas especies recupera cada corrida y que proporción es de cada clase.

Primero se deben leer las tablas principales, se cargan los cuatro archivos `.csv`. A diferencia del archivo de conteos de 6.4, este tiene la taxonomia de cada accesión, asi que permite filtrar por clase y ver los nobres de las especies.

```r
t12_local  <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_local/get_seeds_local/12SV5_filtered_get_seeds_local_output_with_taxonomy.csv")
t12_remota <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_remote_params/get_seeds_local/12SV5_filtered_get_seeds_remote_output_with_taxonomy.csv")
tmm_local  <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_local/get_seeds_local/MiMammal_filtered_get_seeds_local_output_with_taxonomy.csv")
tmm_remota <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_remote_params/get_seeds_local/MiMammal_filtered_get_seeds_remote_output_with_taxonomy.csv")
```

*Resultado:* `read.csv()` no muestra nada en la consola. Quedan cargadas las cuatro tablas con su respectiva asignación,`t12_local`, `t12_remota`, `tmm_local`, `tmm_remota`.

Posteriormente, se cuenta cuántas especies de cada clase recupera cada corrida, primero con 12SV5:

```r
esp12_local  <- unique(subset(t12_local,  !is.na(species))[, c("class", "species")])
esp12_remota <- unique(subset(t12_remota, !is.na(species))[, c("class", "species")])
 
c12_local_clases  <- as.data.frame(table(esp12_local$class,  useNA = "ifany"))
c12_remota_clases <- as.data.frame(table(esp12_remota$class, useNA = "ifany"))
c12_local_clases$Via  <- "Local (mito)"
c12_remota_clases$Via <- "Remota (nt)"
 
clases12 <- rbind(c12_local_clases, c12_remota_clases)
colnames(clases12) <- c("Clase", "Especies", "Via")
clases12$Clase <- as.character(clases12$Clase)
clases12$Clase[is.na(clases12$Clase)] <- "Sin asignar"
clases12
```

y después con MiMammal:

```r
espmm_local  <- unique(subset(tmm_local,  !is.na(species))[, c("class", "species")])
espmm_remota <- unique(subset(tmm_remota, !is.na(species))[, c("class", "species")])
 
cmm_local_clases  <- as.data.frame(table(espmm_local$class,  useNA = "ifany"))
cmm_remota_clases <- as.data.frame(table(espmm_remota$class, useNA = "ifany"))
cmm_local_clases$Via  <- "Local (mito)"
cmm_remota_clases$Via <- "Remota (nt)"
 
clasesmm <- rbind(cmm_local_clases, cmm_remota_clases)
colnames(clasesmm) <- c("Clase", "Especies", "Via")
clasesmm$Clase <- as.character(clasesmm$Clase)
clasesmm$Clase[is.na(clasesmm$Clase)] <- "Sin asignar"
clasesmm
```

**Qué significa el código.** Los pasos tienen el mismo fundamento para ambos sets:
 
- `subset(..., !is.na(species))` quita las filas sin nombre de especie.
- `[, c("class", "species")]` se queda con las columnas de clase y especie, y `unique()` deja una sola fila por especie, aunque tenga varias accesiones.
- `table(..., useNA = "ifany")` cuenta las especies de cada clase, incluidas las que no tienen clase asignada (NA). `as.data.frame()` convierte ese conteo en una tabla.
- La columna `Via` anota de qué vía es cada conteo, y `rbind()` une las dos vías en una sola tabla.
- `colnames()` nombra las columnas y la última línea cambia el NA de la clase por `"Sin asignar"`.

*Resultado:* resumido por clase `clases12` y `clasesmm`.
 
| Clase | 12SV5 local | 12SV5 remota | MiMammal local | MiMammal remota |
|---|---|---|---|---|
| Actinopteri (peces óseos) | 3,845 | 10,114 | 425 | 2,058 |
| Mammalia (mamíferos) | 1,452 | 3,202 | 1,465 | 3,141 |
| Aves | 1,079 | 2,725 | 718 | 1,765 |
| Amphibia (anfibios) | 309 | 3,936 | 158 | 272 |
| Lepidosauria (lagartijas y serpientes) | 150 | 1,508 | 86 | 824 |
| Chondrichthyes (tiburones y rayas) | 219 | 446 | 210 | 361 |
| Sin asignar | 150 | 249 | 70 | 82 |
| Hyperoartia (lampreas) | 17 | 29 | — | — |
| Cladistia (bichires) | 11 | 14 | 1 | — |
| Clases fuera de los vertebrados | 4 (3 clases) | — | 58 (11 clases) | — |
| **Total** | **7,236** | **22,223** | **3,191** | **8,503** |
 
Las clases fuera de los vertebrados, solo en la vía local, son:
- 12SV5, Malacostraca (2), Agaricomycetes (1) y Rhynchonellata (1)
- MiMammal, Insecta (19), Gastropoda (16), Cephalopoda (8), Asteroidea (5), Polyplacophora (3), Palaeonemertea (2) y una especie de Branchiopoda, Enteropneusta, Malacostraca, Pilidiophora y Pycnogonida.

Después se unen las tablas de los dos sets para dibujarlos en la misma gráfica. Las clases con pocas especies no se lograrian diferenciar asi que se unen en el grupo "otra clase".

```r
clases12$Set <- "12SV5"
clasesmm$Set <- "MiMammal"
clases_todas <- rbind(clases12, clasesmm)
 
principales <- c("Actinopteri", "Mammalia", "Aves", "Amphibia",
                 "Lepidosauria", "Chondrichthyes", "Sin asignar")
clases_todas$Clase[!(clases_todas$Clase %in% principales)] <- "Otras clases"
 
totales <- aggregate(Especies ~ Via + Set, data = clases_todas, FUN = sum)
totales
```
 
**Significado del código:**
 
- `clases12$Set <- "12SV5"` agrega a cada tabla la columna `Set`, y `rbind()` las une.
- `principales` guarda los nombres de las clases que se dejan separadas. `%in%` revisa si cada clase está en esa lista y el `!` lo invierte, de modo que a las que **no** están se les cambia el nombre a `"Otras clases"`.
- `aggregate()` suma las especies de todas las clases de cada corrida, quedaria sobre cada barra.

*Resultado:* quedaron las tablas `clases_todas`, con ocho grupos de clase, y `totales`.
 
```text
           Via      Set Especies
1 Local (mito)    12SV5     7236
2  Remota (nt)    12SV5    22223
3 Local (mito) MiMammal     3191
4  Remota (nt) MiMammal     8503
```
 
Cada barra es una corrida y cada color una clase. Todas las barras estan con base en el 100%, asi que lo que se compara es la proporcion de especies de cada clase, el número sobre la barra es el total de especies.

```r
library(ggplot2)
 
g_apilada <- ggplot(clases_todas, aes(x = Via, y = Especies, fill = Clase)) +
  geom_col(position = "fill", width = 0.7) +
  geom_text(data = totales, aes(x = Via, y = 1.04, label = Especies),
            inherit.aes = FALSE, size = 3.5) +
  facet_wrap(~ Set) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Proporción de especies por clase",
       subtitle = "El número sobre cada barra es el total de especies",
       x = NULL, y = "Porcentaje de especies") +
  theme_minimal(base_size = 12)
g_apilada
 
ggsave("C:/rcrux/outputs_pcr_insilico/03_graficas/clases_apiladas.png",
       g_apilada, width = 8, height = 6, dpi = 300, bg = "white")
```

- `ggplot(clases_todas, aes(...))`: indica la tabla y qué columna va en cada parte de la gráfica: `Via` en el eje x, `Especies` en el eje y, y `Clase` en el color de relleno (`fill`).
- `geom_col()`: dibuja las barras. Con `position = "fill"` apila las clases y ajusta cada barra para que sume 100 %; `width = 0.7` es el ancho de la barra.
- `geom_text()`: escribe el total de especies sobre cada barra. Usa la tabla `totales` en lugar de `clases_todas`; `y = 1.04` lo coloca justo arriba del 100 %, e `inherit.aes = FALSE` hace que no use las columnas indicadas en `ggplot()`, porque `totales` no tiene la columna `Clase`.
- `facet_wrap(~ Set)`: divide la gráfica en dos paneles, uno por set.
- `scale_y_continuous(labels = scales::percent)`: escribe el eje y como porcentaje (25 %, 50 %…) en lugar de fracción (0.25, 0.50…).
- `scale_fill_brewer(palette = "Set2")`: usa la paleta de colores `Set2`, que tiene ocho colores, los mismos que grupos de clase.
- `labs()`: pone el título, el subtítulo y el título del eje y; `x = NULL` quita el título del eje x.
- `theme_minimal(base_size = 12)`: fondo sin recuadro y tamaño de letra 12.
`g_apilada` muestra la gráfica en RStudio y `ggsave()` la guarda en el archivo; `bg = "white"` le pone fondo blanco.

*Resultado:* `clases_apiladas.png` quedó en `C:/rcrux/outputs_pcr_insilico/03_graficas/`. 

| | Especies (total) | Especies de mamíferos |
|---|---|---|
| 12SV5 local (`mito`) | 7,236 | 1,452 |
| 12SV5 remota (`nt`) | 22,223 | 3,202 |
| MiMammal local (`mito`) | 3,191 | 1,465 |
| MiMammal remota (`nt`) | 8,503 | 3,141 |

- **Ningún set recupera solo mamíferos:** con los dos se esperan también secuencias de aves, peces, anfibios y reptiles.

---
 
## Comparación de los sets de primers
 
> [!IMPORTANT]
> **Análisis propio.** Desde aquí el código no es parte de rCRUX ni de su README: usa los archivos que rCRUX generó en el Paso 6 para comparar los dos sets de primers, con énfasis en mamíferos.

Con el siguiente código se comparan los dos sets con las siguientes preguntas como guía:

1. ¿Cuántas especies de mamíferos recupera cada set?
2. ¿Qué porcentaje de los mamíferos de la base amplifica cada set?
3. ¿Los dos sets recuperan las mismas especies, y en qué grupos difieren? 
4. ¿Recuperan a los mamíferos del sitio de estudio?

Usa las tablas `t12_local`, `t12_remota`, `tmm_local` y `tmm_remota`.

### 1. Recuperación de especies de mamíferos

rCRUX cuenta todas las especies que amplifica cada set (6.4). Esta sección separa los mamíferos de cada corrida, cuenta sus especies, compara la vía local con la remota y guarda las listas. Se comienza con la separación de mamíferos.

`subset()` se queda solo con las filas cuya columna `class` es `"Mammalia"`. Las filas sin clase asignada (NA) quedan fuera, porque `subset()` trata los valores faltantes como falsos:

```r
m12_local  <- subset(t12_local,  class == "Mammalia")
m12_remota <- subset(t12_remota, class == "Mammalia")
mmm_local  <- subset(tmm_local,  class == "Mammalia")
mmm_remota <- subset(tmm_remota, class == "Mammalia")
```

*Resultado:* quedan cargadas las cuatro tablas con su respectiva asignación,`m12_local`, `m12_remota`, `mmm_local`, `mmm_remota`.

Una vez asignados, se cuentan los nombre de especies disntintas de cada tabla. `na.omit()` quita los registros sin nombre de especie antes de contar, para que el NA no se cuente como una especie más, como si lo hace rCRUX en el paso 6.4.

```r
length(unique(na.omit(m12_local$species)))
length(unique(na.omit(m12_remota$species)))
length(unique(na.omit(mmm_local$species)))
length(unique(na.omit(mmm_remota$species)))
```

*Resultado:* en las dos vías, los dos sets recuperan un número parecido de especies de mamíferos.
 
| | Especies de mamíferos |
|---|---|
| 12SV5 local | 1,452 |
| 12SV5 remota | 3,202 |
| MiMammal local | 1,465 |
| MiMammal remota | 3,141 |

Para comparar las dos vías de cada set se usa `intersect()` que devuelve las especies que aparecen tanto en la vía local como en la remota del mismo set. A diferencia del paso anterior, aquí no se usa `na.omit()`, así que primero se comprueba con `anyNA()` que el resultado no incluye un NA, que se contaría como una especie compartida:

```r
anyNA(intersect(m12_local$species, m12_remota$species))
anyNA(intersect(mmm_local$species, mmm_remota$species))
```

Salida:

```text
[1] FALSE
[1] FALSE
```

*Resultado:* Como no hay NA, se procede a contar las especies compartidas:

```r
length(intersect(m12_local$species, m12_remota$species))
length(intersect(mmm_local$species, mmm_remota$species))
```

*Resultado:* casi todas las especies de mamíferos de la vía local están también en la remota: 1,435 de 1,452 con 12SV5 y 1,459 de 1,465 con MiMammal.

Se guarda una lista por corrida con el orden, la familia y el nombre de cada especie de mamifero. `unique()` deja una sola fila por combinación de orden, familia y especie, aunque la especie tenga varias accesiones.

```
write.csv(unique(m12_local[, c("order", "family", "species")]),  "C:/rcrux/outputs_pcr_insilico/02_mamiferos/mamiferos_12SV5_local.csv",     row.names = FALSE)
write.csv(unique(m12_remota[, c("order", "family", "species")]), "C:/rcrux/outputs_pcr_insilico/02_mamiferos/mamiferos_12SV5_remota.csv",    row.names = FALSE)
write.csv(unique(mmm_local[, c("order", "family", "species")]),  "C:/rcrux/outputs_pcr_insilico/02_mamiferos/mamiferos_MiMammal_local.csv",  row.names = FALSE)
write.csv(unique(mmm_remota[, c("order", "family", "species")]), "C:/rcrux/outputs_pcr_insilico/02_mamiferos/mamiferos_MiMammal_remota.csv", row.names = FALSE)
```

Como resultado se obtienen los cuatro archivos `.csv` con mamiferos en `C:/rcrux/outputs_pcr_insilico/02_mamiferos/`.

La pregunta 1 dice cuántas especies de mamíferos recupera cada set, pero no cuántas había para recuperar. Esta sección cuenta las especies de mamíferos que contiene la base `mito` y calcula qué porcentaje amplifica cada set. Solo se hace para la vía local: la base `nt` no está en la computadora, así que no se puede contar cuántas especies de mamíferos contiene.

Para ello se pide a BLAST+ la lista de accesiones de la base `mito` y se lee en R. Con `system()` se manda el comando al sistema desde R. `blastdbcmd` es el programa de BLAST+ que extrae información de una base: `-db` es la base, `-entry all` pide todas las secuencias, `-outfmt "%a"` pide solo la accesión de cada una y `-out` es el archivo donde se guarda la lista. mientras que con `readLines()` se lee el archivo, una accesión por línea. `length()` cuenta las accesiones y `length(unique())` las accesiones distintas.

```r
system('blastdbcmd -db C:/rcrux/blast_db/mito -entry all -outfmt "%a" -out C:/rcrux/outputs_pcr_insilico/mito_accesiones.txt')
 
acc_mito <- readLines("C:/rcrux/outputs_pcr_insilico/mito_accesiones.txt")
length(acc_mito)
length(unique(acc_mito))
```

*Resultado:* `mito_accesiones.txt` quedó en `C:/rcrux/outputs_pcr_insilico/`, con 18,042 accesiones distintas. Son 8 más que las 18,034 secuencias de la base (6.2), porque algunas secuencias tienen más de una accesión; no afecta el conteo, porque se cuentan nombres de especie distintos.

Posteriormente, se busca la taxonomia de cada accesion con `taxonomizr` con la misma base `accessionTaxa.sql` que usa rCRUX, para que los nombres de especie sean comparables con los de las corridas.
Recordando que `accessionToTaxa()` traduce cada accesión a su identificador taxonómico de NCBI, `getTaxonomy()` devuelve la taxonomía de cada identificador y `as.data.frame()` la convierte en tabla.

```r
library(taxonomizr)
 
taxid_mito <- accessionToTaxa(acc_mito, accession_taxa_sql_path)
tax_mito <- as.data.frame(getTaxonomy(taxid_mito, accession_taxa_sql_path))
```

Como *resultado* queda la tabla `tax_mito`, con una fila por accesión. Las dos funciones dejan avisos `In file.remove(tmp): … 'Permission denied'`, que son inofensivos.

Se separan los maiferos y se cuentan sus especies, después `intersect()` comprueba que todas las especies de mamíferos que recuperó cada set estén en esa lista; si no fuera así, el porcentaje no sería válido.

```r
mito_mam <- subset(tax_mito, class == "Mammalia")
length(unique(na.omit(mito_mam$species)))
 
length(intersect(unique(na.omit(m12_local$species)), mito_mam$species))
length(intersect(unique(na.omit(mmm_local$species)), mito_mam$species))
```
 
Salida:
 
```text
[1] 1561
[1] 1452
[1] 1465
```

*Resultado:* `mito` contiene 1,561 especies de mamíferos, y las 1,452 de 12SV5 y las 1,465 de MiMammal están todas entre ellas. La cobertura de cada set en la vía local es:
 
| Set | Mamíferos recuperados | Mamíferos en `mito` | Cobertura |
|---|---|---|---|
| 12SV5 | 1,452 | 1,561 | 93.0 % |
| MiMammal | 1,465 | 1,561 | 93.9 % |

Para responder si **los dos sets recuperan las mismas especies,** se utilizan las tablas `m12_local`, `m12_remota`, `mmm_local`, `mmm_remota` y el paquete [`eulerr`](https://cran.r-project.org/package=eulerr) (versión 7.1.0), que dibuja los círculos con un área proporcional al número de especies, de modo que el tamaño de cada zona refleja los datos.

Se arma una lista por vía, con las especies de mamíferos de cada set. Con `list()` guarda dos elementos con nombre, uno por set, que seran los que aparecen en el diagrama y `unique(na.omit())` deja los nombres de especie distintos, sin NA.


```r
if (!requireNamespace("eulerr", quietly = TRUE)) install.packages("eulerr", type = "binary")
 
library(eulerr)
 
venn_esp_local  <- list("12SV5"    = unique(na.omit(m12_local$species)),
                        "MiMammal" = unique(na.omit(mmm_local$species)))
venn_esp_remota <- list("12SV5"    = unique(na.omit(m12_remota$species)),
                        "MiMammal" = unique(na.omit(mmm_remota$species)))
```

*Resultado:* no muestra nada en la consola; quedaron las listas `venn_esp_local` y `venn_esp_remota`.

Con `euler()` se calcula, a partir de la lista, cuántas especies hay en cada zona, el tamaño y la posición de los círculos; `plot()` los dibuja. Se hace un diagrama por vía y se guarda cada uno con `png()` y `dev.off()`.

*Vía local*

```r
e_local <- euler(venn_esp_local)
d_local <- plot(e_local,
                quantities = TRUE,
                fills = list(fill = c("#E69F00", "#56B4E9"), alpha = 0.6),
                main = "Especies de mamíferos — vía local (mito)")
d_local
 
png("C:/rcrux/outputs_pcr_insilico/04_venn/venn_especies_local.png",
    width = 1800, height = 1500, res = 300)
print(d_local)
dev.off()
```

*Vía remota*

```r
e_remota <- euler(venn_esp_remota)
d_remota <- plot(e_remota,
                 quantities = TRUE,
                 fills = list(fill = c("#E69F00", "#56B4E9"), alpha = 0.6),
                 main = "Especies de mamíferos — vía remota (nt)")
d_remota
 
png("C:/rcrux/outputs_pcr_insilico/04_venn/venn_especies_remota.png",
    width = 1800, height = 1500, res = 300)
print(d_remota)
dev.off()
```

*Resultado:* `venn_especies_local.png` y `venn_especies_remota.png` quedaron en `C:/rcrux/outputs_pcr_insilico/04_venn/`. La zona compartida ocupa casi toda el área y las exclusivas quedan en un espcio reducido en las orillas, lo que indica que las especies exclusivas de cada set son pocas.

Al imprimir el objeto de `euler()` aparece una tabla con el número de especies de cada zona (columna `original`) y el área con la que se dibujó (`fitted`).
 
```r
e_local
e_remota
```

Salida:

```text
> e_local
               original fitted
12SV5                79     79
MiMammal             92     92
12SV5&MiMammal     1373   1373
 
> e_remota
               original fitted
12SV5               162    162
MiMammal            101    101
12SV5&MiMammal     3040   3040
```

Donde `12SV5` y `MiMammal` son las especies exclusivas de cada set y `12SV5&MiMammal` las compartidas. `original` es el número de especies y `fitted` el área con la que se dibujó la zona, como son iguales, el diagrama representa los datos sin distorsión.

Los diagramas dan números, no nombres. Por eso se arma una tabla con una fila por especie, que dice en qué set aparece, en qué vía, y su orden y familia.

*Vía local*

```r
esp_local <- sort(unique(c(venn_esp_local[["12SV5"]], venn_esp_local[["MiMammal"]])))
tab_esp_local <- data.frame(species  = esp_local,
                            "12SV5"  = esp_local %in% venn_esp_local[["12SV5"]],
                            MiMammal = esp_local %in% venn_esp_local[["MiMammal"]],
                            check.names = FALSE)
tab_esp_local$Via <- "Local (mito)"
```

*Vía remota*

```r
esp_remota <- sort(unique(c(venn_esp_remota[["12SV5"]], venn_esp_remota[["MiMammal"]])))
tab_esp_remota <- data.frame(species  = esp_remota,
                             "12SV5"  = esp_remota %in% venn_esp_remota[["12SV5"]],
                             MiMammal = esp_remota %in% venn_esp_remota[["MiMammal"]],
                             check.names = FALSE)
tab_esp_remota$Via <- "Remota (nt)"
```

Unen las dos vías y  `ifelse()` asigna el grupo: si la especie está en los dos sets, **Ambos sets**; si no, revisa si está en 12SV5 (**Solo 12SV5**) y, si tampoco, queda **Solo MiMammal**. Se escribe `tab_esp[["12SV5"]]` con doble corchete por el mismo nombre que empieza con número. `tax_esp` junta el orden y la familia de cada especie a partir de las cuatro tablas de mamíferos; `unique()` deja una fila por especie. `merge()` agrega esas dos columnas a `tab_esp`, usando el nombre de la especie; `all.x = TRUE` conserva todas las filas de `tab_esp`.

```r
# Unir las dos vías y agregar el grupo
tab_esp <- rbind(tab_esp_local, tab_esp_remota)
tab_esp$grupo <- ifelse(tab_esp[["12SV5"]] & tab_esp[["MiMammal"]],
                        "Ambos sets",
                        ifelse(tab_esp[["12SV5"]],
                               "Solo 12SV5",
                               "Solo MiMammal"))
 
# Agregar el orden y la familia de cada especie
tax_esp <- unique(rbind(m12_local[, c("order", "family", "species")],
                        m12_remota[, c("order", "family", "species")],
                        mmm_local[, c("order", "family", "species")],
                        mmm_remota[, c("order", "family", "species")]))
tab_esp <- merge(tab_esp, tax_esp, by = "species", all.x = TRUE)
 
# Revisar
table(tab_esp$Via, tab_esp$grupo)
```

Salida:
 
```text
               Ambos sets
  Local (mito)       1373
  Remota (nt)        3040
 
               Solo 12SV5
  Local (mito)         79
  Remota (nt)         162
 
               Solo MiMammal
  Local (mito)            92
  Remota (nt)            101
```

Se ordena la tabla y se guarda. `order()` ordena las filas por vía, grupo, orden, familia y especie, y el vector `c(...)` deja las columnas en orden de lectura. `write.csv()` guarda la tabla; `row.names = FALSE` evita que se guarde el número de fila.

```r
tab_esp <- tab_esp[order(tab_esp$Via, tab_esp$grupo, tab_esp$order, tab_esp$family, tab_esp$species),
                   c("Via", "grupo", "order", "family", "species", "12SV5", "MiMammal")]
 
write.csv(tab_esp, "C:/rcrux/outputs_pcr_insilico/04_venn/venn_especies_mamiferos.csv", row.names = FALSE)
```

Para saber donde se concentran las exclusivas, se cuentan las que hay en cada orden y se compara con el total de especies de ese orden en la misma via.

Total de especies por orden en cada vía:

```r
tot_orden <- as.data.frame(table(Via = tab_esp$Via, order = tab_esp$order))
names(tot_orden)[3] <- "total"
```

Exclusivas por orden, en cada vía y grupo:

```r
exclusivas <- subset(tab_esp, grupo != "Ambos sets")
exc_orden <- as.data.frame(table(Via = exclusivas$Via,
                                 grupo = exclusivas$grupo,
                                 order = exclusivas$order))
names(exc_orden)[4] <- "exclusivas"
exc_orden <- subset(exc_orden, exclusivas > 0)
```

Porcentaje de exclusivas sobre el total del orden:

```r
exc_orden <- merge(exc_orden, tot_orden, by = c("Via", "order"))
exc_orden$porcentaje <- round(100 * exc_orden$exclusivas / exc_orden$total, 1)
exc_orden <- exc_orden[order(exc_orden$Via, exc_orden$grupo, -exc_orden$exclusivas), ]
exc_orden
```

Exclusivas por familia, de la que tiene más a la que tiene menos:

```r
exc_familia <- as.data.frame(table(Via = exclusivas$Via,
                                   grupo = exclusivas$grupo,
                                   family = exclusivas$family))
names(exc_familia)[4] <- "exclusivas"
exc_familia <- subset(exc_familia, exclusivas > 0)
exc_familia <- exc_familia[order(exc_familia$Via, exc_familia$grupo, -exc_familia$exclusivas), ]
```

Se guardan las dos tablas, por orden y por familia.

```r
write.csv(exc_orden,   "C:/rcrux/outputs_pcr_insilico/04_venn/exclusivas_por_orden.csv",   row.names = FALSE)
write.csv(exc_familia, "C:/rcrux/outputs_pcr_insilico/04_venn/exclusivas_por_familia.csv", row.names = FALSE)
```

*Resultado:*`exclusivas_por_orden.csv` y `exclusivas_por_familia.csv` quedaron en `C:/rcrux/outputs_pcr_insilico/04_venn/`. Las exclusivas se concentran en pocos órdenes, así que no son diferencias al azar:

| Orden (familia principal) | Vía | Solo 12SV5 | Solo MiMammal |
|---|---|---|---|
| Chiroptera (Vespertilionidae, 55) | local | 55 de 175 | 3 |
| Chiroptera (Vespertilionidae 35, Phyllostomidae 15, Pteropodidae 8) | remota | 64 de 582 | 4 |
| Rodentia (Cricet0idae, 31) | remota | 51 de 964 | 11 |
| Primates (Cercopithecidae, 36) | local | 4 | 36 de 269 |
| Primates (Cercopithecidae, 36) | remota | 12 | 36 de 369 |
| Perissodactyla (Equidae 9, Rhinocerotidae 4, Tapiridae 3) | local | 0 | 16 de 16 (100 %) |
| Macroscelidea (Macroscelididae, 15) | local | 0 | 15 de 17 |
| Proboscidea (Elephantidae, 9) | remota | 0 | 12 de 12 |
| Pilosa | remota | 0 | 12 de 16 |

En conclusión, las de **12SV5** son sobre todo murciélagos (Chiroptera): en `mito`, MiMammal no amplifica el 31 % de las especies de murciélagos, todas de la familia Vespertilionidae. En `nt` se suman los roedores de la familia Cricetidae. Las de **MiMammal** son órdenes que 12SV5 amplifica poco o nada: en `mito`, 12SV5 no amplifica ninguna de las 16 especies de Perissodactyla (Equidae, Rhinocerotidae, Tapiridae), y en `nt` ninguna de las 12 de Proboscidea. Esto explica que, en la vía local, *Equus* solo lo recupere MiMammal.

Los diagramas de la pregunta 3 muestran que los dos sets comparten casi todos los mamiferos. Lo que podría inclinar la elección son las especies exclusivas. Para responder a la pregunta **4. Recuperan a los mamíferos del sitio de estudio** esta sección revisa, en cada vía y con que set se recupera, las opciones son: Ambos sets, Solo 12SV5, Solo MiMammal o No recuperada.

Para ello se necesita cargar el archivo `C:/rcrux/especies_sitio.csv`: la lista de especies del sitio (aves, reptiles y mamíferos), guardada desde Excel, con el nombre científico, el nombre común y la clase en las tres primeras columnas.

Comenzando con la carga del archivo:

```r
sitio <- read.csv("C:/rcrux/especies_sitio.csv", fileEncoding = "UTF-8-BOM")
```

Se dejan las tres primeras columnas y se quitan los espacios sobrantes de los nombres cientificos para que coincidan con los de rCRUX:

```r
sitio <- sitio[, 1:3]
names(sitio) <- c("especie", "nombre_comun", "clase")
sitio$especie <- trimws(sitio$especie)
```

Con`subset(..., clase == "MAMMALIA")` ayuda que solo queden los mamíferos y `sub(" .*", "", ...)` borra todo lo que va después del primer espacio, así que deja la primera palabra del nombre: el género.

```r
sitio <- subset(sitio, clase == "MAMMALIA")[, c("especie", "nombre_comun")]
sitio$genero <- sub(" .*", "", sitio$especie)
sitio
```

*Resultado:* la tabla `sitio` tiene los 21 mamíferos, con su nombre común y su género.

Después, se busca cada especie del sitio en `tab_esp` y se anota su grupo en la vía local y en la remota.
 
```r
grupo_local  <- subset(tab_esp, Via == "Local (mito)")[, c("species", "grupo")]
grupo_remota <- subset(tab_esp, Via == "Remota (nt)")[, c("species", "grupo")]
names(grupo_local)  <- c("especie", "local")
names(grupo_remota) <- c("especie", "remota")
```

Con `all.x = TRUE` conserva las especies del sitio que no están en `tab_esp`; quedan con NA, despues el NA se cambia por No recuperada:

```r
cruce_esp <- merge(sitio, grupo_local,  by = "especie", all.x = TRUE)
cruce_esp <- merge(cruce_esp, grupo_remota, by = "especie", all.x = TRUE)
cruce_esp$local[is.na(cruce_esp$local)]   <- "No recuperada"
cruce_esp$remota[is.na(cruce_esp$remota)] <- "No recuperada"
cruce_esp
```

Y se cuenta cuántas especies del sitio caen en cada grupo:

```r
table(cruce_esp$local)
table(cruce_esp$remota)
```

Salida (`cruce_esp`, ordenada por grupo):
 
| Especie | Nombre común | Local (`mito`) | Remota (`nt`) |
|---|---|---|---|
| *Bos taurus* | Vaca | Ambos sets | Ambos sets |
| *Canis latrans* | Coyote | Ambos sets | Ambos sets |
| *Canis lupus* | Perro | Ambos sets | Ambos sets |
| *Didelphis virginiana* | Tlacuache | Ambos sets | Ambos sets |
| *Herpailurus yagouaroundi* | Yaguarundí | Ambos sets | Ambos sets |
| *Homo sapiens* | Hombre | Ambos sets | Ambos sets |
| *Lynx rufus* | Lince | Ambos sets | Ambos sets |
| *Mephitis mephitis* | Zorrillo | Ambos sets | Ambos sets |
| *Odocoileus virginianus* | Venado cola blanca | Ambos sets | Ambos sets |
| *Procyon lotor* | Mapache | Ambos sets | Ambos sets |
| *Puma concolor* | Puma | Ambos sets | Ambos sets |
| *Urocyon cinereoargenteus* | Zorro | Ambos sets | Ambos sets |
| *Equus asinus* | Burro | **Solo MiMammal** | Ambos sets |
| *Equus caballus* | Caballo | **Solo MiMammal** | Ambos sets |
| *Oryctolagus cuniculus* | Conejo | **Solo MiMammal** | Ambos sets |
| *Bassariscus astutus* | Cacomixtle | No recuperada | Ambos sets |
| *Nasua narica* | Coatí | No recuperada | Ambos sets |
| *Peromyscus mekisturus* | Ratón | No recuperada | Ambos sets |
| *Mephitis macroura* | Zorrillo | No recuperada | No recuperada |
| *Dasypus unknown* | Armadillo | No recuperada | No recuperada |
| *Rattus unknown* | Rata | No recuperada | No recuperada |


*Resultado:* en la vía remota **los dos sets** recuperan las mismas 18 de las 21 especies del sitio. En la local, **MiMammal** recupera además el burro, el caballo y el conejo, los tres domésticos. *Ninguna especie del sitio es exclusiva de 12SV5.*

Dos registros del sitio están identificados solo a género (*Dasypus* y *Rattus*), así que no pueden coincidir por especie. Para ellos se revisa si cada set recupera al menos una especie del mismo género.

Se agrega a `tab_esp` el género de cada especie:

```r
tab_esp$genero <- sub(" .*", "", tab_esp$species)
```

`aggregate()` junta las especies por vía y género, y con `FUN = any` da `TRUE` en cada set si al menos una especie de ese género está en el set.

```r
generos <- aggregate(cbind(`12SV5`, MiMammal) ~ Via + genero, data = tab_esp, FUN = any)
```
 
Se dejan solo los géneros del sitio y se les asigna el grupo:
 
```r
cruce_gen <- subset(generos, genero %in% sitio$genero)
cruce_gen$grupo <- ifelse(cruce_gen[["12SV5"]] & cruce_gen[["MiMammal"]], "Ambos sets",
                          ifelse(cruce_gen[["12SV5"]], "Solo 12SV5", "Solo MiMammal"))
cruce_gen <- cruce_gen[order(cruce_gen$Via, cruce_gen$grupo, cruce_gen$genero), ]
cruce_gen
```
 
Salida `cruce_gen`:
 
| | Ambos sets | Solo MiMammal | Ausente |
|---|---|---|---|
| Vía local | 15 géneros | *Equus*, *Oryctolagus* | *Bassariscus* |
| Vía remota | los 18 géneros | — | — |
 
Por último, `setdiff()` devuelve los géneros del sitio que no aparecen en ninguna vía con ningún set:
 
```r
setdiff(unique(sitio$genero), cruce_gen$genero)
```
 
Salida:
 
```text
character(0)
```
 
*Resultado:* *Dasypus* y *Rattus* quedan en **Ambos sets** en las dos vías: los dos sets recuperan especies de esos géneros. En la vía local, *Equus* y *Oryctolagus* solo los recupera MiMammal, y *Bassariscus* no lo recupera ningún set. **Ningún género del sitio queda sin recuperar en las dos vías.**

Se guardan las tablas con los resultados obtenidos:

```r
write.csv(cruce_esp, "C:/rcrux/outputs_pcr_insilico/05_cruce_sitio/cruce_sitio_especies.csv", row.names = FALSE)
write.csv(cruce_gen, "C:/rcrux/outputs_pcr_insilico/05_cruce_sitio/cruce_sitio_generos.csv",  row.names = FALSE)
```

Por último, se revisa si las especies del sitio que la vía local no recupera están en `mito`, con la lista `mito_mam` de la pregunta 2:
 
```r
no_local <- subset(cruce_esp, local == "No recuperada")$especie
no_local[no_local %in% mito_mam$species]
```
 
Salida:
 
```text
character(0)
```

## Conclusión
 
Con la PCR *in silico*, los dos sets son equivalentes para el objetivo del estudio. Recuperan un número parecido de especies de mamíferos (1,452 y 1,465 en `mito`; 3,202 y 3,141 en `nt`), amplifican casi todos los mamíferos de `mito` (93.0 % y 93.9 %) y comparten la mayoría de las especies (1,373 en la vía local y 3,040 en la remota). Con los mamíferos del sitio pasa lo mismo: en `nt` recuperan las mismas 18 de las 21 especies, y la única ventaja de MiMammal en `mito` (el burro, el caballo y el conejo, los tres domésticos) desaparece en `nt`. Las especies del sitio que la vía local no recupera no están en `mito`, así que su ausencia se debe a la base y no a los primers.
 
Las diferencias entre los sets se concentran en grupos que no son el centro del estudio. 12SV5 amplifica murciélagos que MiMammal no (el 31 % de los de `mito`, todos Vespertilionidae) y, en `nt`, roedores cricétidos; MiMammal amplifica grupos donde 12SV5 falla casi por completo, como Perissodactyla en `mito` y Proboscidea en `nt`. Como la lista del sitio no incluye murciélagos y solo tiene dos roedores, la equivalencia para el sitio no cubre los grupos donde 12SV5 tiene ventaja.
 
Esta equivalencia es una predicción y puede sobreestimar cuánto coinciden los sets en campo. En aguajes de zonas semiáridas, Schenekar et al. (2026) usaron estos mismos dos sets: MiMammal detectó más especies de mamíferos (37 contra 29), pero 12SV5 detectó cinco que MiMammal no, y solo el 47 % de las especies las detectaron los dos. Por eso recomiendan usar *"at least two primer pairs to maximize species detection"*.
 
Como solo se puede secuenciar un set, la elección se hará en el laboratorio: todas las muestras se amplificarán con los dos sets y se secuenciará el que amplifique en más muestras, con la banda del tamaño esperado y sin bandas inespecíficas. Si los dos amplifican de forma parecida, se secuenciará MiMammal, porque una mayor proporción de lo que amplifica son mamíferos (6.5) y en campo detectó más especies de mamíferos (Schenekar et al., 2026). Si 12SV5 amplifica en más muestras, se secuenciará 12SV5, porque su amplicón es más corto y los fragmentos cortos son más abundantes en muestras ambientales (Schenekar et al., 2026).
 
**Referencias:**
- Curd, E. E., Gal, L., Gallego, R., Silliman, K., Nielsen, S., & Gold, Z. (2024). rCRUX: A rapid and versatile tool for generating metabarcoding reference libraries in R. *Environmental DNA*, 6, e489. https://doi.org/10.1002/edn3.489
- Riaz, T., Shehzad, W., Viari, A., Pompanon, F., Taberlet, P., & Coissac, E. (2011). ecoPrimers: inference of new DNA barcode markers from whole genome sequence analysis. *Nucleic Acids Research*, 39(21), e145. https://doi.org/10.1093/nar/gkr732
- Schenekar, T., Baxter, J., Sedlmayr, I., Gladitsch, J., Mahlangu, S., & Mwale, M. (2026). Unveiling vertebrate biodiversity in arid and semi-arid terrestrial ecosystems through eDNA metabarcoding at savanna waterholes. *Evolutionary Applications*, 19, e70200. https://doi.org/10.1111/eva.70200
- Ushio, M., Fukuda, H., Inoue, T., Makoto, K., Kishida, O., Sato, K., Murata, K., Nikaido, M., Sado, T., Sato, Y., Takeshita, M., Iwasaki, W., Yamanaka, H., Kondoh, M., & Miya, M. (2017). Environmental DNA enables detection of terrestrial mammals from forest pond water. *Molecular Ecology Resources*, 17, e63–e75. https://doi.org/10.1111/1755-0998.12690
