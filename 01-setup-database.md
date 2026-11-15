# Preparación del entorno de rCRUX

En estos cinco pasos se prepara todo lo que rCRUX necesita antes de correr una PCR *in silico*:

1. **Rtools:** la herramienta de compilación que Windows necesita para instalar paquetes de R.
2. **rCRUX y sus dependencias:** el paquete y los paquetes de R que requiere.
3. **BLAST+ (*Basic Local Alignment Search Tool*) :** los programas de NCBI (*National Center for Biotechnology Information*) con los que rCRUX busca los primers y extrae secuencias.
4. **`accessionTaxa.sql`:** la base de taxonomía que traduce los números de acceso de NCBI a su clasificación taxonómica.
5. **`mito`:** la base BLAST de secuencias mitocondriales en la que se buscan los primers.

Script asociado: `scripts/01-setup-database.R`.

---

## Paso 1. Verificar el entorno de R y Rtools

**Qué hace:** comprueba la versión de R y que Windows tenga una herramienta de compilación (`make`). Algunos paquetes de R se instalan desde código fuente y necesitan compilarse en la computadora. En Windows, esa herramienta viene en Rtools, y la versión debe corresponder a la de R (Rtools44 es la que va con R 4.4.x).

El README de rCRUX no menciona Rtools porque es un requisito de Windows y no del paquete.

### 1.1 Verificar la versión de R y buscar `make`

Primero se consulta qué versión de R está instalada con `R.version.string`. `Sys.which("make")` busca ese programa en el PATH (la lista de carpetas donde Windows busca programas) y si devuelve una cadena vacía, no hay herramienta de compilación.

```r
R.version.string

Sys.which("make")
```

*Resultado:* R `"R version 4.4.0 (2024-04-24 ucrt)"`. `Sys.which("make")` devolvió cadena vacía, es decir, Rtools no estaba instalado.

### 1.2 Instalar Rtools44 y volver a verificar

Se descargó Rtools44 (`rtools44-6459-6401.exe`, ~431 MB) desde [CRAN](https://cran.r-project.org/bin/windows/Rtools/rtools44/), se aceptaron los permisos de instalación y se reinició RStudio por completo. Después se repitió la verificación:

```r
Sys.which("make")
```

*Resultado:*"C:\\rtools44\\usr\\bin\\make.exe"

R 4.4.0 no tenía herramienta de compilación. Con Rtools44 instalado, `make` quedó detectado.

> *Nota, si `make` sigue sin aparecer*
>
> | Si pasa esto… | Haz lo siguiente |
> |---|---|
> | `Sys.which("make")` sigue devolviendo cadena vacía después de instalar Rtools | Cierra y abre RStudio por completo (no basta con reiniciar la sesión de R) y repite la verificación |
> | Persiste tras reiniciar | Confirma que instalaste Rtools44 y no otra versión. |

---

## Paso 2. Instalar rCRUX y sus dependencias

**Qué hace:** Instala rCRUX desde el repositorio de sus autores en GitHub (`CalCOFI/rCRUX`), porque no está en CRAN. `install_github()` no resolvió por sí sola cinco dependencias: tres de Bioconductor (`GenomeInfoDb`, `Biostrings`, `phyloseq`) y dos de GitHub (`phylotools`, `psadd`). Por eso se instalan primero, se verifica cada grupo y al final se instala rCRUX.

**Intentos que fallaron**

| # | Qué se intentó | Qué pasó | Qué se concluyó |
|---|---|---|---|
| 1 | `devtools::install_github("CalCOFI/rCRUX", build_vignettes = TRUE)` | Error de credenciales: Windows tenía guardado un PAT (Personal Access Token) de GitHub no válido | El repositorio es público y no necesita autenticación. Se borró la credencial con `gitcreds::gitcreds_delete()` (opción 2) |
| 2 | El mismo comando, ya sin la credencial | Tras ~50 min: `dependencies 'phylotools', 'psadd', 'phyloseq' are not available for package 'rCRUX'` | `phyloseq`, `Biostrings` y `GenomeInfoDb` son de Bioconductor. `phylotools` y `psadd` no estaban en CRAN en ese momento, aunque rCRUX los declara como dependencias. Hay que instalarlos aparte y antes. Declarar bien estas dependencias es un pendiente abierto en el repositorio de rCRUX ([issue #22](https://github.com/CalCOFI/rCRUX/issues/22) y [issue #23](https://github.com/CalCOFI/rCRUX/issues/23), 2023) |

### 2.1 Cargar devtools

`devtools` contiene `install_github()`, la función que instala rCRUX en 2.4. Se carga con `library()`.

```r
# install.packages("devtools")   # solo si no está instalado
library(devtools)
```

Los dos comandos siguientes solo son necesarios si `install_github()` falla por credenciales (intento 1). `gitcreds_delete()` muestra la credencial guardada y un menú; se eligió la opción 2 (*Delete these credentials*).

```r
install.packages("gitcreds")
gitcreds::gitcreds_delete()
```

### 2.2 Instalar las dependencias de Bioconductor

[Bioconductor](https://bioconductor.org/) es un repositorio distinto de CRAN y tiene su propio instalador, `BiocManager`. Esta línea instala `BiocManager` solo si falta:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

Los paquetes son grandes, así que el tiempo de espera de descarga se sube a 600 segundos. El cambio solo vale para la sesión actual. 

`update = FALSE` y `ask = FALSE` impiden que Bioconductor actualice otros paquetes o haga preguntas. Si de todos modos pregunta por actualizaciones, se responde `n`.

```r
options(timeout = 600)
BiocManager::install(c("GenomeInfoDb", "Biostrings", "phyloseq"),
                     update = FALSE, ask = FALSE)
```

Después, se comprueba la instalación. `%in%` busca cada nombre en la lista de paquetes instalados y devuelve `TRUE` si lo encuentra. Con esto se conrfima la instalación exitosa.

```r
c("GenomeInfoDb", "Biostrings", "phyloseq") %in% rownames(installed.packages())
```

*Resultado:* `TRUE TRUE TRUE`.

### 2.3 Instalar las dependencias de GitHub

`phylotools` y `psadd` se instalan desde el repositorio de cada autor con `remotes`, que se instala junto con `devtools`.

```r
remotes::install_github("helixcn/phylotools")
remotes::install_github("cpauvert/psadd")
```

Se verifican de la misma forma:

```r
c("phylotools", "psadd") %in% rownames(installed.packages())
```

*Resultado:* `TRUE TRUE`.

### 2.4 Instalar rCRUX

Con las cinco dependencias confirmadas se instala rCRUX con el comando del README oficial. Es el paso largo: tardó unos **54 minutos** y no debe interrumpirse.

```r
devtools::install_github("CalCOFI/rCRUX", build_vignettes = TRUE)
```

Luego se carga el paquete, se listan y cuentan sus objetos exportados, para confirmar que se puede usar.

```r
library(rCRUX)
ls("package:rCRUX")
length(ls("package:rCRUX"))
```

Salida de `ls("package:rCRUX")`:

```text
[1] "%>%"
"blast_datatable"
[3] "blast_seeds"
"combine_fasta_and_taxonomy"
[5] "compare_ref_db"
"convert_to_qiime2"
[7] "derep_and_clean_db"
"expand_multi_taxids"
[9] "filter_primer_hits"
"get_seeds_local"
[11] "get_seeds_remote"                        "get_taxonomy_from_accession"
[13] "has_blast_plus_binaries"                    "iterative_primer_search"
[15] "output_to_taxonomy_file"
"run_blastdbcmd"
[17] "run_blastdbcmd_blastn_and_aggregate_resuts"
"run_blastn"
[19] "run_primer_blastn"
"save_state"
```

`length()` devolvió **20 objetos**: 19 funciones de rCRUX más el operador `%>%`.

Por último, se registran la versión y el commit exactos. El repositorio sigue cambiando, y el commit es lo que permite reinstalar esta misma versión.

```r
packageVersion("rCRUX")
packageDescription("rCRUX")$RemoteSha
```

*Resultado:* versión `0.1.0.0` (en el archivo DESCRIPTION aparece como `0.1.0.000`), commit `4c0aa9046d3636262d3f75b20d2a17f85381872c`.

**rCRUX quedó instalado y funcional**, El siguiente requisito es BLAST+ (Paso 3).


> *Nota*
>
> | Si pasa esto… | Haz lo siguiente |
> |---|---|
> | Al cargar `devtools` aparece `package 'usethis' was built under R version 4.4.3` | Es un aviso inofensivo: el paquete se compiló con un R un poco más nuevo que el instalado (4.4.0) |
> | Aparece *"`install_github()` was deprecated in devtools 2.5.0"* | Es solo un aviso y la función sigue funcionando. Se usa `devtools` porque es lo que indica el README de rCRUX |
> | Aparece *"Skipping install… the SHA1 (4c0aa904) has not changed since last install"* | No es un error: rCRUX ya está instalado en esa versión. Agrega `force = TRUE` solo si quieres reinstalar (otros ~54 min) |
> | `devtools` da problemas | `install.packages("remotes")` y luego `remotes::install_github("CalCOFI/rCRUX")` *(no probado)* |

---

## Paso 3. Instalar BLAST+

**Qué hace.** Instala NCBI BLAST+ y comprueba que R lo encuentre. rCRUX ejecuta por debajo dos programas de BLAST+: `blastn`, que busca los primers en la base de datos, y `blastdbcmd`, que extrae secuencias. rCRUX no los instala, y sin ellos `get_seeds_local()` no puede correr.

### 3.1 Descargar BLAST+ y comprobar que el archivo no está dañado

BLAST+ se descarga desde el [FTP de NCBI](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/). El archivo depende del sistema operativo. Para Windows se descargaron `ncbi-blast-2.17.0+-win64.exe` y su archivo de verificación `ncbi-blast-2.17.0+-win64.exe.md5`.

Antes de instalar se comprueba que la descarga llegó completa. `readLines()` lee el código MD5 que NCBI calculó para el instalador y `tools::md5sum()` calcula el del archivo descargado. Si coinciden, el archivo está íntegro. Si no, hay que descargarlo otra vez.

```r
readLines("C:/Users/HP/Downloads/ncbi-blast-2.17.0+-win64.exe.md5")
tools::md5sum("C:/Users/HP/Downloads/ncbi-blast-2.17.0+-win64.exe")
```

Salida:

```text
[1] "b78629d1739d1a70bb0907a07f7ea207  ncbi-blast-2.17.0+-win64.exe"

C:/Users/HP/Downloads/ncbi-blast-2.17.0+-win64.exe 
                "b78629d1739d1a70bb0907a07f7ea207" 
```

*Resultado:* los dos códigos coinciden, así que el instalador está completo.

### 3.2 Instalar BLAST+

La instalación se hace fuera de R, ejecutando `ncbi-blast-2.17.0+-win64.exe`.

### 3.3 Verificar que R encuentra BLAST+

`has_blast_plus_binaries()` es una función de rCRUX que busca `blastn` y `blastdbcmd` y les pide su versión. Si los encuentra, devuelve `TRUE`.

```r
rCRUX::has_blast_plus_binaries()
```

Salida:

```text
[1] TRUE
attr(,"checked_binaries")
             blastn          blastdbcmd 
  "blastn: 2.17.0+" "BLASTD~1: 2.17.0+" 
```

*Resultado:* `BLASTD~1` es el nombre corto con el que Windows muestra `blastdbcmd`. No es un error.
BLAST+ 2.17.0+ quedó instalado y R encuentra `blastn` y `blastdbcmd`.

### Decisión
Se usó BLAST+ **2.17.0+**, una versión posterior al rango que verificaron los autores de rCRUX: *"Version 2.10.1+ through 2.13.0 are verified compatible with rCRUX"* ([README de rCRUX](https://github.com/CalCOFI/rCRUX)). 

Que funciona se confirmó aquí con `has_blast_plus_binaries()` y, más adelante, con la ejecución correcta de `get_seeds_local()`. El siguiente requisito es la base de taxonomía `accessionTaxa.sql` (Paso 4).


> *Nota*
>
> | Si pasa esto… | Haz lo siguiente |
> |---|---|
> | `has_blast_plus_binaries()` no devuelve `TRUE` | Cierra RStudio por completo, vuelve a abrirlo y repite la verificación. Si sigue igual, comprueba que el instalador corresponda a tu sistema operativo (`win64.exe` en Windows) |
> | Los dos códigos MD5 no coinciden | Borra el instalador, descárgalo de nuevo y repite la comprobación |

---

## Paso 4. Construir accessionTaxa.sql

**Qué hace:** Construye la base de datos local de taxonomía que rCRUX necesita para traducir números de accesión (por ejemplo `NC_012920.1`) a su clasificación taxonómica (especie, género, familia, etc.). 

Es una base SQLite creada con el paquete `taxonomizr` a partir de dos fuentes públicas de NCBI:

- el archivo de taxonomía `taxdump` (nombres y jerarquía de los taxones)
- dos tablas `accession2taxid`, `nucl_gb` y `nucl_wgs`, que relacionan cada acceso con su taxón.

El archivo termina en `.sql`, sin embargo, `accessionTaxa.sql` no es un script de texto con comandos SQL, es un archivo binario de base de datos SQLite, que se consulta desde R con las funciones de taxonomizr. No se descarga ya construido, se genera en la propia computadora a partir de:

- Taxonomía (`taxdump.tar.gz`): https://ftp.ncbi.nih.gov/pub/taxonomy/
- Tablas acceso → taxid (`accession2taxid`): https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/
- Paquete que la construye (taxonomizr): https://github.com/sherrillmix/taxonomizr

**Por qué no se usó la vía oficial.** La vía oficial es `prepareDatabase()`, que descarga y construye todo de una vez. Esa descarga automática, que usa FTP, falló cuatro veces seguidas por errores de conexión (por ejemplo, `Transferred a partial file [ftp.ncbi.nih.gov]: server did not report OK, got 426`). Por eso los archivos se descargaron a mano por HTTPS y la base se construyó corriendo por separado las funciones que `prepareDatabase()` usa por dentro. El resultado es el mismo tipo de base que produce el método oficial, y esta vía está prevista en el README de rCRUX.

> *Nota: ¿Que significa FTP y HTTPS?*
> *FTP y HTTPS.* Son dos protocolos, es decir, dos formas en que la computadora pide y recibe archivos de un servidor. Aparecen en la tabla de intentos y en los pasos 4.3 y 4.4.
> - **FTP** (*File Transfer Protocol*, protocolo de transferencia de archivos) es un protocolo antiguo hecho solo para transferir archivos. Usa conexiones separadas para los comandos y para los datos, por eso se corta con más facilidad en conexiones inestables. El código `426` del intento 1 es una respuesta propia de FTP: la transferencia se abortó. Es el protocolo que usa `prepareDatabase()`.
> - **HTTPS** (*Hypertext Transfer Protocol Secure*, protocolo seguro de transferencia de hipertexto) es el que usan los navegadores para abrir páginas web. Descarga por una sola conexión cifrada.
>
> El servidor de NCBI se llama `ftp.ncbi.nih.gov` por razones históricas, pero también responde por `https://`, que es como se descargaron los archivos.
 
 ---

### 4.1 Crear las carpetas de trabajo

Todo el proyecto vive en `C:/rcrux`, con una carpeta por tipo de contenido:

- `temp_downloads` guarda archivos temporales de SQLite durante la construccion (paso 4.5)
- `taxonomizr_db` guarda las descargas y la base de taxonomía

`showWarnings = FALSE` evita un aviso si la carpeta ya existe, y `recursive = TRUE` crea también `C:/rcrux` si todavía no existe.

```r
dir.create("C:/rcrux/temp_downloads", showWarnings = FALSE)

dir.create("C:/rcrux/taxonomizr_db", recursive = TRUE)
```

### 4.2 Cargar taxonomizr y definir dónde quedará la base

Primero se comprueba que `taxonomizr` esté instalado (debe devolver `TRUE`), se carga y se define la ruta final del archivo `.sql`. Esa ruta se guarda en un objeto porque se usa varias veces. Si se cierra R, la variable se pierde y hay que volver a definirla.

```r
"taxonomizr" %in% rownames(installed.packages())

library(taxonomizr)

accession_taxa_sql_path <- "C:/rcrux/taxonomizr_db/accessionTaxa.sql"
```

*Resultado:* `TRUE`. Al cargar `taxonomizr` aparece el aviso `package 'taxonomizr' was built under R version 4.4.3`: el paquete se compiló con una versión de R un poco más nueva que la instalada (4.4.0). No impide nada.

Para saber qué hace `prepareDatabase()` por dentro se muestra su código. En R, el nombre de una función sin paréntesis muestra su código sin ejecutarla:
 
```r
taxonomizr::prepareDatabase
```
 
El código llama, en este orden, a cinco funciones de taxonomizr:
 
| Función dentro de `prepareDatabase()` | Qué hace | Cómo se hizo |
|---|---|---|
| `getNamesAndNodes()` | Descarga y extrae el taxdump (`names.dmp`, `nodes.dmp`) | A mano desde el navegador + `untar()` (paso 4.3) |
| `read.names.sql()` | Carga los nombres en la base | Paso 4.5 |
| `read.nodes.sql()` | Carga la jerarquía en la base | Paso 4.5 |
| `getAccession2taxid()` | Descarga las tablas `accession2taxid` | Paso 4.4, con `protocol = "https"` |
| `read.accession2taxid()` | Carga la relación acceso → taxón en la base | Paso 4.5 |
 
Así se comprueba que hacerlo a mano equivale al método oficial. Aquí `read.nodes.sql()` se corrió antes que `read.names.sql()`; el orden no importa porque cada una llena su propia tabla.

### 4.3 Descargar la taxonomía de NCBI

`taxdump.tar.gz` contiene la jerarquía taxonómica (`nodes.dmp`) y los nombres (`names.dmp`). Como la descarga automática falló cuatro veces, se descargó desde el navegador por [HTTPS](https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz) y se guardó en `C:/rcrux/taxonomizr_db/`. `untar()` lo descomprime y `list.files()` muestra qué quedó en la carpeta.

```r
untar("C:/rcrux/taxonomizr_db/taxdump.tar.gz",
      exdir = "C:/rcrux/taxonomizr_db")
list.files("C:/rcrux/taxonomizr_db")
```
 
Salida (este comando se volvió a correr con el Paso 4 ya terminado, por eso también aparecen la base y las tablas de 4.4):
 
```text
 [1] "accessionTaxa.sql"           "citations.dmp"
 [3] "delnodes.dmp"                "division.dmp"
 [5] "gc.prt"                      "gencode.dmp"
 [7] "images.dmp"                  "merged.dmp"
 [9] "names.dmp"                   "nodes.dmp"
[11] "nucl_gb.accession2taxid.gz"  "nucl_wgs.accession2taxid.gz"
[13] "readme.txt"                  "taxdump.tar.gz"
```
 
Los dos archivos que se usan de aquí son `names.dmp` y `nodes.dmp`.

### 4.4 Descargar las tablas accession2taxid por HTTPS
 
`getAccession2taxid()` es la función que `prepareDatabase()` llama por dentro, y por defecto descarga por FTP. Aquí se ejecuta sola con `protocol = "https"`.

`types` pide las dos tablas de nucleótidos (los mismos valores por defecto de taxonomizr), y con `resume = TRUE`, si la descarga se corta basta con repetir el comando y sigue donde se quedó. 

Como los archivos pesan varios GB, el tiempo de espera se sube a 7200 s (120 min). Este cambio solo vale para la sesión actual.

```r
options(timeout = 7200)
taxonomizr::getAccession2taxid(
  outDir   = "C:/rcrux/taxonomizr_db",
  types    = c("nucl_gb", "nucl_wgs"),
  protocol = "https",
  resume   = TRUE
)
```
 
Salida:
 
```text
Download status: 1 done; 0 in progress. Total size: 2568.34 Mb (100%)... done!
Downloaded file: https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
Modified: 2026-08-01 03:21:28
Status: 200
Download status: 1 done; 0 in progress. Total size: 6423.63 Mb (100%)... done!
Downloaded file: https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
Modified: 2026-08-01 04:37:45
Status: 206
[1] "C:/rcrux/taxonomizr_db/nucl_gb.accession2taxid.gz"
[2] "C:/rcrux/taxonomizr_db/nucl_wgs.accession2taxid.gz"
```

### 4.5 Construir la base de datos

Con los cuatro archivos listos se ejecutan, en orden, las tres funciones que `prepareDatabase()` corre por dentro. Las tres escriben en el mismo archivo `accessionTaxa.sql`:
 
- `read.nodes.sql()` carga la jerarquía taxonómica (qué taxón pertenece a cuál).
- `read.names.sql()` carga los nombres de cada taxón.
- `read.accession2taxid()` carga la relación acceso → taxón. Es la más pesada.

`extraSqlCommand` hace que SQLite escriba sus archivos temporales en `temp_downloads`, una carpeta en un disco con espacio suficiente. La documentación de taxonomizr recomienda esta opción cuando la carpeta temporal por defecto tiene poco espacio.
 
```r
read.nodes.sql(
  nodeFile = "C:/rcrux/taxonomizr_db/nodes.dmp",
  sqlFile   = accession_taxa_sql_path
)
 
read.names.sql(
  nameFile = "C:/rcrux/taxonomizr_db/names.dmp",
  sqlFile   = accession_taxa_sql_path
)
 
read.accession2taxid(
  taxaFiles = c("C:/rcrux/taxonomizr_db/nucl_gb.accession2taxid.gz",
                "C:/rcrux/taxonomizr_db/nucl_wgs.accession2taxid.gz"),
  sqlFile   = accession_taxa_sql_path,
  extraSqlCommand = "PRAGMA temp_store_directory = 'C:/rcrux/temp_downloads'"
)
```
 
`read.accession2taxid()` muestra dos fases, *"Reading in values. This may take a while."* y *"Adding index. This may also take a while."*, y terminó sin interrupciones. En este equipo la construcción completa tardó unas **2 h 14 min**.

### 4.6 Verificar la base
 
Se comprueba que el archivo existe, cuánto pesa y que funciona. Para la prueba se busca un acceso conocido: `NC_012920.1`, el genoma mitocondrial humano de referencia, que debe devolver *Homo sapiens*.
 
```r
file.exists(accession_taxa_sql_path)
file.info(accession_taxa_sql_path)$size / 1e9
getTaxonomy(accessionToTaxa("NC_012920.1", accession_taxa_sql_path), accession_taxa_sql_path)
```
 
Salida:
 
```text
[1] TRUE
[1] 113.438
     domain      phylum     class     
9606 "Eukaryota" "Chordata" "Mammalia"
     order      family      genus 
9606 "Primates" "Hominidae" "Homo"
     species       
9606 "Homo sapiens"
Hubo 33 avisos (use warnings() para verlos)
```
 
La prueba recorrió las tres tablas: el acceso se traduce al taxid 9606 y el taxid a la clasificación completa.
 
*Resultado:* `accessionTaxa.sql` quedó construida (113.4 GB) con la taxonomía de NCBI del 2026-08-07 y las tablas `accession2taxid` del 2026-08-01, y traduce los accesos correctamente.
 
### Decisión

Como la base traduce los accesos correctamente, la ruta `C:/rcrux/taxonomizr_db/accessionTaxa.sql` es la que reciben las funciones de rCRUX en el argumento `accession_taxa_sql_path`.
 
Los archivos descargados (taxdump, `.dmp` y las dos tablas `.accession2taxid.gz`) se conservaron, para no tener que descargarlos de nuevo si hay que reconstruir la base. El README de rCRUX recomienda borrarlos cuando la base se construye a mano (10.07 GB), según los autores, conservarlos puede generar un aviso en rCRUX *(no observado hasta este paso)*.
 
> *Nota*
>
> | Si pasa esto… | Haz lo siguiente |
> |---|---|
> | `getTaxonomy()` muestra `Hubo 33 avisos`; con `warnings()` se ven como `In file.remove(tmp): no fue posible abrir el archivo 'C:\Users\HP\AppData\Local\Temp\Rtmp…\file…', motivo 'Permission denied'` | Inofensivo: `getTaxonomy()` crea archivos temporales para la consulta y Windows no deja borrarlos porque siguen abiertos. Aparecen en cada consulta y el resultado es correcto |


**Intentos que fallaron**
 
| # | Qué se intentó | Qué pasó | Qué se concluyó |
|---|---|---|---|
| 1 | `prepareDatabase(accession_taxa_sql_path)`, la función oficial que descarga y construye todo de una vez | Falló la descarga de `taxdump.tar.gz`: `Transferred a partial file [ftp.ncbi.nih.gov]: server did not report OK, got 426` | Error de conexión FTP. El avance quedó guardado en la carpeta temporal de R |
| 2 | El mismo comando con `resume = FALSE`, para descargar desde cero | La descarga llegó al 100 %, pero taxonomizr detectó que el archivo no coincidía con el de NCBI: `Local: 2c384486b007d967e9d125bc267f53f2` / `Remote: 13f279c3e479a0437e06e6b16c2541df` | El archivo se dañó en el camino. taxonomizr revisa el MD5 por su cuenta y no lo usó |
| 3 | El mismo comando en una sesión nueva de R | Se cortó casi al empezar (4 %): `download interrupted ... Download failed` | `options(timeout)` se reinicia en cada sesión y no se había vuelto a poner |
| 4 | `Sys.setenv(TMPDIR = "C:/rcrux/temp_downloads")`, borrar los restos con "taxdump" en esa carpeta y repetir con `resume = FALSE` | `transfer closed with 6819615 bytes remaining to read` | Cuarto fallo seguido por FTP: se dejó de usar la descarga automática|

---

## Paso 5. Descargar la base BLAST `mito`
 
**Qué hace:** descarga y valida la base de datos de secuencias en la que rCRUX hace la PCR *in silico*.
 
**Qué pide rCRUX.** El README de rCRUX indica que se necesita una base de nucleótidos local en formato BLAST (*"rCRUX requires a local blast-formatted nucleotide database"*). La usan dos funciones del pipeline:
 
- `get_seeds_local()`: la PCR *in silico*. Busca en la base las secuencias que contienen los dos primers.
- `blast_seeds()`: el paso siguiente. Usa esas secuencias como semilla para buscar en la base secuencias parecidas a las que les falta uno o ambos primers.

Como base, el README propone la `nt` completa de NCBI (~242 GB) y su descarga con `update_blastdb.pl`, una herramienta que NCBI incluye con BLAST+ (*"NCBI provides a tool (perl script) for downloading databases as part of the blast+ package"*). Aquí se descargó la base preformada `mito` de NCBI (*NCBI Genomic Mitochondrial Reference Sequences*, 194 MB) en lugar de `nt`, que contiene genomas mitocondriales de referencia y se usa **solo para la PCR *in silico** (`get_seeds_local()`); `blast_seeds()` queda para una fase posterior.

Usar `mito` para la PCR *in silico* tiene precedente en el artículo de rCRUX. Para COI, los autores corrieron `get_seeds_local()` con tres bases distintas, una de ellas la *"NCBI mito blast database"*, y después corrieron `blast_seeds()` contra `nt` y `mito` (Curd et al., 2024, Data S1 y Tabla S1).

**Por qué no se usó `update_blastdb.pl`.** Es un script de Perl, y Windows no trae Perl. Solo automatiza la descarga de los archivos `.tar.gz` que NCBI publica en https://ftp.ncbi.nlm.nih.gov/blast/db/, así que se descargaron a mano desde el navegador; el MD5 de 5.1 confirma que es el mismo archivo.
 
### 5.1 Descargar la base y comprobar que no está dañada
 
Se crea la carpeta donde vivirá la base. Después se descargan desde el navegador `mito.tar.gz` y su archivo de verificación `mito.tar.gz.md5` (https://ftp.ncbi.nlm.nih.gov/blast/db/mito.tar.gz, publicado el 2026-08-04) y se guardan en esa carpeta. Igual que con BLAST+ (paso 3.1), se compara el MD5 que publica NCBI con el del archivo descargado.
 
```r
dir.create("C:/rcrux/blast_db")
 
hash_esperado <- readLines("C:/rcrux/blast_db/mito.tar.gz.md5")
print(hash_esperado)
 
library(tools)
hash_obtenido <- md5sum("C:/rcrux/blast_db/mito.tar.gz")
print(hash_obtenido)
```
 
Salida:
 
```text
[1] "bf58a9a1426cc41f2d2ad3b9f97ee0c1  mito.tar.gz"
     C:/rcrux/blast_db/mito.tar.gz 
"bf58a9a1426cc41f2d2ad3b9f97ee0c1" 
```
 
*Resultado:* los dos códigos coinciden, así que la descarga está completa.
 
### 5.2 Descomprimir y validar la base
 
`untar()` extrae los archivos de la base y `list.files()` muestra qué quedó en la carpeta.
 
```r
untar("C:/rcrux/blast_db/mito.tar.gz", exdir = "C:/rcrux/blast_db")
list.files("C:/rcrux/blast_db")
```
 
Salida:
 
```text
 [1] "mito.ndb"  "mito.nhr"  "mito.nin"  "mito.nnd"  "mito.nni"  "mito.nog"
 [7] "mito.nos"  "mito.not"  "mito.nsq"  "mito.ntf"  "mito.nto"
[12] "mito.tar.gz"  "mito.tar.gz.md5"
[14] "taxdb.btd"  "taxdb.bti"  "taxonomy4blast.sqlite3"
```
 
Los 11 archivos `mito.*` forman la base en formato BLASTDB versión 5. La lista de 8 archivos que menciona el README de rCRUX (con `.nhd` y `.nhi`) corresponde a una versión anterior del formato. `taxdb.*` y `taxonomy4blast.sqlite3` vienen incluidos en `mito.tar.gz`.
 
Después se comprueba que BLAST+ reconoce la base. `blastdbcmd -info` lee la base y muestra su descripción:
 
```r
system('blastdbcmd -db "C:/rcrux/blast_db/mito" -info')
```
 
Salida:
 
```text
Database: NCBI Genomic Mitochondrial Reference Sequences
	18,034 sequences; 594,169,998 total bases
 
Date: Aug 3, 2026  12:48 AM	Longest sequence: 2,290,019 bases
 
BLASTDB Version: 5
 
Volumes:
	C:\rcrux\blast_db\mito
[1] 0
```
 
Por último se extrae una secuencia, igual que la prueba que propone el README de rCRUX con `blastdbcmd -entry`. El acceso del README (`MN937193.1`) es de GenBank y no forma parte de `mito`, así que se usa el genoma mitocondrial humano de referencia, el mismo acceso que se usó para probar `accessionTaxa.sql` (paso 4.6):
 
```r
system('blastdbcmd -db "C:/rcrux/blast_db/mito" -entry NC_012920.1 -range 1-60')
```
 
Salida:
 
```text
>NC_012920.1:1-60 Homo sapiens mitochondrion, complete genome
GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTT
[1] 0
```
 
### Resultado

La base `mito` quedó descargada, íntegra y funcional: 18,034 secuencias mitocondriales de referencia, versión del 3 de agosto de 2026, de la que BLAST+ lee la descripción y extrae secuencias.
 
### Decisión

Como la base quedó íntegra y BLAST+ la reconoce, la ruta `C:/rcrux/blast_db/mito` es la que recibe `get_seeds_local()` en el argumento `blast_db_path`.

**Limitación.** `mito` solo incluye genomas mitocondriales completos de referencia (RefSeq), así que las especies que en NCBI solo tienen fragmentos de 12S, 16S o COI no están en esta base. En el artículo de rCRUX, `get_seeds_local()` con `mito` devolvió 24,214 resultados de COI de 14,090 especies, contra 227,899 resultados de 39,545 especies con las bases de EMBL (Curd et al., 2024, Tabla S2). Por eso los números de la PCR *in silico* son una estimación mínima de lo que cada primer puede amplificar, y una especie ausente puede deberse a la base y no al primer.

> *Nota*
>
> | Si pasa esto… | Haz lo siguiente |
> |---|---|
> | Después de cada `system()` aparece `[1] 0` | No es un error: es el código de salida del programa, y `0` significa que terminó bien |
> | Los dos códigos MD5 no coinciden | Borra `mito.tar.gz`, descárgalo de nuevo y repite la comprobación |

---
 
## Vía remota con `get_seeds_remote()`
 
**Qué hace:** rCRUX tiene dos funciones para la PCR *in silico*. `get_seeds_local()` busca en una base BLAST guardada en la computadora (Paso 5). `get_seeds_remote()` no necesita base local: manda los primers a Primer-BLAST, la herramienta web de NCBI, y lee la página de resultados (*"submits a web form directly to NCBI's primer blast tool"*, Curd et al., 2024). 

Los autores la proponen para probar primers *"without downloading the large (>340 GB) nt database"*, aunque advierten que está sujeta a los límites del servidor de NCBI y señalan a `get_seeds_local()` como *"the preferred implementation in the rCRUX environment"*.

El argumento `organism` indica los grupos en los que se busca: `1476529` son los vertebrados sin mandíbula (lampreas y mixinos) y `7776` los vertebrados con mandíbula (todos los demás, incluidos los mamíferos).

El README de rCRUX indica correrla así:
 
```r
forward_primer_seq = "TAGAACAGGCTCCTCTAG"
reverse_primer_seq = "TTAGATACCCCACTATGC"
output_directory_path <- "/my/directory/12S_V5F1_remote_111122"
metabarcode_name <- "12S_V5F1"
accession_taxa_sql_path <- "/my/directory/accessionTaxa.sql"
 
get_seeds_remote(forward_primer_seq,
                 reverse_primer_seq,
                 metabarcode_name,
                 output_directory_path,
                 accession_taxa_sql_path,
                 organism = c("1476529", "7776"),
                 return_table = FALSE)
```

El ejemplo del README no corre tal como viene: el código de rCRUX en el commit `4c0aa904` tiene tres errores, y el ejemplo necesita dos ajustes más. Las cinco correcciones se describen abajo y ninguna cambia el análisis (primers, bases, filtros o parámetros), solo hacen que el código haga lo que los autores describen.

Para identificar los errores, el ejemplo se corrió tal cual, cambiando solo las rutas, en una sesión limpia de R y con rCRUX en el commit `4c0aa904`, el más reciente del repositorio. Cada vez que apareció un error se corrigió solo ese error y se volvió a correr. El código está en `scripts/01_setup.R`, sección "Nota: get_seeds_remote()", dentro de un bloque `if (FALSE)` para que no se ejecute al correr el script completo, porque cada consulta a NCBI tarda de 4 a 20 minutos.

### Correcciones
 
**Corrección 1.** En `R/iterative_primer_search.r`, a `inherits()` le falta el argumento `what`:
 
> `if (inherits(response) == "try-error") {` — `R/iterative_primer_search.r`
 
Error que produce: `Error en inherits(response): el argumento "what" está ausente, sin valor por omisión`. Corrección: `inherits(response, "try-error")`, y lo mismo con `parsed`. El mismo error está reportado en el repositorio de rCRUX con el ejemplo de la documentación ([issue #37](https://github.com/CalCOFI/rCRUX/issues/37)).
 
**Corrección 2.** En la función `primer_search()`, el cálculo de `primers` y `num_primers` está dentro de un `if` que no se ejecuta en el uso normal (cuando sí se dan los dos primers):
 
> ```r
> if (missing(forward) || missing(reverse)) {
>   BLAST_primer()
>   primers = enumerate_primers(forward, reverse)
>   num_primers = nrow(primers)
> }
> ```
> — `R/modifiedPrimerTreeFunctions.R`
 
Error que produce: `Output is not a data.frame` y el aviso `objeto 'num_primers' no encontrado`. Corrección: sacar esas dos líneas del `if`. Ese archivo está señalado por los autores como pendiente de limpieza ([issue #17](https://github.com/CalCOFI/rCRUX/issues/17)).
 
**Corrección 3.** En `filter_primer_hits()`, `.data$forward_primer_seq` busca una columna con ese nombre en la tabla de resultados, y esa columna no existe:
 
> `nchar(.data$forward_primer_seq) - nchar(.data$reverse_primer_seq))` — `R/filter_primer_hits.r`
 
Error que produce: `Column 'forward_primer_seq' not found in '.data'`. Corrección: usar los primers que recibe la función, `nchar(forward_primer_seq)`. Los propios autores dejaron la duda anotada en el mismo archivo:
 
> *"That means we either need to pass the forward and reverse primers to this function or we need to pass their lengths here. Or maybe that is supposed to happen in another function?"* — `R/filter_primer_hits.r`
 
**Corrección 4.** El ejemplo del README pasa los argumentos sin nombre y en un orden distinto al de la función: pone `metabarcode_name` antes que `output_directory_path`, y la función los define al revés:
 
> ```r
> get_seeds_remote <- function(forward_primer_seq,
>                              reverse_primer_seq,
>                              output_directory_path,
>                              metabarcode_name,
> ```
> — `R/get_seeds_remote.R`
 
Con el orden del README, el nombre del marcador se usa como carpeta de salida (aparece `Output directory: 12S_V5F1/get_seeds_local`). Corrección: llamar a la función con los argumentos nombrados. No es el único error reportado en el README: el [issue #31](https://github.com/CalCOFI/rCRUX/issues/31) señala que el README nombra `blast_seeds_local()` y `blast_seeds_remote()` donde debería decir `get_seeds_local()` y `get_seeds_remote()`.
 
**Corrección 5.** `get_seeds_remote()` crea la subcarpeta de resultados sin `recursive = TRUE`:
 
> `dir.create(out, showWarnings = FALSE)` — `R/get_seeds_remote.R`
 
Si la carpeta de salida no existe, la subcarpeta no se crea y la función falla al guardar: `no se puede abrir la conexión … No such file or directory`. Corrección: crear antes la carpeta de salida.
 
### Código
 
Las correcciones 1 a 3 se aplican durante la sesión de R con `assignInNamespace()`, que reemplaza la función dentro del paquete sin modificar los archivos instalados. **Hay que volver a aplicarlas cada vez que se abre R, después de `library(rCRUX)`**.
 
Las correcciones 1 y 3 cambian el texto de la función original: `deparse()` convierte la función en texto, `gsub()` reemplaza el fragmento con error y `eval(parse())` la vuelve a convertir en función.
 
```r
library(rCRUX)
 
# Corrección 1: iterative_primer_search()
txt <- deparse(rCRUX::iterative_primer_search)
txt <- gsub('inherits\\((\\w+)\\) == "try-error"', 'inherits(\\1, "try-error")', txt)
fixed_iterative_primer_search <- eval(parse(text = txt))
environment(fixed_iterative_primer_search) <- asNamespace("rCRUX")
assignInNamespace("iterative_primer_search", fixed_iterative_primer_search, ns = "rCRUX")
```
 
La corrección 2 reescribe `primer_search()` completa, igual a la original pero con el cálculo de `primers` y `num_primers` fuera del `if`:
 
```r
# Corrección 2: primer_search(), con el cálculo de num_primers fuera del if
fixed_primer_search <- function(forward, reverse, num_permutations = 25,
                                ..., .parallel = FALSE, .progress = 'none') {
  if (missing(forward) || missing(reverse)) {
    BLAST_primer()
  }
  primers <- enumerate_primers(forward, reverse)
  num_primers <- nrow(primers)
  if (num_primers > num_permutations) {
    warning(immediate. = TRUE, 'primers have ', num_primers,
            ' possible combinations due to degenerate bases, sampling ',
            num_permutations, " primers, use 'num_permutations' to change")
    primers <- primers[sample.int(num_primers, num_permutations, replace = FALSE), ]
  }
  message('BLASTing ', nrow(primers), ' primer combinations')
  plyr::alply(.data = primers, .margins = 1, .expand = FALSE,
              .parallel = .parallel, .progress = .progress,
              function(row) {
                BLAST_primer(row$forward, row$reverse, ...)
              })
}
environment(fixed_primer_search) <- asNamespace("rCRUX")
assignInNamespace("primer_search", fixed_primer_search, ns = "rCRUX")
```
 
```r
# Corrección 3: filter_primer_hits()
txt <- deparse(rCRUX::filter_primer_hits)
txt <- gsub(".data$forward_primer_seq", "forward_primer_seq", txt, fixed = TRUE)
txt <- gsub(".data$reverse_primer_seq", "reverse_primer_seq", txt, fixed = TRUE)
fixed_filter_primer_hits <- eval(parse(text = txt))
environment(fixed_filter_primer_hits) <- asNamespace("rCRUX")
assignInNamespace("filter_primer_hits", fixed_filter_primer_hits, ns = "rCRUX")
```
 
Con las tres funciones corregidas, se corre el ejemplo del README con las rutas de este proyecto, creando antes la carpeta de salida (corrección 5) y con los argumentos nombrados (corrección 4):
 
```r
forward_primer_seq      <- "TAGAACAGGCTCCTCTAG"
reverse_primer_seq      <- "TTAGATACCCCACTATGC"
output_directory_path   <- "C:/rcrux/outputs_remote/12S_V5F1_remote"
metabarcode_name        <- "12S_V5F1"
accession_taxa_sql_path <- "C:/rcrux/taxonomizr_db/accessionTaxa.sql"
 
dir.create(output_directory_path, recursive = TRUE, showWarnings = FALSE)
 
get_seeds_remote(forward_primer_seq      = forward_primer_seq,
                 reverse_primer_seq      = reverse_primer_seq,
                 output_directory_path   = output_directory_path,
                 metabarcode_name        = metabarcode_name,
                 accession_taxa_sql_path = accession_taxa_sql_path,
                 organism = c("1476529", "7776"),
                 return_table = FALSE)
 
warnings()
```
 
`warnings()` muestra los avisos de la corrida; sirve para comprobar que no se perdió ningún grupo de `organism`.
 
### Resultado
 
La función guardó sus tres archivos de salida, pero con el aviso `parse_primer_hits returned an object that is not a dataframe` y `This occurred while processing organism 7776`. El resultado tuvo 102 hits y 30 especies, 101 de ellos de la clase Hyperoartia (lampreas): funcionó para `1476529`, pero **perdió todo el grupo `7776`** sin detenerse.
 
Para saber si el problema estaba en la respuesta de NCBI, se consultó directamente a NCBI solo para `7776` con `primer_search()`, la función interna con la que rCRUX hace la consulta, y se leyó la respuesta con `parse_primer_hits()`. Este código no forma parte del flujo de rCRUX; solo revisa la página que devuelve NCBI:
 
```r
raw_response_list <- rCRUX:::primer_search(forward  = forward_primer_seq,
                                           reverse  = reverse_primer_seq,
                                           organism = "7776",
                                           primer_specificity_database = "nt")
r1 <- raw_response_list[[1]]
html_text <- httr::content(r1, as = "text", encoding = "UTF-8")
nchar(html_text)
length(gregexpr("<pre", html_text)[[1]])
p <- rCRUX:::parse_primer_hits(r1)
nrow(p)
```
 
La página trajo los resultados (448,893 caracteres y 1,000 etiquetas `<pre>`) y se leyó sin error: 1,000 resultados, por ejemplo *Bos grunniens*, *Elacatinus prochilos* y *Anser anser*. La misma consulta que antes no se pudo leer, esta vez sí: **la respuesta de NCBI varía entre consultas**. Además, que sean exactamente 1,000 indica un tope de resultados, que el README de rCRUX advierte:
 
> *"When using default parameters only 1047 hits are returned from NCBI's primer blast (run 11-11-22). Returns hit sizes and contents are variable depending on parameters, random blast sampling, and database updates."* — [README de rCRUX](https://github.com/CalCOFI/rCRUX)
 
Los tiempos que imprime rCRUX (*"BLAST alignment completed in 13.97 seconds"*) están en minutos aunque digan *seconds*: esa consulta tardó unos 14 minutos según `Sys.time()`.

### Ventajas y limitaciones frente a `get_seeds_local()`
 
| | `get_seeds_local()` + `mito` | `get_seeds_remote()` |
|---|---|---|
| Dónde busca | `mito`: 18,034 genomas mitocondriales de referencia (Paso 5) | `nt` completa, en los servidores de NCBI |
| Especies que en NCBI solo tienen fragmentos de 12S, 16S o COI | No las encuentra | Sí las encuentra (por ejemplo, accesos como `OZ482865.1` o `XR_014572079.1`, que no son genomas de referencia) |
| Funciona tal como viene | Sí | No: requiere las cinco correcciones |
| Resultados estables entre corridas | Sí | No: la misma consulta a veces se lee y a veces no |
| Tope de resultados | No | Unos 1,000 por consulta con los valores por defecto. El ejemplo avanzado del README lo sube con `HITSIZE = '1000000'` |
| Necesita base local | Sí | No |

### Decisión
 
`get_seeds_local()` con `mito` es la vía principal para comparar los primers, porque funciona como la publican los autores y da resultados estables y sin tope. `get_seeds_remote()` no es confiable por sí sola: necesita cinco correcciones, puede perder grupos enteros sin detenerse y su número de resultados tiene tope. Su ventaja es que busca en `nt`, así que encuentra especies que no están en `mito`.
 
Por eso cada set de primers se corre con las dos vías y se comparan los resultados: la remota sirve para notar diferencias y para revisar si una especie no detectada con `mito` falta por la base o por el primer. En cada corrida remota se aplican las correcciones y se revisan los avisos (`warnings()`) para confirmar que no se perdió ningún grupo.
 
---
 
## Referencias
 
Curd, E. E., Gal, L., Gallego, R., Silliman, K., Nielsen, S., & Gold, Z. (2024). rCRUX: A rapid and versatile tool for generating metabarcoding reference libraries in R. *Environmental DNA*, 6, e489. https://doi.org/10.1002/edn3.489

---
