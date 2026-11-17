### PIPELINE rCRUX ###

# Fuente: README - CalCOFI/rCRUX
# https://github.com/CalCOFI/rCRUX

# Este script hace la PCR in silico (Paso 6) con dos sets de primers:
# 12SV5 (Riaz et al., 2011) y MiMammal-U (Ushio et al., 2017),
# y después compara los dos sets con énfasis en mamíferos.

# Requiere los Pasos 1-5 terminados (01-setup-database.md y 01_setup.R).
# NO correr todo de una vez: el paso 6.3 necesita internet y conviene
# revisar su resultado.

# La explicación de cada bloque está en 02-pcr-in-silico.md


# Paso 6. PCR in silico ---####

# 6.1 Cargar rCRUX, rutas y carpetas de resultados ---####

library(rCRUX)

accession_taxa_sql_path <- "C:/rcrux/taxonomizr_db/accessionTaxa.sql"
blast_db_path <- "C:/rcrux/blast_db/mito"

# rCRUX crea carpetas con dir.create() sin recursive = TRUE,
# así que las carpetas que contienen cada corrida deben existir antes
dir.create("C:/rcrux/outputs_pcr_insilico")

dir.create("C:/rcrux/outputs_pcr_insilico/01_rcrux")        # 6.2 y 6.3
dir.create("C:/rcrux/outputs_pcr_insilico/02_mamiferos")    # pregunta 1
dir.create("C:/rcrux/outputs_pcr_insilico/03_graficas")     # 6.5
dir.create("C:/rcrux/outputs_pcr_insilico/04_venn")         # pregunta 3
dir.create("C:/rcrux/outputs_pcr_insilico/05_cruce_sitio")  # pregunta 4


# 6.2 Vía local: get_seeds_local() con la base mito ---####

# Valores por defecto de rCRUX, excepto max_to_blast = 1
# (README: un primer a la vez para computadoras con 16 GB de RAM)

# 12SV5 (Riaz et al., 2011)
get_seeds_local(forward_primer_seq = "TAGAACAGGCTCCTCTAG",
                reverse_primer_seq = "TTAGATACCCCACTATGC",
                metabarcode_name = "12SV5",
                output_directory_path = "C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_local",
                accession_taxa_sql_path = accession_taxa_sql_path,
                blast_db_path = blast_db_path,
                max_to_blast = 1)
# 7,341 accesiones

# MiMammal (Ushio et al., 2017)
get_seeds_local(forward_primer_seq = "GGGTTGGTAAATTTCGTGCCAGC",
                reverse_primer_seq = "CATAGTGGGGTATCTAATCCCAGTTTG",
                metabarcode_name = "MiMammal",
                output_directory_path = "C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_local",
                accession_taxa_sql_path = accession_taxa_sql_path,
                blast_db_path = blast_db_path,
                max_to_blast = 1)
# 3,270 accesiones


# 6.3 Vía remota: get_seeds_remote() con Primer-BLAST de NCBI ---####

# Necesita internet.

# a) Correcciones 1-3 al código de rCRUX.
# Solo duran mientras está abierta la sesión de R: se corren en cada
# sesión nueva, después de library(rCRUX) y antes de get_seeds_remote()
# (ver 01-setup-database.md, sección "Vía remota con get_seeds_remote()").

# Corrección 1: iterative_primer_search(),
# inherits(x) == "try-error" -> inherits(x, "try-error")
txt <- deparse(rCRUX::iterative_primer_search)
txt <- gsub('inherits\\((\\w+)\\) == "try-error"',
            'inherits(\\1, "try-error")',
            txt)
fixed_iterative_primer_search <- eval(parse(text = txt))
environment(fixed_iterative_primer_search) <- asNamespace("rCRUX")
assignInNamespace("iterative_primer_search",
                  fixed_iterative_primer_search,
                  ns = "rCRUX")

# Corrección 2: primer_search(), num_primers fuera del if
fixed_primer_search <- function(forward, reverse, num_permutations = 25,
                                ..., .parallel = FALSE, .progress = 'none') {
  if (missing(forward) || missing(reverse)) {
    BLAST_primer()
  }
  primers <- enumerate_primers(forward,
                               reverse)

  num_primers <- nrow(primers)

  if (num_primers > num_permutations) {
    warning(immediate. = TRUE,
            'primers have ',
            num_primers,
            ' possible combinations due to degenerate bases, sampling ',
            num_permutations,
            " primers, use 'num_permutations' to change")
    primers <- primers[sample.int(num_primers,
                                  num_permutations,
                                  replace =
                                    FALSE), ]
  }
  message('BLASTing ', nrow(primers), ' primer combinations')
  plyr::alply(.data = primers,
              .margins = 1,
              .expand = FALSE,
              .parallel = .parallel,
              .progress = .progress,
              function(row) {
                BLAST_primer(row$forward, row$reverse, ...)
              })
}
environment(fixed_primer_search) <- asNamespace("rCRUX")
assignInNamespace("primer_search", fixed_primer_search, ns = "rCRUX")

# Corrección 3: filter_primer_hits(), usar los primers que recibe la función
txt <- deparse(rCRUX::filter_primer_hits)
txt <- gsub(".data$forward_primer_seq",
            "forward_primer_seq",
            txt,
            fixed = TRUE)
txt <- gsub(".data$reverse_primer_seq",
            "reverse_primer_seq",
            txt,
            fixed = TRUE)
fixed_filter_primer_hits <- eval(parse(text = txt))
environment(fixed_filter_primer_hits) <- asNamespace("rCRUX")
assignInNamespace("filter_primer_hits", fixed_filter_primer_hits, ns = "rCRUX")

# Parámetros de Primer-BLAST del ejemplo del README (get_seeds_remote):
# HITSIZE, evalue, word_size, MAX_TARGET_PER_TEMPLATE y NUM_TARGETS_WITH_PRIMERS.
# Correcciones 4 (argumentos con nombre) y 5 (crear la carpeta antes) incluidas.

# b) 12SV5 (Riaz et al., 2011)
dir.create("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_remote_params")
get_seeds_remote(forward_primer_seq = "TAGAACAGGCTCCTCTAG",
                 reverse_primer_seq = "TTAGATACCCCACTATGC",
                 output_directory_path = "C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_remote_params",
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
# Deben aparecer dos "BLAST alignment completed" y solo avisos file.remove(tmp);
# si aparece "Timeout was reached", repetir el comando

# c) MiMammal (Ushio et al., 2017)
# minimum_length y MAX_TARGET_SIZE: README, tabla "Pre-Made Version
# Controlled Reference Databases", fila MiFish Universal
dir.create("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_remote_params")
get_seeds_remote(forward_primer_seq = "GGGTTGGTAAATTTCGTGCCAGC",
                 reverse_primer_seq = "CATAGTGGGGTATCTAATCCCAGTTTG",
                 output_directory_path = "C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_remote_params",
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


# 6.4 Conteos de rCRUX por corrida ---####

# Archivos que rCRUX escribe al final de cada corrida con write.csv(tax_rank_sum, ...)
# (R/get_seeds_local.R y R/get_seeds_remote.R)

c12_local  <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_local/get_seeds_local/12SV5_filtered_get_seeds_local_unique_taxonomic_rank_counts.csv")
c12_remota <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_remote_params/get_seeds_local/12SV5_filtered_get_seeds_remote_unique_taxonomic_rank_counts.csv")
cmm_local  <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_local/get_seeds_local/MiMammal_filtered_get_seeds_local_unique_taxonomic_rank_counts.csv")
cmm_remota <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_remote_params/get_seeds_local/MiMammal_filtered_get_seeds_remote_unique_taxonomic_rank_counts.csv")

conteos <- rbind(c12_local, c12_remota, cmm_local, cmm_remota)
rownames(conteos) <- c("12SV5 local", "12SV5 remota", "MiMammal local", "MiMammal remota")
conteos
# species: 7,237 / 22,224 / 3,192 / 8,504 (rCRUX cuenta el NA como una especie)


# 6.5 Gráfica de los resultados ---####

# Tablas principales (taxonomía de cada accesión)
t12_local  <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_local/get_seeds_local/12SV5_filtered_get_seeds_local_output_with_taxonomy.csv")
t12_remota <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/12SV5_remote_params/get_seeds_local/12SV5_filtered_get_seeds_remote_output_with_taxonomy.csv")
tmm_local  <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_local/get_seeds_local/MiMammal_filtered_get_seeds_local_output_with_taxonomy.csv")
tmm_remota <- read.csv("C:/rcrux/outputs_pcr_insilico/01_rcrux/MiMammal_remote_params/get_seeds_local/MiMammal_filtered_get_seeds_remote_output_with_taxonomy.csv")

# Especies por clase: 12SV5
# Una fila por especie (con su clase); se quitan las filas sin nombre de especie
esp12_local  <- unique(subset(t12_local,  !is.na(species))[, c("class", "species")])
esp12_remota <- unique(subset(t12_remota, !is.na(species))[, c("class", "species")])

c12_local_clases  <- as.data.frame(table(esp12_local$class,  useNA = "ifany"))
c12_remota_clases <- as.data.frame(table(esp12_remota$class, useNA = "ifany"))
c12_local_clases$Via  <- "Local (mito)"
c12_remota_clases$Via <- "Remota (nt)"

clases12 <- rbind(c12_local_clases, c12_remota_clases)
colnames(clases12) <- c("Clase", "Especies", "Via")
clases12$Clase <- as.character(clases12$Clase)
clases12$Clase[is.na(clases12$Clase)] <- "Sin asignar"   # tortugas, cocodrilos, etc.
clases12

# Especies por clase: MiMammal
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

# Unir los dos sets; las clases con pocas especies van a "Otras clases"
clases12$Set <- "12SV5"
clasesmm$Set <- "MiMammal"
clases_todas <- rbind(clases12, clasesmm)

principales <- c("Actinopteri", "Mammalia", "Aves", "Amphibia",
                 "Lepidosauria", "Chondrichthyes", "Sin asignar")
clases_todas$Clase[!(clases_todas$Clase %in% principales)] <- "Otras clases"

# Total de especies de cada corrida (el número sobre cada barra)
totales <- aggregate(Especies ~ Via + Set, data = clases_todas, FUN = sum)
totales
# 7,236 / 22,223 / 3,191 / 8,503

# Gráfica de barras apiladas (cada barra suma 100 %)
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


# Comparación de los sets de primers ---####
# ANÁLISIS PROPIO: no es parte de rCRUX ni de su README.
# Usa las tablas t12_local, t12_remota, tmm_local y tmm_remota de 6.5.


# 1. Recuperación de especies de mamíferos ---####

# Solo mamíferos; las filas con class NA quedan fuera
# (?subset: "missing values are taken as false")
m12_local  <- subset(t12_local,  class == "Mammalia")
m12_remota <- subset(t12_remota, class == "Mammalia")
mmm_local  <- subset(tmm_local,  class == "Mammalia")
mmm_remota <- subset(tmm_remota, class == "Mammalia")

# Especies de mamíferos (nombres únicos, sin NA)
length(unique(na.omit(m12_local$species)))    # 1,452
length(unique(na.omit(m12_remota$species)))   # 3,202
length(unique(na.omit(mmm_local$species)))    # 1,465
length(unique(na.omit(mmm_remota$species)))   # 3,141

# Especies compartidas entre la vía local y la remota de cada set
# intersect() no quita los NA: primero se comprueba que no haya
anyNA(intersect(m12_local$species, m12_remota$species))   # FALSE
anyNA(intersect(mmm_local$species, mmm_remota$species))   # FALSE

length(intersect(m12_local$species, m12_remota$species))   # 1,435
length(intersect(mmm_local$species, mmm_remota$species))   # 1,459

# Listas de mamíferos por corrida (orden, familia, especie)
write.csv(unique(m12_local[, c("order", "family", "species")]),  "C:/rcrux/outputs_pcr_insilico/02_mamiferos/mamiferos_12SV5_local.csv",     row.names = FALSE)
write.csv(unique(m12_remota[, c("order", "family", "species")]), "C:/rcrux/outputs_pcr_insilico/02_mamiferos/mamiferos_12SV5_remota.csv",    row.names = FALSE)
write.csv(unique(mmm_local[, c("order", "family", "species")]),  "C:/rcrux/outputs_pcr_insilico/02_mamiferos/mamiferos_MiMammal_local.csv",  row.names = FALSE)
write.csv(unique(mmm_remota[, c("order", "family", "species")]), "C:/rcrux/outputs_pcr_insilico/02_mamiferos/mamiferos_MiMammal_remota.csv", row.names = FALSE)


# 2. Porcentaje de los mamíferos de la base que amplifica cada set ---####
# Solo vía local: la base nt no está en la computadora

# Accesiones de la base mito (BLAST+)
system('blastdbcmd -db C:/rcrux/blast_db/mito -entry all -outfmt "%a" -out C:/rcrux/outputs_pcr_insilico/mito_accesiones.txt')

acc_mito <- readLines("C:/rcrux/outputs_pcr_insilico/mito_accesiones.txt")
length(acc_mito)           # 18,042
length(unique(acc_mito))   # 18,042 (la base tiene 18,034 secuencias)

# Taxonomía de cada accesión, con la misma base que usa rCRUX
library(taxonomizr)

taxid_mito <- accessionToTaxa(acc_mito, accession_taxa_sql_path)
tax_mito <- as.data.frame(getTaxonomy(taxid_mito, accession_taxa_sql_path))
# avisos file.remove(tmp): inofensivos

# Mamíferos en mito y comprobación de que los recuperados están entre ellos
mito_mam <- subset(tax_mito, class == "Mammalia")
length(unique(na.omit(mito_mam$species)))                                 # 1,561

length(intersect(unique(na.omit(m12_local$species)), mito_mam$species))   # 1,452
length(intersect(unique(na.omit(mmm_local$species)), mito_mam$species))   # 1,465
# Cobertura: 12SV5 93.0 %, MiMammal 93.9 %


# 3. Especies compartidas entre sets y grupos donde difieren ---####

# eulerr: diagramas de Euler (área proporcional al número de especies)
# En Windows se instala la versión compilada (la de código fuente pide Rust)
if (!requireNamespace("eulerr", quietly = TRUE)) install.packages("eulerr", type = "binary")

library(eulerr)

# Listas de especies de mamíferos por set, una por vía
venn_esp_local  <- list("12SV5"    = unique(na.omit(m12_local$species)),
                        "MiMammal" = unique(na.omit(mmm_local$species)))
venn_esp_remota <- list("12SV5"    = unique(na.omit(m12_remota$species)),
                        "MiMammal" = unique(na.omit(mmm_remota$species)))

# Diagrama, vía local
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

# Diagrama, vía remota
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

# Número de especies de cada zona (columna original)
e_local    # 12SV5 79 / MiMammal 92 / 12SV5&MiMammal 1373
e_remota   # 12SV5 162 / MiMammal 101 / 12SV5&MiMammal 3040

# Tabla de especies: en qué set aparece cada una
# Vía local
esp_local <- sort(unique(c(venn_esp_local[["12SV5"]], venn_esp_local[["MiMammal"]])))
tab_esp_local <- data.frame(species  = esp_local,
                            "12SV5"  = esp_local %in% venn_esp_local[["12SV5"]],
                            MiMammal = esp_local %in% venn_esp_local[["MiMammal"]],
                            check.names = FALSE)
tab_esp_local$Via <- "Local (mito)"

# Vía remota
esp_remota <- sort(unique(c(venn_esp_remota[["12SV5"]], venn_esp_remota[["MiMammal"]])))
tab_esp_remota <- data.frame(species  = esp_remota,
                             "12SV5"  = esp_remota %in% venn_esp_remota[["12SV5"]],
                             MiMammal = esp_remota %in% venn_esp_remota[["MiMammal"]],
                             check.names = FALSE)
tab_esp_remota$Via <- "Remota (nt)"

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

# Revisar: deben ser los mismos números que e_local y e_remota
table(tab_esp$Via, tab_esp$grupo)

# Ordenar y guardar
tab_esp <- tab_esp[order(tab_esp$Via, tab_esp$grupo, tab_esp$order, tab_esp$family, tab_esp$species),
                   c("Via", "grupo", "order", "family", "species", "12SV5", "MiMammal")]

write.csv(tab_esp, "C:/rcrux/outputs_pcr_insilico/04_venn/venn_especies_mamiferos.csv", row.names = FALSE)

# Dónde se concentran las exclusivas
# Total de especies por orden en cada vía
tot_orden <- as.data.frame(table(Via = tab_esp$Via, order = tab_esp$order))
names(tot_orden)[3] <- "total"

# Exclusivas por orden, en cada vía y grupo
exclusivas <- subset(tab_esp, grupo != "Ambos sets")
exc_orden <- as.data.frame(table(Via = exclusivas$Via,
                                 grupo = exclusivas$grupo,
                                 order = exclusivas$order))
names(exc_orden)[4] <- "exclusivas"
exc_orden <- subset(exc_orden, exclusivas > 0)

# Porcentaje de exclusivas sobre el total del orden
exc_orden <- merge(exc_orden, tot_orden, by = c("Via", "order"))
exc_orden$porcentaje <- round(100 * exc_orden$exclusivas / exc_orden$total, 1)
exc_orden <- exc_orden[order(exc_orden$Via, exc_orden$grupo, -exc_orden$exclusivas), ]
exc_orden

# Exclusivas por familia, de la que tiene más a la que tiene menos
exc_familia <- as.data.frame(table(Via = exclusivas$Via,
                                   grupo = exclusivas$grupo,
                                   family = exclusivas$family))
names(exc_familia)[4] <- "exclusivas"
exc_familia <- subset(exc_familia, exclusivas > 0)
exc_familia <- exc_familia[order(exc_familia$Via, exc_familia$grupo, -exc_familia$exclusivas), ]

write.csv(exc_orden,   "C:/rcrux/outputs_pcr_insilico/04_venn/exclusivas_por_orden.csv",   row.names = FALSE)
write.csv(exc_familia, "C:/rcrux/outputs_pcr_insilico/04_venn/exclusivas_por_familia.csv", row.names = FALSE)


# 4. Mamíferos del sitio de estudio ---####

# "UTF-8-BOM": Excel guarda el .csv con una marca invisible al inicio
sitio <- read.csv("C:/rcrux/especies_sitio.csv", fileEncoding = "UTF-8-BOM")

# Tres primeras columnas con nombres cortos; sin espacios sobrantes
sitio <- sitio[, 1:3]
names(sitio) <- c("especie", "nombre_comun", "clase")
sitio$especie <- trimws(sitio$especie)

# Solo mamíferos; el género es la primera palabra del nombre
sitio <- subset(sitio, clase == "MAMMALIA")[, c("especie", "nombre_comun")]
sitio$genero <- sub(" .*", "", sitio$especie)
sitio
# 21 mamíferos

# Cruce por especie: grupo de cada especie del sitio en cada vía
grupo_local  <- subset(tab_esp, Via == "Local (mito)")[, c("species", "grupo")]
grupo_remota <- subset(tab_esp, Via == "Remota (nt)")[, c("species", "grupo")]
names(grupo_local)  <- c("especie", "local")
names(grupo_remota) <- c("especie", "remota")

# all.x = TRUE conserva las especies sin coincidencia; su NA pasa a "No recuperada"
cruce_esp <- merge(sitio, grupo_local,  by = "especie", all.x = TRUE)
cruce_esp <- merge(cruce_esp, grupo_remota, by = "especie", all.x = TRUE)
cruce_esp$local[is.na(cruce_esp$local)]   <- "No recuperada"
cruce_esp$remota[is.na(cruce_esp$remota)] <- "No recuperada"
cruce_esp

table(cruce_esp$local)    # Ambos sets 12 / Solo MiMammal 3 / No recuperada 6
table(cruce_esp$remota)   # Ambos sets 18 / No recuperada 3

# Cruce por género (para Dasypus y Rattus, identificados solo a género)
tab_esp$genero <- sub(" .*", "", tab_esp$species)

# TRUE si el set recupera al menos una especie del género
generos <- aggregate(cbind(`12SV5`, MiMammal) ~ Via + genero, data = tab_esp, FUN = any)

cruce_gen <- subset(generos, genero %in% sitio$genero)
cruce_gen$grupo <- ifelse(cruce_gen[["12SV5"]] & cruce_gen[["MiMammal"]], "Ambos sets",
                          ifelse(cruce_gen[["12SV5"]], "Solo 12SV5", "Solo MiMammal"))
cruce_gen <- cruce_gen[order(cruce_gen$Via, cruce_gen$grupo, cruce_gen$genero), ]
cruce_gen

# Géneros del sitio que ningún set recupera en ninguna vía
setdiff(unique(sitio$genero), cruce_gen$genero)   # character(0)

# Guardar
write.csv(cruce_esp, "C:/rcrux/outputs_pcr_insilico/05_cruce_sitio/cruce_sitio_especies.csv", row.names = FALSE)
write.csv(cruce_gen, "C:/rcrux/outputs_pcr_insilico/05_cruce_sitio/cruce_sitio_generos.csv",  row.names = FALSE)

# ¿Las especies que la vía local no recupera están en mito?
no_local <- subset(cruce_esp, local == "No recuperada")$especie
no_local[no_local %in% mito_mam$species]   # character(0): ninguna está en mito


# Conclusión ---####
#
# Con la PCR in silico los dos sets son equivalentes para el objetivo del estudio:
# recuperan un número parecido de mamíferos, amplifican el 93.0 % y el 93.9 % de
# los mamíferos de mito y recuperan las mismas 18 de las 21 especies del sitio en nt.
# Las diferencias están en grupos que no son el centro del estudio
# (murciélagos para 12SV5; Perissodactyla y Proboscidea para MiMammal).
#
# Como solo se puede secuenciar un set, todas las muestras se amplificarán con los
# dos y se secuenciará el que amplifique en más muestras (ver 02-pcr-in-silico.md,
# Conclusión).
