# Metabarcoding de eDNA para mamíferos terrestres

Este repositorio documenta el flujo de trabajo bioinformático para detectar mamíferos terrestres en la UMA de Santa Cruz Nuevo, Puebla, a través del análisis de metabarcoding de DNA ambiental (eDNA).

---

Autora: Brenda Duana Hernández.

Universidad Nacional Autónoma de México

Posgrado en Ciencias Biológicas

---

## Contenido del repositorio

El flujo de trabajo se divide en cuatro grandes pasos:

**1. Preparación del entorno de rCRUX y PCR *in silico* con base local `mito`**, R.
-  [`01-setup-database.md`](01-setup-database.md)
-  [`01-setup-database.R`](01-setup-database.R)
-  [`02-pcr-in-silico.md`](02-pcr-in-silico.md)
-  [`02-pcr-in-silico.R`](02-pcr-in-silico.R)

**2. Base de referencia con `blast_seeds()` contra la base local `core_nt`**, en servidor, R [Análisis por generar]

**3. Procesamiento de lecturas crudas de secuenciación**, bash y R [Flujo generado con datos de prueba].
-  [`03-data.md`](03-data.md)
-  [`04-quality-control.md`](04-quality-control.md)
-  [`05-cutadapat.md`](05-cutadapat.md)
-  [`06-dada2.md`](06-dada2.md)

**4. Análisis estadísticos**, R [Análisis por generar].

El paso 1 y 2 está basado en el flujo de trabajo de [rCRUX](https://github.com/CalCOFI/rCRUX) (Curd et al., 2024), que permite construir una base de secuencias de referencia del mismo fragmento que amplifican los primers a partir de secuencias provenientes de NCBI, para asignar taxonomía a las secuencias de eDNA. 

Dentro de este repositorio se creo un [`diccionario de comandos`](diccionario-de-comandos.md) que fueron necesarios para al ejecución del código.

> *Cuando un paso se aparta de la documentación oficial de una herramienta, se cita textualmente la fuente (README, código o artículo de referencia)*

---

## Estructura del repositorio

Cada `.md` explica un paso o un grupo de pasos y su script está en la carpeta `scripts/` con el mismo número. En cada uno se describe:

1. **Qué hace:** Para qué sirve el paso dentro del flujo.
2. **Código:** Cada comando o grupo de comandos va precedido de una explicación de por qué se usa.
3. **Resultado:** La salida que se obtuvo.
4. **Decisión:** Qué se decidió a partir del resultado y qué deja para el siguiente paso.

---
 
## Software y versiones
 
| Herramienta | Versión |
|---|---|
| R | 4.4.0 (ucrt) |
| Rtools | 44 |
| rCRUX | 0.1.0.000 (commit `4c0aa9046d3636262d3f75b20d2a17f85381872c`) |
| BLAST+ | 2.17.0+ |
 
---

### Uso de inteligencia artificial

La verificación del código y la corrección de script que no funciona tal como lo documentan los autores se apoyaron en Claude (Anthropic, 2026), un modelo de lenguaje. Todo el código se ejecutó y verificó por la autora.

---

## Referencias

* Anthropic. (2026). *Claude* (versión Opus 5.5) [Modelo de lenguaje de gran tamaño]. https://claude.ai
* Curd, E. E., Gal, L., Gallego, R., Silliman, K., Nielsen, S., & Gold, Z. (2024). rCRUX: A rapid and versatile tool for generating metabarcoding reference libraries in R. *Environmental DNA*, 6, e489. https://doi.org/10.1002/edn3.489

