# 🧬 Taller de Análisis Genético - Antropología Molecular

¡Bienvenidos al repositorio de la parte teórica de la clase de Antropología Molecular! 

En este espacio alojaremos todos los scripts, protocolos y bases de datos necesarios para llevar a cabo nuestros análisis genéticos. El objetivo de este repositorio es que aprendan a construir un *pipeline* bioinformático reproducible: desde la descarga de secuencias crudas hasta la visualización de filogenias y estructura poblacional lista para publicación.

---

## 📂 Estructura del Repositorio

Para mantener el orden y asegurar la reproducibilidad de nuestros análisis, este repositorio está dividido en los siguientes directorios:

```text
📦 AntropologiaMolecular_Analisis
 ┣ 📂 Data
 ┃ ┣ 📂 Raw          # Datos crudos descargados de bases públicas (FASTA, VCF)
 ┃ ┗ 📂 Clean        # Datos filtrados y alineamientos listos para R
 ┣ 📂 Scripts
 ┃ ┣ 📜 00_Introduccion_Filogeneticos.R # Introducción al análisis de secuencias (ape)
 ┃ ┣ 📜 01_Bioinfo_Retrieval.sh    # Comandos de terminal para descarga de datos
 ┃ ┣ 📜 02_Phylogeny_Networks.R    # Script para modelos, árboles (ape) y redes
 ┃ ┗ 📜 03_PopStructure.R          # Script para PCA, DAPC y AMOVA (adegenet, poppr)
 ┣ 📂 Results
 ┃ ┣ 📂 Trees        # Árboles filogenéticos exportados (PDF/PNG)
 ┃ ┗ 📂 Plots        # Gráficos de estructura poblacional y redes de haplotipos
 ┗ 📜 README.md      # Este documento
