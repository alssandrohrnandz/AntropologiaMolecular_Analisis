# 🧬 Taller de Análisis Genético - Antropología Molecular

¡Bienvenidos al repositorio de la parte teórica de la clase de Antropología Molecular! 

En este espacio alojaremos todos los scripts, protocolos y bases de datos necesarios para llevar a cabo nuestros análisis genéticos. El objetivo de este repositorio es que aprendan a construir un *pipeline* bioinformático reproducible: desde la descarga de secuencias crudas hasta la visualización de filogenias y estructura poblacional lista para publicación.

---

## 📂 Estructura del Repositorio

Para mantener el orden y asegurar la reproducibilidad de nuestros análisis, este repositorio está dividido en los siguientes directorios:

```text
📦 AntropologiaMolecular_Analisis
 ┣ 📂 Examenes
 ┃ ┗ 📝 00_Examen_Practico.R  
 ┣ 📂 Data
 ┃ ┣ 📂 Raw          # Datos crudos descargados de bases públicas (FASTA, VCF)
 ┃ ┗ 📂 Clean        # Datos filtrados y alineamientos listos para R
 ┣ 📂 Scripts
 ┃ ┣ 📜 00_Introduccion_Filogeneticos.R # Introducción al análisis de secuencias (ape)
 ┃ ┣ 📜 01_Distancias_Geneticas.R    # Conversión y cálculo de la distancia Fst
 ┃ ┣ 📜 02_Estructura_Poblacional.R    # Script para AMOVA
 ┃ ┗ 📜 03_PopStructure.R          # Script para PCA, DAPC y AMOVA (adegenet, poppr)
 ┣ 📂 Results
 ┃ ┣ 📂 00_Introduccion        # Árboles filogenéticos exportados (PDF/PNG)
 ┃ ┣ 📂 01_Distancias_Geneticas        # Gráfico de heatmap de distancias genéticas
 ┃ ┗ 📂 02_Estructura_Poblacional #Grafico de Pastel donde se observa dónde existe mayor diversidad con AMOVA
 ┗ 📜 README.md      # Este documento
