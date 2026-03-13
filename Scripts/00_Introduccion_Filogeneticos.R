# ==============================================================================
# TALLER DE ANÁLISIS GENÉTICO: CONCEPTOS CLAVE EN FILOGENIA
# Author: Alessandro Lopez-Hernandez
# Date: March 4th 2026
# ==============================================================================
# Correr a partir de aquí en R Studio con Ctrl+Enter línea por línea
# Instalación de paquetes
# En caso de que no lo tengan isntalado se ejecutarán en automatico
if (!require("ape")) install.packages("ape")
if (!require("phangorn")) install.packages("phangorn")

library(ape) #Ejecutar para cargar
library(phangorn) #Ejecutar para cargar

# Cargamos datos de ejemplo (15 secuencias de ADN mitocondrial)
set.seed(4598)
data(woodmouse)

# ==============================================================================
# MODELOS DE SUSTITUCIÓN 
# ==============================================================================
# Calculamos las distancias genéticas entre las secuencias de mtDNA
# usando dos modelos de sustitución distintos. Por default, solo tenemos
# K80 (Kimura 80) mejor conocido como "Kimura-2 parameters distance" (Kimura M., 1980).
# En esta práctica, averigua el todos los modelos existentes que se aplica en dist.dna(). 
# Aquí deberan sustituir "K80" por otros modelos. Modificarlo. 
dist_jc  <- dist.dna(woodmouse, model = "JC69")
dist_gtr <- dist.dna(woodmouse, model = "K80")

# ENTREGABLE 1: Comparación de distancias
par(mfrow = c(1, 2))
plot(dist_jc, dist_gtr, main = "Corrección de Distancias",
     # Cambiar de acuerdo al modelo escogido en vez de K80 
     # El resto dejarlo como está
     xlab = "Distancia Simple (JC69)", ylab = "Distancia Compleja (K80)", 
     pch = 19, col = rgb(0.2, 0.4, 0.6, 0.5))
abline(0, 1, col = "red", lwd = 2) 
# Nota: Los puntos sobre la línea roja indican donde el modelo simple subestima la evolución.

# ==============================================================================
# ALGORITMOS DE CONSTRUCCIÓN (NJ vs UPGMA)
# ==============================================================================
# Analizar cuidadosamente el algoritmo de construcción de árboles y las diferencias
# entre NJ y UPGMA
tree_nj    <- nj(dist_gtr)    # Neighbor-Joining (No asume reloj molecular)
tree_upgma <- upgma(dist_gtr) # UPGMA (Asume reloj molecular estricto)

# ENTREGABLE 2: Comparación de Topologías
par(mfrow = c(1, 2))
plot(tree_nj, main = "Neighbor-Joining (NJ)", sub = "Sin reloj molecular")
plot(tree_upgma, main = "UPGMA", sub = "Con reloj molecular")

# ==============================================================================
# MÉTRICAS DE SOPORTE (BOOTSTRAP)
# ==============================================================================

# Función para que el bootstrap sepa qué algoritmo usar (NJ en este caso)
# En esta parte se debe sustituir K80 por alguno de los modelos disponibles
# que pudieron averiguar
fun_nj <- function(x) nj(dist.dna(x, model = "K80"))

# Corremos 100 réplicas de Bootstrap (puede tardar unos segundos)
set.seed(123) # Para que todos tengan el mismo resultado
bs_values <- boot.phylo(tree_nj, woodmouse, fun_nj, B = 100)

# ENTREGABLE 3: Árbol con valores de confianza
plot(tree_nj, main = "Árbol Final con Soporte de Bootstrap")
# Añadimos los valores a los nodos (solo si son mayores a 50%)
nodelabels(bs_values, adj = c(1.2, -0.5), frame = "n", cex = 0.8,
           col = ifelse(bs_values > 70, "darkgreen", "red"))

# ==============================================================================
# FINAL DEL SCRIPT
# ==============================================================================