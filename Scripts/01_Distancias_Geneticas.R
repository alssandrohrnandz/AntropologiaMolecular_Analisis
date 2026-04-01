# ==============================================================================
# CÁLCULO DE DISTANCIAS GENÉTICAS Y MAPA DE CALOR (Fst)
# Author: Alessandro Lopez-Hernandez
# Date: March 15th 2026
# Correr a partir de aquí en R Studio con Ctrl+Enter línea por línea
# Instalación de paquetes
# En caso de que no lo tengan isntalado se ejecutarán en automatico
# ==============================================================================
if(!require("hierfstat")) install.packages("hierfstat")
if (!require("adegenet")) install.packages("adegenet")

library(hierfstat)
library(adegenet)

# hierfstat usa una estructura de datos ligeramente distinta a adegenet.
# Usamos esta función puente para transformar nuestro objeto 'genind'.
# Usaremos los mismos datos que 'adegenet' de eHGDP
data(eHGDP)
poblaciones_originales <- pop(eHGDP)
levels(poblaciones_originales) <- eHGDP$other$popInfo$Region
pop(eHGDP) <- poblaciones_originales
# TAREA 1: Probar con distintas poblaciones (máximo 3)
poblaciones_deseadas <- c("CENTRAL_SOUTH_ASIA", "MIDDLE_EAST", "AMERICA", "OCEANIA","CENTRAL_SOUTH_ASIA", "AFRICA") 
genotipos_humanos <- eHGDP
# Descomentar la siguiente línean para trabajar
# genotipos_humanos <- eHGDP[pop(eHGDP) %in% poblaciones_deseadas, ]
pop(genotipos_humanos) <- droplevels(pop(genotipos_humanos))

datos_hier <- genind2hierfstat(genotipos_humanos)

# Cálculo de la Matriz Fst por pares (Pairwise Fst)
# Utilizamos el estimador clásico de Weir & Cockerham (1984)
# NOTA: Al tener 377 marcadores y 1000 individuos, esto puede tardar unos 5 minutos
# dependiendo del poder de procesamiento de tu computadora
matriz_fst <- pairwise.WCfst(datos_hier)

# Exploraremos los datos crudos en la consola
# Redondeamos a 3 decimales para no saturar la pantalla
cat("\n--- Matriz de Fst por Pares ---\n")
print(round(matriz_fst, 3))

# ------------------------------------------------------------------------------
# VISUALIZACIÓN: MAPA DE CALOR (HEATMAP) EN R BASE
# ------------------------------------------------------------------------------

# Paso A: Limpieza de la diagonal
# La diagonal compara a la población consigo misma (ej. AFRICA vs AFRICA), 
# por lo que el cálculo arroja 'NA' (Not Available). Los cambiamos a 0 para graficar.
matriz_fst[is.na(matriz_fst)] <- 0

# Paso B: Crear una paleta de colores
# De blanco (0 = idénticos) pasando por naranja, a rojo oscuro (1 = muy diferentes)
# PUEDES PERSONALIZAR LOS COLORES pero te sugiero un máximo de 3 para no
# saturar la pantalla
paleta_colores <- colorRampPalette(c("white","gold","firebrick"))(100)

# Paso C: Generar el Gráfico
# Generamos primero el archivo en blanco para colocar nuestro plo después
png(filename = "Results/01_Distancias_Geneticas/HeatMap_Fst.png", width = 900, height = 900, res = 120)
# El plot se generará en el lienzo en blanco
heatmap(matriz_fst, 
        Rowv = NA, Colv = NA,      # Apaga los dendrogramas para mantenerlo simple
        col = paleta_colores,      # Aplica nuestra paleta
        scale = "none",            # Evita que R normalice los datos por fila
        main = "Mapa de Calor de Diferenciación Genética (Fst)")

# "Apagamos" el lienzo y podremos visualizar el resultado en la carpeta de destino
dev.off()
