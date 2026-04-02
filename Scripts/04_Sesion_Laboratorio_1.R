## -----------------------------------------------------------------------------
options(warn = -1)  # oculta warnings

# Cargamos el paquete "ape"
library(ape)
# Descargamos la secuencia de interés utilizando su número de acceso
secuencia_1 <- read.GenBank("U37731.1")
# Visualizamos la secuencia descargada
print(secuencia_1)


## -----------------------------------------------------------------------------
# Descargamos varias secuencias utilizando sus números de acceso
secuencias <- read.GenBank(c("U37731.1", "U37730.1", "U37733.1", "U37737.1", "U37752.1",
                             "U37751.1", "U37739.1", "U37743.1")) 
# Visualizamos las secuencias descargadas
print(secuencias)


## -----------------------------------------------------------------------------
#| message: false
#| warning: false
    # Instalación de Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    # Instalación del paquete msa
BiocManager::install("msa")


## -----------------------------------------------------------------------------
#| message: false
#| warning: false
suppressPackageStartupMessages(library(msa))
suppressPackageStartupMessages(library(Biostrings))

# Descargamos las secuencias utilizando sus números de acceso
secuencias_crudas <- read.GenBank(c("U37731.1", "U37730.1", "U37733.1", "U37737.1", "U37752.1",
                             "U37751.1", "U37739.1", "U37743.1"))

# Este paso es escencial para convertir las secuencias al formato que el paquete msa requiere. El paquete ape devuelve las secuencias en un formato "DNAbin", por lo que es necesario transformarlas a un formato compatible con el paquete msa, como "DNAStringSet". Para esto, podemos usar la función as.character para transformar las secuencias en un formato compatible con el paquete msa.
secuencias_letras <- as.character(secuencias_crudas) # Transformamos las secuencias en un formato compatible con el paquete msa

# Esto convierte c("A", "T", "G") en "ATG" para cada individuo
secuencias_texto <- sapply(secuencias_letras, paste, collapse = "")

# Convertimos el texto limpio al formato que msa exige
secuencias_biostrings <- DNAStringSet(secuencias_texto)

# Alineamos las secuencias utilizando la función msa
alineamiento <- msa(secuencias_biostrings, method = "Muscle")

# Visualizamos el alineamiento
print(alineamiento)


## -----------------------------------------------------------------------------
# Guardamos el alineamiento en un archivo FASTA
alineamiento_texto <- as(alineamiento, "DNAStringSet")

# El archivo se va a guardar donde sea que tengas tu proyecto de R. Si quieres guardarlo en una carpeta específica, solo tienes que especificar la ruta completa en el argumento filepath. 
writeXStringSet(alineamiento_texto, filepath = "../Results/alineamiento.fasta")


## -----------------------------------------------------------------------------
# Transformamos el alineamiento al formato "DNAbin"
alineamiento_dnabin <- as.DNAbin(alineamiento) 

# Calculamos las distancias genéticas entre las secuencias utilizando la función dist.dna
distancias_geneticas <- dist.dna(alineamiento_dnabin, model = "K80") # Puedes cambiar el modelo de sustitución de nucleótidos si lo deseas

# Visualizamos las distancias genéticas
print(distancias_geneticas)


## -----------------------------------------------------------------------------
# Generamos un heatmap de las distancias genéticas
heatmap(as.matrix(distancias_geneticas), symm = TRUE, col = colorRampPalette(c("white", "red"))(100), margins = c(5, 5))


## -----------------------------------------------------------------------------
# Generamos un dendrograma de las distancias genéticas
dendrograma <- hclust(distancias_geneticas, method = "average")
plot(dendrograma, main = "Dendrograma de Distancias Genéticas", xlab = "Secuencias", ylab = "Distancia Genética")


## -----------------------------------------------------------------------------
# Cargamos el paquete "phangorn" para usar la función upgma
suppressPackageStartupMessages(library(phangorn))
# Construimos el árbol filogenético utilizando el método de UPGMA
tree_upgma <- upgma(distancias_geneticas)
# Visualizamos el árbol filogenético
plot(tree_upgma, main = "UPGMA")


## -----------------------------------------------------------------------------
# Construimos el árbol filogenético utilizando el método de Neighbor-Joining
tree_nj <- nj(distancias_geneticas)
# Visualizamos el árbol filogenético
plot(tree_nj, main = "Neighbor-Joining")


## -----------------------------------------------------------------------------
fun_nj <- function(x) nj(dist.dna(x, model = "K80"))
# Generamos un árbol filogenético con soporte estadístico utilizando el método de bootstrap

bs_values <- boot.phylo(tree_nj, alineamiento_dnabin, fun_nj, B = 100)# B es el número de réplicas de bootstrap
# Visualizamos el árbol filogenético con soporte estadístico
plot(tree_nj, main = "Árbol Final con Soporte de Bootstrap")
nodelabels(bs_values, adj = c(1.2, -0.5), frame = "n", cex = 0.8,
           col = ifelse(bs_values > 70, "darkgreen", "red"))  # Agregamos los valores de bootstrap a los nodos del árbol


## -----------------------------------------------------------------------------
# Guardamos el árbol filogenético generado con el método de UPGMA en formato New
write.tree(tree_upgma, file = "../Results/arbol_upgma.newick")
# Guardamos el árbol filogenético generado con el método de Neighbor-Joining en formato Newick
write.tree(tree_nj, file = "../Results/arbol_nj.newick")

