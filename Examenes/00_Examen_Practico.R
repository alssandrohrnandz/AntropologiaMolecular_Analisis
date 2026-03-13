# ==============================================================================
# EXAMEN PRÁCTICO DE R:CONCEPTOS CLAVE EN FILOGENIA
# Author: Alessandro Lopez-Hernandez
# Date: March 13th 2026
# ==============================================================================
# Correr a partir de aquí en R Studio con Ctrl+Enter línea por línea
# Instalación de paquetes
# En caso de que no lo tengan isntalado se ejecutarán en automatico
# TOTAL: 50 puntos

## CARGA DE PAQUETES ##
if (!require("adegenet")) install.packages("adegenet")
if (!require("ape")) install.packages("ape")
if (!require("phangorn")) install.packages("phangorn")
set.seed(5198)

library(ape)
library(phangorn)
library(adegenet)

# =============================================================================
# MODULO 1: PCA
# Instrucciones: ejecuta el codigo siguiente.
# Estos datos fueron extraidos del consorcio internacional
# "1000 Human Genomes Project" y tomaron representantes
# de los cinco continentes. Como investigadora, te tocará analizar
# unicamente 4 regiones.
# 1. Utiliza la siguiente función para cargar tus datos genéticos:

data(eHGDP)

# Después usa esta función que te permite explorar los objetos
# tipo .GENID para ver las regiones que existen

unique(eHGDP$other$popInfo$Region)

# Se te desplegara una lista de regiones que de donde provienen
# estos 1,1000 individuos. Escoge tres regiones de tu elección
# teniendo cuidado de mantener su nombre. Por ejemplo, si 
# decidiste usar a personas de la región "EAST_ASIA", usa el mismo
# nombre. Ejecuta la siguiente parte del codigo y escoge:

poblaciones_originales <- pop(eHGDP)
levels(poblaciones_originales) <- eHGDP$other$popInfo$Region
pop(eHGDP) <- poblaciones_originales

#Sustituye aqui por las regiones deseadas:
poblaciones_deseadas <- c("CENTRAL_SOUTH_ASIA", "MIDDLE_EAST", "AMERICA") 
#Filtramos:
genotipos_humanos <- eHGDP[pop(eHGDP) %in% poblaciones_deseadas, ]
pop(genotipos_humanos) <- droplevels(pop(genotipos_humanos))

#Generamos el PCA
obj_pca <- dudi.pca(tab(genotipos_humanos, NA.method="mean"), 
                    center = TRUE, scale = FALSE, scannf = FALSE, nf = 3)
#Obtenemos las coordenadas con la que graficaremos
coordenadas <- obj_pca$li
grupos <- pop(genotipos_humanos)

# s.class dibuja los puntos y el plot
s.class(dfxy = coordenadas, 
        fac = grupos, 
        xax = 1, yax = 2,  # Graficar PC1 vs PC2
        #Colores editables. Puedes usarlos en inglés o en sistema hexagesimal de colores
        col = c("firebrick", "dodgerblue", "forestgreen", "gold"),
        cpoint = 1.5)      
#Tu título
title("Estructura Poblacional Humana de 3 poblaciones")

# ==============================================================================
# Responde las siguientes preguntas
# 1. ¿Cómo salió tu PCA? Generar la imagen y entregarla (5 puntos)
# 2. ¿Qué tanta relación tienen tus poblaciones escogidas? Argumenta (5 puntos)
# 3. ¿Qué significa cada punto en el PCA? ¿Qué representa la distancia entre
#       los puntos? (5 puntos)
# 4. Ejecuta el siguiente codigo y responde ¿Cuánta varianza explica 
#       el Componente 1 y 2? Entregar la imagen resultante (10 puntos)

eigenvalores <- obj_pca$eig
porcentajes <- (eigenvalores / sum(eigenvalores)) * 100
head(round(porcentajes, 2), 10)
posiciones_barras <- barplot(porcentajes[1:10], 
                             main = "Varianza Explicada por cada PC (Scree Plot)",
                             xlab = "Componentes Principales (PCs)",
                             ylab = "Porcentaje de Varianza (%)",
                             names.arg = paste0("PC", 1:10),
                             col = "steelblue", #Colores editables
                             border = "black",
                             ylim = c(0, max(porcentajes) + 8), 
                             las = 1)
lines(x = posiciones_barras, 
      y = porcentajes[1:10], 
      type = "b", pch = 19, col = "red", lwd = 2)

text(x = posiciones_barras, 
     y = porcentajes[1:10], 
     labels = paste0(round(porcentajes[1:10], 1), "%"), # Redondea a 1 decimal y añade el "%"
     pos = 3,     # pos = 3 le indica a R que ponga el texto "arriba" del punto Y
     cex = 0.85,  # Tamaño de letra editable
     col = "black",
     font = 2)
# ===============================================================================
# CONCEPTOS CLAVE EN FILOGENIA
# ===============================================================================

# Cargamos nuestros datos
data("woodmouse")

# Calculamos la matriz de distancias genéticas usando el modelo K80
dist_mat <- dist.dna(woodmouse, model = "K80")
# Construimos el árbol Neighbor-Joining
nj_tree <- nj(dist_mat)
#Aqui pedimos estrictamente que genere un arbol enrazando al raton 
# llamado aqui como "No305", actuando como grupo externo. 
rooted_tree <- root(nj_tree, outgroup = "No305", resolve.root = TRUE)
# Graficamos ambos árboles para la comparación visual
par(mfrow = c(1, 2))
plot(nj_tree, main = "Árbol NJ Sin Enraizar (K80)", type = "unrooted")
plot(rooted_tree, main = "Árbol NJ Enraizado (K80)")

# 1. Generación de la imagen como entregable (5 puntos)
# 2. En este árbol observamos que se enraizó el sujeto "No305"
#       Observa detenidamente al sujeto "No1103" (No1103S en el otro arbol), 
#       ¿Cómo cambió su posición topológica (quién es su clado hermano) antes 
#       y después de aplicar el enraizamiento? ¿Por qué ocurre esta 
#       ilusión óptica en los árboles sin enraizar?
# 3. Ejecuta el siguiente código para generar un árbol de máxima
#       parsimonia. 10 puntos

woodmouse_phyDat <- as.phyDat(woodmouse)
arbol_inicial <- rtree(n = 15, tip.label = labels(woodmouse_phyDat))
tree_pars <- optim.parsimony(arbol_inicial, woodmouse_phyDat)

par(mfrow = c(1, 2))
plot(rooted_tree, main = "Árbol NJ Enraizado (K80)")
edgelabels(round(nj_tree$edge.length, 3), bg = "lightgreen", cex = 0.7)

plot(tree_pars, main= "Árbol de Máxima Parsimonia")
edgelabels(tree_pars$edge.length, bg = "lightblue", cex = 0.7)

# 4. El significado del largo de las ramas es totalmente distinto
#       entre árboles. ¿Cuál es la razón por la que el árbol NJ
#       tiene decimales? ¿Por qué en el de MP tenemos números enteros?
#       (5 puntos)
