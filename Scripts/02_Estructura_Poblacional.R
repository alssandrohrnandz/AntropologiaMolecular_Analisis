# ==============================================================================
# ANÁLISIS DE VARIANZA MOLECULAR (AMOVA)
# Author: Alessandro Lopez-Hernandez
# Date: March 15th 2026
# Correr a partir de aquí en R Studio con Ctrl+Enter línea por línea
# Instalación de paquetes
# En caso de que no lo tengan isntalado se ejecutarán en automatico
# ==============================================================================
# Cargar el paquete poppr (Especializado en genética de poblaciones)
# install.packages("poppr")

library(poppr)
library(adegenet)
# Será necesarioc argar nuevamente nuestros datos del eHHDP
data(eHGDP)
poblaciones_originales <- pop(eHGDP)
levels(poblaciones_originales) <- eHGDP$other$popInfo$Region
pop(eHGDP) <- poblaciones_originales
# En esta ocasion usaremos solo 6 poblaciones, tres de las cuales
# deben ser asiaticas
poblaciones_deseadas <- c("CENTRAL_SOUTH_ASIA", "MIDDLE_EAST", "CENTRAL_SOUTH_ASIA", "AMERICA", "OCEANIA", "AFRICA") 
genotipos_humanos <- eHGDP[pop(eHGDP) %in% poblaciones_deseadas, ]
pop(genotipos_humanos) <- droplevels(pop(genotipos_humanos))

# Definir la jerarquía (Estratos)
# poppr necesita saber explícitamente cómo están agrupados los datos.
# Creamos un "estrato" llamado Continente basado en la población activa.
strata(genotipos_humanos) <- data.frame(Continente = pop(genotipos_humanos))

# Ejecutar el AMOVA
# Le pedimos que analice la varianza en función del estrato "~Continente".
# Nota: En procesadores estándar esto puede tardar un minuto. 
amova_resultado <- poppr.amova(genotipos_humanos, ~Continente, 
                               missing = "ignore", # Ignora alelos faltantes
                               method = "ade4")    # Motor matemático

# Ver el resultado estadístico en consola
cat("\n--- Tabla de Resultados AMOVA ---\n")
print(amova_resultado)
# Extraer los porcentajes de varianza para graficar
# El objeto guarda los componentes de varianza en $varcomp
varianzas_crudas <- amova_resultado$componentsofcovariance$Sigma
var_entre <- varianzas_crudas[1]
var_dentro <- sum(varianzas_crudas[-1])
porcentajes_finales <- c(var_entre, var_dentro) / sum(varianzas_crudas) * 100

# ------------------------------------------------------------------------------
# VISUALIZACIÓN: GRÁFICO DE PASTEL DEL AMOVA
# ------------------------------------------------------------------------------

# Nombres legibles para el gráfico
etiquetas <- c("Entre Continentes", "Dentro de Continentes e Individuos")

png(filename = "Results/02_Estructura_Poblacional/AMOVA_pastel.png", width = 900, height = 900, res = 120)

pie(porcentajes_finales, 
    labels = paste0(etiquetas, "\n", round(porcentajes_finales, 1), "%"), 
    col = c("firebrick", "dodgerblue"), 
    border = "white",
    cex = 0.9,
    main = "Partición de la Varianza Genética Humana")

dev.off()
