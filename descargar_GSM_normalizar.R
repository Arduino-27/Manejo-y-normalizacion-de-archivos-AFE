#27 de julio 2024
#Creado por Carmen Trejo 
#Descarga archivos GSM y normaliza entre microarreglos 
setwd("~/Documentos")
#Instalar paquetes 
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")
BiocManager::install("limma")

# Verificar y cargar el paquete
if (require(BiocManager, quietly = TRUE)) {
  cat("El paquete está instalado y se ha cargado correctamente.\n")
} else {
  cat("El paquete no está instalado.\n")
}

#cargar los paquetes
library(GEOquery)
library(limma)

# Verificar y cargar el paquete
if (require(GEOquery, quietly = TRUE)) {
  cat("El paquete está instalado y se ha cargado correctamente.\n")
} else {
  cat("El paquete no está instalado.\n")
}

if (require(limma, quietly = TRUE)) {
  cat("El paquete está instalado y se ha cargado correctamente.\n")
} else {
  cat("El paquete no está instalado.\n")
}

# Descargar archivo GSM_miRNAs de RB 
gsm_mirna <- getGEO("GSE208677", GSEMatrix = TRUE)
# Extraer los datos de la muestra
gsm_data_mirna <- exprs(gsm_mirna[[1]])

# Mostrar las primeras filas de los datos
head(gsm_data_mirna)

# Verificar los datos
summary(gsm_data_mirna)

# Limpiar los datos (por ejemplo, eliminar filas con NA)
gsm_data_clean_mirna <- na.omit(gsm_data_mirna)

# Normalización por quantiles
norm_data_mirna <- normalizeBetweenArrays(gsm_data_clean_mirna, method = "quantile")

# Mostrar las primeras filas de los datos normalizados
head(norm_data_mirna)
mirnas <- norm_data_mirna

# Visualización con boxplot antes y después de la normalización
par(mfrow=c(1,2))
boxplot(gsm_data_clean_mirna, main = "Antes de la normalización")
boxplot(mirnas, main = "Después de la normalización")

#Guarda los datos normalizados 
write.csv(mirnas, "datos_normalizados_miRNAs.csv", row.names = TRUE)


# Descargar un archivo GSM_mRNAs RNAs MENSAJEROS 
gsm_mensajeros <- getGEO("GSE208143", GSEMatrix = TRUE)
# Extraer los datos de la muestra
gsm_data_mensajeros <- exprs(gsm_mensajeros[[1]])

# Mostrar las primeras filas de los datos
head(gsm_data_mensajeros)

# Verificar los datos
summary(gsm_data_mensajeros)

# Limpiar los datos (por ejemplo, eliminar filas con NA)
gsm_data_clean_mensajeros <- na.omit(gsm_data_mensajeros)

# Normalización por quantiles
norm_data_mensajeros <- normalizeBetweenArrays(gsm_data_clean_mensajeros, method = "quantile")

# Mostrar las primeras filas de los datos normalizados
head(norm_data_mensajeros)
mensajeros <- norm_data_mensajeros
# Visualización con boxplot antes y después de la normalización
par(mfrow=c(1,2))
boxplot(gsm_data_clean_mensajeros, main = "Antes de la normalización")
boxplot(mensajeros, main = "Después de la normalización")

#Guarda los datos normalizados 
write.csv(mensajeros, "datos_normalizados_mRNAs.csv", row.names = TRUE)



