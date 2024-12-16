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

annotation <- vroom::vroom("GSM6338047_26_R2.txt", skip = 9, .name_repair = janitor::make_clean_names)
length (grep("^A_", annotation$probe_name))

#Boxplots (Comparación antes y después de la normalización)
# Abrir un dispositivo gráfico para PDF
pdf("boxplots_normalization_mensajeros.pdf", width = 8, height = 6)
par(mfrow=c(1,2))
boxplot(gsm_data_clean_mensajeros, main = "Antes de la normalización", outline=FALSE, las=2)
boxplot(mensajeros, main = "Después de la normalización", outline=FALSE, las=2)
dev.off()

#Guarda los datos normalizados 
write.csv(mensajeros, "datos_normalizados_mRNAs.csv", row.names = TRUE)

# Vector de nombres
vec_nombres <- rownames(mensajeros)

# Data frame con nombres y valores
df <- data.frame(
  nombres = annotation$probe_name,
  valores = annotation$systematic_name
)

# Convertir el vector a un data frame
df_vec <- data.frame(nombres = vec_nombres)

# Hacer el left join usando merge()
resultado <- merge(df_vec, df, by = "nombres", all.x = TRUE)
resultado_f <- resultado[which(!duplicated(resultado$nombres)),]
# annotation
gconvert <- gconvert(query=resultado_f$valores, target="HGNC", mthreshold=1, filter_na=FALSE)
class(gconvert[,1])<-"integer"
gconvert<- gconvert[sort.list(gconvert[,1]),]
genesymbol <- gconvert[,5]
#genesymbol<-ifelse(genesymbol=="N/A", as.vector(rownames(edata)), as.vector(genesymbol))

keep <- which (!is.na(genesymbol))
genesymbol <- genesymbol[keep]
mensajeros <- mensajeros[keep,]
dim(mensajeros)

###  Medians collapse  ###
colMedians=function(mat,na.rm=TRUE) {return(apply(mat,2,median,na.rm=na.rm))}
precolaps <- data.frame(genesymbol, mensajeros)
colapsed <- do.call(rbind, lapply(split(precolaps,precolaps[,1]), function(chunk) {colMedians(chunk[,-1])}))

#Guarda los datos normalizados, anotados y colapsados por mediana 
write.csv(colapsed, "mensajeros_normalizados_anotados_mediana.csv", row.names = TRUE)




