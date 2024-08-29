#26 de agosto de 2024 
#Creado por Carmen Trejo y Dr. Hugo Tovar 
#Normalizacion de archivos AFE de mirnas de retinoblastoma
#Basado en AgiMicroRna de Pedro Lopez-Romero; April 30, 2024
setwd("/Users/hachepunto/Library/CloudStorage/OneDrive-Personal/datafolder/vponce/carmen")

#####################################################
library(limma)
library(AgiMicroRna)

#Cargar lista en txt con el nombre de las muestras y convertirlo en un dataframe 
targets <- vroom::vroom("targets.txt")
targets <- as.data.frame(targets)
row.names(targets) <- targets[[1]]
targets <- targets[-1]

dd <- read.agiMicroRna(targets,
columns=list(TGS="gTotalGeneSignal",
TPS="gTotalProbeSignal",
meanS="gMedianSignal",
procS="gProcessedSignal"),
other.columns=list(IsGeneDetected="gIsGeneDetected",
IsSaturated="gIsSaturated",
IsFeatNonUnifOF="gIsFeatNonUnifOL",
IsFeatPopnOL="gIsFeatPopnOL",
BGKmd="gBGMedianSignal"),
annotation = c( "ControlType", "ProbeName","SystematicName"),
verbose=TRUE)

dd=readMicroRnaAFE(targets,verbose=TRUE)

pdf("qcplots.pdf")
qcPlots(dd,offset=5,
MeanSignal=TRUE,
ProcessedSignal=FALSE,
TotalProbeSignal=FALSE,
TotalGeneSignal=FALSE,
BGMedianSignal=FALSE,
BGUsed=FALSE,
targets)
dev.off()

eset_mirna <- esetMicroRna(dd, targets)
expr_mirna <- exprs(eset_mirna)

annotation <- vroom::vroom("GSM6364088_53_R1.txt", skip = 9, .name_repair = janitor::make_clean_names)

rownames(expr_mirna) <- annotation$systematic_name


is_mirna <- grep("^hsa", rownames(expr_mirna))

mirna_exprMat <- expr_mirna[is_mirna,]


ddTGS=tgsMicroRna(dd,
half=TRUE,
makePLOT=FALSE,
verbose=FALSE)

####################################################









# Descargar un archivo GSM_miRNAs 
mirna_raw <- getGEO("GSE208677", GSEMatrix = TRUE)
# Extraer los datos de la muestra
mirna_exp <- exprs(mirna_raw[[1]])

# Mostrar las primeras filas de los datos
head(gsm_data)

# Verificar los datos
summary(gsm_data)

# Limpiar los datos (por ejemplo, eliminar filas con NA)
gsm_data_clean <- na.omit(gsm_data)

# Normalización por quantiles
norm_data <- normalizeBetweenArrays(gsm_data_clean, method = "quantile")

# Mostrar las primeras filas de los datos normalizados
head(norm_data)

# Visualización con boxplot antes y después de la normalización
par(mfrow=c(1,2))
boxplot(gsm_data_clean, main = "Antes de la normalización")
boxplot(norm_data, main = "Después de la normalización")

#Guarda los datos normalizados 
write.csv(norm_data, "datos_normalizados_miRNAs.csv", row.names = TRUE)


# Descargar un archivo GSM_mRNAs 
gsm <- getGEO("GSE208143", GSEMatrix = TRUE)
# Extraer los datos de la muestra
gsm_data <- exprs(gsm[[1]])

# Mostrar las primeras filas de los datos
head(gsm_data)

# Verificar los datos
summary(gsm_data)

# Limpiar los datos (por ejemplo, eliminar filas con NA)
gsm_data_clean <- na.omit(gsm_data)

# Normalización por quantiles
norm_data <- normalizeBetweenArrays(gsm_data_clean, method = "quantile")

# Mostrar las primeras filas de los datos normalizados
head(norm_data)

# Visualización con boxplot antes y después de la normalización
par(mfrow=c(1,2))
boxplot(gsm_data_clean, main = "Antes de la normalización")
boxplot(norm_data, main = "Después de la normalización")

#Guarda los datos normalizados 
write.csv(norm_data, "datos_normalizados_mRNAs.csv", row.names = TRUE)



