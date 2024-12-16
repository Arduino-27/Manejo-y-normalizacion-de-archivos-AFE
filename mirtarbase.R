bd <- read.csv("miRTarBase.csv", header = TRUE, sep = ",")
mirnasF <- read.csv("mirna_mayor_menor_0_001.csv", header = TRUE, sep = ",")

# Filtrar los miRNAs en bd que coincidan con los de mirnasF
miRNAs_filtrados <- bd[bd$miRNA %in% mirnasF$miRNAs, ]
# Seleccionar solo las columnas miRNA y Target.Gene
miRNAs_filtradosOrdenados<- miRNAs_filtrados[, c("miRNA", "Target.Gene", "Experiments")]

# Cargar el paquete dplyr
library(dplyr)

# Eliminar duplicados para cada miRNA y su Target.Gene
miRNAs_targets_unicos <- miRNAs_filtradosOrdenados %>%
  distinct(miRNA, Target.Gene, .keep_all = TRUE)

#Guardar miRNAs_targets_unicos
write.csv(miRNAs_targets_unicos, "mirnasCoincidenMirtarbaseBDarray.csv", row.names = FALSE)

#Buscar coincidencias mirnas y genes con BSIF 
mirtarbase<- miRNAs_filtrados[, c("miRNA", "Target.Gene")]
mirtarbaseUnicos <- mirtarbase %>%
  distinct(miRNA, Target.Gene, .keep_all = TRUE)

bsif <- read.csv("BSIF_significativo_0.001.csv", header = TRUE, sep = ",")

# Reemplazar puntos por guiones en la columna de miRNA de bsif
bsif$miRNA <- gsub("\\.", "-", bsif$miRNA)

# Verificar los primeros valores después del reemplazo
head(bsif$miRNA)

coincidencias <- merge(
  mirtarbaseUnicos, 
  bsif, 
  by.x = c("miRNA", "Target.Gene"), 
  by.y = c("miRNA", "Gene")
)

colnames(mirtarbaseUnicos)
colnames(bsif)

write.csv(coincidencias, "coincidenciasBsifMirtarbase.csv", row.names = FALSE)

