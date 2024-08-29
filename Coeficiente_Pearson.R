#27 de junio 2024
#Creado por Dr. Enrique Lemus 
#Calcular correlacion entre dos matrices mRNAs y miRNAs del proyecto GSE208677 y GSE208143

setwd("~/Documentos")
miRNA <- read.csv("datos_normalizados_mirnas.csv") 
mRNA <- read.csv("datos_normalizados_mRNAs.csv")

# Leemos la tabla de datos - las columnas en A y B deben estar en el MISMO orden 
# y ser el mismo número

A=read.table(file="datos_normalizados_mRNAs.csv", header = TRUE, row.names = 1, sep=",")
B=read.table(file="datos_normalizados_mirnas.csv", header = TRUE, row.names = 1,sep=",")

# Haces una matriz vacía
cor_matrix <- matrix(NA, nrow = nrow(A), ncol = nrow(B))
pval_matrix <- matrix(NA, nrow = nrow(A), ncol = nrow(B))

# Calculas las correlaciones
for (i in 1:nrow(A)) {
  for (j in 1:nrow(B)) {
    row_A <- as.numeric(A[i,])
    row_B <- as.numeric(B[j,])
    
    #Calculas la correlación 
    result <- cor.test(row_A, row_B)
    
    #Lo pones en una matriz
    cor_matrix[i, j] <- result$estimate
    
    pval_matrix[i,j] <-result$p.value
  }
}
# Ponerle los nombres a columnas y renglones

rownames(cor_matrix) <-rownames(A)
colnames(cor_matrix) <-rownames(B)

rownames(pval_matrix) <-rownames(A)
colnames(pval_matrix) <-rownames(B)

# Escribir tus matrices en tablas CSV

write.csv(cor_matrix, file="pearson.csv")
write.csv(pval_matrix, file="pvalues.csv")

