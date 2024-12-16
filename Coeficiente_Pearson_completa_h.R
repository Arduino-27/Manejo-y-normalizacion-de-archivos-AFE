# 20 de septiembre de 2024
# Carmen Trejo y Hugo Tovar
# Calcular correlacion entre dos matrices mRNAs y miRNAs
# usando paralelización

# radian


setwd("/scratch/home/hachepunto/carmen/coeficiente") # en fenix gpu

library(tictoc)
library(parallel)
library(Hmisc)

# Número de núcleos a usar
n_cores <- 50

# Leemos la tabla de datos - las columnas en A y B deben estar en el MISMO orden 
# y ser el mismo número

m=read.table(file="mensajeros_normalizados_anotados_mediana.csv", header = TRUE, row.names = 1, sep=",")
mi=read.table(file="matriz_mirnas_colapsada.csv", header = TRUE, row.names = 1,sep=",")


#      ____             _            __           
#     / __ \____ ______(_)__  ____  / /____  _____
#    / /_/ / __ `/ ___/ / _ \/ __ \/ __/ _ \/ ___/
#   / ____/ /_/ / /__/ /  __/ / / / /_/  __(__  ) 
#  /_/    \__,_/\___/_/\___/_/ /_/\__/\___/____/  
#                                                 


# pacientes
A <- as.matrix(m[,1:25])
B <- as.matrix(mi[,1:25])


# Haces una matriz vacía
cor_matrix <- matrix(NA, nrow = nrow(A), ncol = nrow(B))
pval_matrix <- matrix(NA, nrow = nrow(A), ncol = nrow(B))



# Función para calcular la correlación y el p-valor entre una fila de A y todas las filas de B
correlate_rows <- function(i) {
  row_A <- as.numeric(A[i, ])
  
  # Inicializar resultados temporales
  cor_results <- numeric(nrow(B))
  pval_results <- numeric(nrow(B))
  
  for (j in 1:nrow(B)) {
    row_B <- as.numeric(B[j, ])
    
    # Calcula la correlación de Pearson
    cor_results[j] <- cor(row_A, row_B, method = "pearson")
    
    # Calcula el p-valor usando rcorr
    pval_results[j] <- rcorr(cbind(row_A, row_B))$P[1, 2]
  }
  
  return(list(cor = cor_results, pval = pval_results))
}


# Medir el tiempo de ejecución
tic("Tiempo de ejecución")
# Usar mclapply para paralelizar el proceso de correlación entre las filas de A y B
resultados <- mclapply(1:nrow(A), correlate_rows, mc.cores = n_cores)
toc() # Termina el temporizador
# Tiempo de ejecución: 399.268 sec elapsed

# Almacenar los resultados en las matrices de correlaciones y p-valores
for (i in 1:nrow(A)) {
  cor_matrix[i, ] <- resultados[[i]]$cor
  pval_matrix[i, ] <- resultados[[i]]$pval
}


# Ponerle los nombres a columnas y renglones

rownames(cor_matrix) <-rownames(A)
colnames(cor_matrix) <-rownames(B)

rownames(pval_matrix) <-rownames(A)
colnames(pval_matrix) <-rownames(B)

# Escribir tus matrices en tablas CSV

write.csv(cor_matrix, file="resultados/P_pearson.csv")
write.csv(pval_matrix, file="resultados/P_pvalues.csv")



#     ______            __             __         
#    / ____/___  ____  / /__________  / /__  _____
#   / /   / __ \/ __ \/ __/ ___/ __ \/ / _ \/ ___/
#  / /___/ /_/ / / / / /_/ /  / /_/ / /  __(__  ) 
#  \____/\____/_/ /_/\__/_/   \____/_/\___/____/  
#                                                 

C <- as.matrix(m[,26:30])
D <- as.matrix(mi[,26:30])

cor_matrixC <- matrix(NA, nrow = nrow(C), ncol = nrow(D))
pval_matrixC <- matrix(NA, nrow = nrow(C), ncol = nrow(D))


correlate_rowsC <- function(i) {
  row_C <- as.numeric(C[i, ])
  
  # Inicializar resultados temporales
  cor_results <- numeric(nrow(D))
  pval_results <- numeric(nrow(D))
  
  for (j in 1:nrow(D)) {
    row_D <- as.numeric(D[j, ])
    
    # Calcula la correlación de Pearson
    cor_results[j] <- cor(row_C, row_D, method = "pearson")
    
    # Calcula el p-valor usando rcorr
    pval_results[j] <- rcorr(cbind(row_C, row_D))$P[1, 2]
  }
  
  return(list(cor = cor_results, pval = pval_results))
}

# Medir el tiempo de ejecución
tic("Tiempo de ejecución")
# Usar mclapply para paralelizar el proceso de correlación entre las filas de A y B
resultadosC <- mclapply(1:nrow(C), correlate_rowsC, mc.cores = n_cores)
toc() # Termina el temporizador
# Tiempo de ejecución: 190.171 sec elapsed

# Almacenar los resultados en las matrices de correlaciones y p-valores
for (i in 1:nrow(C)) {
  cor_matrixC[i, ] <- resultadosC[[i]]$cor
  pval_matrixC[i, ] <- resultadosC[[i]]$pval
}




# Ponerle los nombres a columnas y renglones

rownames(cor_matrixC) <-rownames(C)
colnames(cor_matrixC) <-rownames(D)

rownames(pval_matrixC) <-rownames(C)
colnames(pval_matrixC) <-rownames(D)

# Escribir tus matrices en tablas CSV

write.csv(cor_matrixC, file="resultados/C_pearson.csv")
write.csv(pval_matrixC, file="resultados/C_pvalues.csv")

# guardar sesión
save.image("resultados/Coeficiente_Pearson_completa_h.RData")
