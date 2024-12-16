#Script para leer, contabilizar y eliminar NAs de las columnas 
#Carmen Trejo 
#26 de septiembre 2024


#Leer Las matrices de controles y contabilizar los NA
datos <- read.csv("C_pearson.csv", header = TRUE, sep = ",")
columnas_con_na <- sum(apply(datos, 2, function(col) any(is.na(col))))

datos_pvalue_C <- read.csv("C_pvalues.csv", header = TRUE, sep = ",")
columnas_con_na_pvalue <- sum(apply(datos_pvalue_P, 2, function(col) any(is.na(col))))

#Leer Las matrices de pacientes y contabilizar los NA

datos_p <- read.csv("P_pearson.csv", header = TRUE, sep = ",")
columnas_con_na_p_pearson <- sum(apply(datos_p, 2, function(col) any(is.na(col))))

datos_pvalue_Pa <- read.csv("P_pvalues.csv", header = TRUE, sep = ",")
columnas_con_na_pvalue_p <- sum(apply(datos_pvalue_Pa, 2, function(col) any(is.na(col))))

mensajeros_matriz <- read.csv("mensajeros_normalizados_anotados_mediana.csv", header = TRUE, sep = ",")
mirnas_matriz <- read.csv("matriz_mirnas_colapsada.csv", header = TRUE, sep = ",")

#Eliminar columnas con NA 
# Eliminar filas con NA primero
datos_limpios_C_pearson <- na.omit(datos_p)

# Luego eliminar columnas con NA
# Crear una nueva matriz o data frame sin columnas con NA
datos_sin_na_columnas <- datos[, colSums(is.na(datos)) == 0]
datos_sin_na_columnas_p_c <- datos_pvalue_C[, colSums(is.na(datos_pvalue_C)) == 0]
datos_sin_na_columnas_pacientes <- datos_p[, colSums(is.na(datos_p)) == 0]
datos_sin_na_columnas_pvalue_pacientes <- datos_pvalue_Pa[, colSums(is.na(datos_pvalue_Pa)) == 0]

# Guardar la matriz en un archivo CSV
write.csv(datos_sin_na_columnas, file = "Coeficiente_pearson_controles.csv", row.names = FALSE)
write.csv(datos_sin_na_columnas_p_c, file = "pvalue_controles.csv", row.names = FALSE)
write.csv(datos_sin_na_columnas_pacientes, file = "pearson_pacientes.csv", row.names = FALSE)
write.csv(datos_sin_na_columnas_pvalue_pacientes, file = "pvalue_pacientes.csv", row.names = FALSE)

controles <- datos_sin_na_columnas 
pacientes <- datos_sin_na_columnas_pacientes

#Matriz n>0 y n<1 CONTROLES
sum(controles > 0 & controles < 1)
# Paso 1: Separar la primera columna (suponiendo que contiene los nombres de los genes)
nombres_genes <- controles[, 1]  # Almacena la primera columna

# Paso 2: Aplicar la condición para reemplazar con NA en el resto de las columnas
submatriz_rangoC <- controles[, -1]  # Excluye la primera columna (nombres de genes)
submatriz_rangoC[!(submatriz_rangoC > 0 & submatriz_rango < 1)] <- NA

# Paso 3: Mostrar la submatriz resultante
print(submatriz_rangoC)

# Paso 4: Pegar la columna de nombres de genes a la submatriz filtrada
submatriz_finalC <- cbind(nombres_genes, submatriz_rangoC)

# Paso 5: Mostrar la submatriz final
print(submatriz_finalC)

#Matriz n<0 y n>-1 CONTROLES
sum(controles < 0 & controles > -1)
# Paso 1: Separar la primera columna (suponiendo que contiene los nombres de los genes)
nombres_genes <- controles[, 1]  # Almacena la primera columna

# Paso 2: Aplicar la condición para reemplazar con NA en el resto de las columnas
submatriz_rangoCn <- controles[, -1]  # Excluye la primera columna (nombres de genes)
submatriz_rangoCn[!(submatriz_rangoCn < 0 & submatriz_rango > -1)] <- NA

# Paso 3: Mostrar la submatriz resultante
print(submatriz_rangoCn)

# Paso 4: Pegar la columna de nombres de genes a la submatriz filtrada
submatriz_finalCn <- cbind(nombres_genes, submatriz_rangoCn)

# Paso 5: Mostrar la submatriz final
print(submatriz_finalCn)

##################################
#Matriz n>0 y n<1 PACIENTES
sum(pacientes > 0 & pacientes < 1)
# Paso 1: Separar la primera columna (suponiendo que contiene los nombres de los genes)
nombres_genes <- pacientes[, 1]  # Almacena la primera columna

# Paso 2: Aplicar la condición para reemplazar con NA en el resto de las columnas
submatriz_rangoP <- pacientes[, -1]  # Excluye la primera columna (nombres de genes)
submatriz_rangoP[!(submatriz_rangoP > 0 & submatriz_rangoP < 1)] <- NA

# Paso 3: Mostrar la submatriz resultante
print(submatriz_rangoP)

# Paso 4: Pegar la columna de nombres de genes a la submatriz filtrada
submatriz_finalP <- cbind(nombres_genes, submatriz_rangoP)

# Paso 5: Mostrar la submatriz final
print(submatriz_finalP)

#Matriz n<0 y n>-1 Pacientes 
sum(pacientes < 0 & pacientes > -1)
# Paso 1: Separar la primera columna (suponiendo que contiene los nombres de los genes)
nombres_genes <- pacientes[, 1]  # Almacena la primera columna

# Paso 2: Aplicar la condición para reemplazar con NA en el resto de las columnas
submatriz_rangoPn <- pacientes[, -1]  # Excluye la primera columna (nombres de genes)
submatriz_rangoPn[!(submatriz_rangoPn < 0 & submatriz_rangoPn > -1)] <- NA

# Paso 3: Mostrar la submatriz resultante
print(submatriz_rangoPn)

# Paso 4: Pegar la columna de nombres de genes a la submatriz filtrada
submatriz_finalPn <- cbind(nombres_genes, submatriz_rangoPn)

# Paso 5: Mostrar la submatriz final
print(submatriz_finalPn)

#######################

#########Ejercicio coeficientes negativos n <-0.5 y n >-1
#Controles 
#Matriz n< -0.5 y n>-1 CONTROLES
sum(controles < -0.6 & controles > -0.9)
# Paso 1: Separar la primera columna (suponiendo que contiene los nombres de los genes)
nombres_genes <- controles[, 1]  # Almacena la primera columna

# Paso 2: Aplicar la condición para reemplazar con NA en el resto de las columnas
submatriz_rango_control <- controles[, -1]  # Excluye la primera columna (nombres de genes)
submatriz_rango_control[!(submatriz_rango_control < -0.6 & submatriz_rango_control > -0.9)] <- NA

# Paso 4: Pegar la columna de nombres de genes a la submatriz filtrada
submatriz_final_controles <- cbind(nombres_genes, submatriz_rango_control)

#5 de nov de 2024
#Contabilizar pvalues 
pvalues <- read.csv("pvalue_pacientes.csv", header = TRUE, sep = ",")

#p values <0.05
sum(pvalues < 0.05)
nombres_genes_pvalues <- pvalues[,1]
submatriz_rangoPn <- pvalues[, -1]  # Excluye la primera columna (nombres de genes)
submatriz_rangoPn[!(submatriz_rangoPn < 0.05)] <- NA
print(submatriz_rangoPn)
submatriz_finalPn_0.05 <- cbind(nombres_genes_pvalues, submatriz_rangoPn)
print(submatriz_finalPn_0.05)

#p values <0.01
sum(pvalues > 0.01)
nombres_genes_pvalues <- pvalues[,1]
submatriz_rangoPn <- pvalues[, -1]  # Excluye la primera columna (nombres de genes)
submatriz_rangoPn[!(submatriz_rangoPn < 0.01)] <- NA
print(submatriz_rangoPn)
submatriz_finalPn_0.01 <- cbind(nombres_genes_pvalues, submatriz_rangoPn)
print(submatriz_finalPn_0.01)

#p values <0.001
sum(pvalues > 0.001)
nombres_genes_pvalues <- pvalues[,1]
submatriz_rangoPn <- pvalues[, -1]  # Excluye la primera columna (nombres de genes)
submatriz_rangoPn[!(submatriz_rangoPn < 0.001)] <- NA
print(submatriz_rangoPn)
submatriz_finalPn_0.001 <- cbind(nombres_genes_pvalues, submatriz_rangoPn)
print(submatriz_finalPn_0.001)

str(submatriz_finalPn)
str(submatriz_finalPn_0.05)

#Creacion de B_SIF de pvalue 0.05 
#convertir en matriz
matrix_p <- as.matrix(submatriz_finalPn_0.05 [,-1])
rownames(matrix_p) <- submatriz_finalPn_0.05$nombres_genes_pvalues

matrix_c <- as.matrix(submatriz_finalPn[-1])
row.names(matrix_c) <- submatriz_finalPn$nombres_genes

#Valores negativos y significativos 
matrix_c_filtrada <- matrix_c
matrix_c_filtrada[is.na(matrix_p)] <- NA

#convertir a formato de 3 columnas mirna,gen,correlacion 
library(tidyr)
library(tidyverse)
# Convertir la matriz filtrada en un data.frame con nombres
B_long <- as.data.frame(matrix_c_filtrada) %>%
  # Cambiar a formato largo
  pivot_longer(
    cols = everything(),
    names_to = "miRNA",
    values_to = "correlation"
  ) %>%
  # Agregar la columna de genes (que corresponde al nombre de las filas de la matriz)
  mutate(Gene = rep(rownames(matrix_c_filtrada), ncol(matrix_c_filtrada))) %>%
  # Filtrar valores NA
  filter(!is.na(correlation))

# Reordenar columnas en formato SIF
B_SIF <- B_long[, c("miRNA", "Gene", "correlation")]

# Resultado
print(head(B_SIF))

library(dplyr)

# Ordenar el dataframe B_SIF por la columna 'correlation' en orden ascendente
B_SIF_ordenado <- B_SIF %>% arrange(correlation)
write.csv(B_SIF_ordenado, file = "BSIF_significativo_0.05.csv", row.names = FALSE)
##################################################################
#Creacion de B_SIF de pvalue 0.01 
#convertir en matriz
matrix_p <- as.matrix(submatriz_finalPn_0.01 [,-1])
rownames(matrix_p) <- submatriz_finalPn_0.01$nombres_genes_pvalues

matrix_c <- as.matrix(submatriz_finalPn[-1])
row.names(matrix_c) <- submatriz_finalPn$nombres_genes

#Valores negativos y significativos 
matrix_c_filtrada <- matrix_c
matrix_c_filtrada[is.na(matrix_p)] <- NA

#convertir a formato de 3 columnas mirna,gen,correlacion 
library(tidyr)
library(tidyverse)
# Convertir la matriz filtrada en un data.frame con nombres
B_long <- as.data.frame(matrix_c_filtrada) %>%
  # Cambiar a formato largo
  pivot_longer(
    cols = everything(),
    names_to = "miRNA",
    values_to = "correlation"
  ) %>%
  # Agregar la columna de genes (que corresponde al nombre de las filas de la matriz)
  mutate(Gene = rep(rownames(matrix_c_filtrada), ncol(matrix_c_filtrada))) %>%
  # Filtrar valores NA
  filter(!is.na(correlation))

# Reordenar columnas en formato SIF
B_SIF_2 <- B_long[, c("miRNA", "Gene", "correlation")]

# Resultado
print(head(B_SIF_2))

library(dplyr)

# Ordenar el dataframe B_SIF por la columna 'correlation' en orden ascendente
B_SIF_ordenado2 <- B_SIF_2 %>% arrange(correlation)
write.csv(B_SIF_ordenado2, file = "BSIF_significativo_0.01.csv", row.names = FALSE)
############################################
#Creacion de B_SIF de pvalue 0.001
#convertir en matriz
matrix_p <- as.matrix(submatriz_finalPn_0.001 [,-1])
rownames(matrix_p) <- submatriz_finalPn_0.001$nombres_genes_pvalues

matrix_c <- as.matrix(submatriz_finalPn[-1])
row.names(matrix_c) <- submatriz_finalPn$nombres_genes

#Valores negativos y significativos 
matrix_c_filtrada <- matrix_c
matrix_c_filtrada[is.na(matrix_p)] <- NA

#convertir a formato de 3 columnas mirna,gen,correlacion 
library(tidyr)
library(tidyverse)
# Convertir la matriz filtrada en un data.frame con nombres
B_long <- as.data.frame(matrix_c_filtrada) %>%
  # Cambiar a formato largo
  pivot_longer(
    cols = everything(),
    names_to = "miRNA",
    values_to = "correlation"
  ) %>%
  # Agregar la columna de genes (que corresponde al nombre de las filas de la matriz)
  mutate(Gene = rep(rownames(matrix_c_filtrada), ncol(matrix_c_filtrada))) %>%
  # Filtrar valores NA
  filter(!is.na(correlation))

# Reordenar columnas en formato SIF
B_SIF_3 <- B_long[, c("miRNA", "Gene", "correlation")]

# Resultado
print(head(B_SIF))

library(dplyr)

# Ordenar el dataframe B_SIF por la columna 'correlation' en orden ascendente
B_SIF_ordenado_3 <- B_SIF_3 %>% arrange(correlation)
write.csv(B_SIF_ordenado_3, file = "BSIF_significativo_0.001.csv", row.names = FALSE)

######################################
signiicativo_0.05 <- read_csv("BSIF_significativo_0.05.csv")
mirnas_unicos <- unique(signiicativo_0.05$miRNA)
numero_mirnas_unicos <- length(mirnas_unicos)
genes_unicos <- unique(signiicativo_0.05$Gene)
numero_genes_unicos <- length(genes_unicos)

# Dividir el dataframe en una lista de sub-dataframes por miRNA
miRNAs_organizados <- split(signiicativo_0.05,signiicativo_0.05$miRNA)
# Ejemplo: acceder al subconjunto de un miRNA específico
# Por ejemplo, para ver los datos del primer miRNA:
head(miRNAs_organizados[[1]])

# Si deseas guardar el resultado organizado en un archivo (opcional):
# Exportar cada sub-dataframe a un archivo separado
lapply(names(miRNAs_organizados), function(miRNA) {
  write.csv(miRNAs_organizados[[miRNA]], paste0("miRNA_", miRNA, ".csv"), row.names = FALSE)
})

write.csv(mirnas_unicos, file = "nombre_mirnas.csv")

# Contar cuántos genes están asociados a cada miRNA
conteo_genes <- aggregate(Gene ~ miRNA, data = B_SIF, FUN = length)
# Cambiar el nombre de las columnas para mayor claridad
colnames(conteo_genes) <- c("miRNAs", "num_genes")
write.csv(conteo_genes, "miRNAs_num_genes_0.05.csv", row.names = FALSE)

######################################
signiicativo_0.01 <- read.csv("BSIF_significativo_0.01.csv")
mirnas_unicos_0.01 <- unique(signiicativo_0.01$miRNA)
numero_mirnas_unicos_0.01 <- length(mirnas_unicos_0.01)
genes_unicos_0.01 <- unique(signiicativo_0.01$Gene)
numero_genes_unicos_0.01 <- length(genes_unicos_0.01)

# Dividir el dataframe en una lista de sub-dataframes por miRNA
miRNAs_organizados_0.01 <- split(signiicativo_0.01,signiicativo_0.01$miRNA)
# Ejemplo: acceder al subconjunto de un miRNA específico
# Por ejemplo, para ver los datos del primer miRNA:
head(miRNAs_organizados[[1]])

# Si deseas guardar el resultado organizado en un archivo (opcional):
# Exportar cada sub-dataframe a un archivo separado
lapply(names(miRNAs_organizados_0.01), function(miRNA) {
  write.csv(miRNAs_organizados_0.01[[miRNA]], paste0("miRNA_", miRNA, ".csv"), row.names = FALSE)
})

write.csv(mirnas_unicos_0.01, file = "nombre_mirnas.csv")

# Contar cuántos genes están asociados a cada miRNA
conteo_genes_0.01 <- aggregate(Gene ~ miRNA, data = B_SIF_2, FUN = length)
# Cambiar el nombre de las columnas para mayor claridad
colnames(conteo_genes_0.01) <- c("miRNAs", "num_genes")
write.csv(conteo_genes_0.01, "miRNAs_num_genes_0.01.csv", row.names = FALSE)

######################################
signiicativo_0.001 <- read.csv("BSIF_significativo_0.001.csv")
mirnas_unicos_0.001 <- unique(signiicativo_0.001$miRNA)
numero_mirnas_unicos_0.001 <- length(mirnas_unicos_0.001)
genes_unicos_0.001 <- unique(signiicativo_0.001$Gene)
numero_genes_unicos_0.001 <- length(genes_unicos_0.001)

# Dividir el dataframe en una lista de sub-dataframes por miRNA
miRNAs_organizados_0.001 <- split(signiicativo_0.001,signiicativo_0.001$miRNA)
# Ejemplo: acceder al subconjunto de un miRNA específico
# Por ejemplo, para ver los datos del primer miRNA:
head(miRNAs_organizados[[1]])

# Si deseas guardar el resultado organizado en un archivo (opcional):
# Exportar cada sub-dataframe a un archivo separado
lapply(names(miRNAs_organizados_0.001), function(miRNA) {
  write.csv(miRNAs_organizados_0.001[[miRNA]], paste0("miRNA_", miRNA, "miRNAs_num_genes_0.01.csv"), row.names = FALSE)
})

write.csv(mirnas_unicos_0.001, file = "nombre_mirnas.csv")

# Contar cuántos genes están asociados a cada miRNA
conteo_genes_0.001 <- aggregate(Gene ~ miRNA, data = B_SIF_3, FUN = length)
# Cambiar el nombre de las columnas para mayor claridad
colnames(conteo_genes_0.001) <- c("miRNAs", "num_genes")
write.csv(conteo_genes_0.001, "miRNAs_num_genes_0.001.csv", row.names = FALSE)
# Ordenar el dataframe de mayor a menor
df_ordenado_0.001 <- conteo_genes_0.001 %>%
  arrange(desc(num_genes))
write.csv(df_ordenado_0.001, "mayor_menor_0_001.csv", row.names = FALSE)

######################################
signiicativo_0.001 <- read.csv("BSIF_significativo_0.001.csv")
genes_unicos_0.001 <- unique(signiicativo_0.001$Gene)
numero_genes_unicos_0.001 <- length(genes_unicos_0.001)
# Dividir el dataframe en una lista de sub-dataframes por miRNA
genes_organizados_0.001 <- split(signiicativo_0.001,signiicativo_0.001$Gene)
# Ejemplo: acceder al subconjunto de un miRNA específico
# Por ejemplo, para ver los datos del primer miRNA:
head(genes_organizados_0.001[[1]])
# Si deseas guardar el resultado organizado en un archivo (opcional):
# Exportar cada sub-dataframe a un archivo separado
lapply(names(genes_organizados_0.001), function(Gene) {
  write.csv(genes_organizados_0.001[[Gene]], paste0("Gene_", Gene, "Gene_num_genes_0.01.csv"), row.names = FALSE)
})

write.csv(genes_unicos_0.001, file = "nombre_gene.csv")
# Contar cuántos genes están asociados a cada miRNA
conteo_mirnas_0.001 <- aggregate(miRNA ~ Gene, data = B_SIF_3, FUN = length)
# Cambiar el nombre de las columnas para mayor claridad
colnames(conteo_mirnas_0.001) <- c("Gene", "num_mirnas")
write.csv(conteo_mirnas_0.001, "genes_num_mirnas_0.001.csv", row.names = FALSE)
# Ordenar el dataframe de mayor a menor
df_ordenado.genes_0.001 <- conteo_mirnas_0.001 %>%
  arrange(desc(num_mirnas))
write.csv(df_ordenado.genes_0.001, "mayor_menor.genes_0_001.csv", row.names = FALSE)
