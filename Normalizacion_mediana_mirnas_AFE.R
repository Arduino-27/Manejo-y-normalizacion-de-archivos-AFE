#28 de agosto 2024
#Autor Carmen Trejo Y Dr. Hugo Tovar Torres
#Normalizacion de archivos AFE 
#Los archivos AFE descargados de GEO deben estar en la misma carpeta que este script

library(limma)
library(AgiMicroRna)
ls("package:AgiMicroRna")
help(package = "AgiMicroRna")
browseVignettes(package = "AgiMicroRna")

#Hacer un dataframe de la lista del nombre de los archivos, con tratamiento, repeticiones y el sujeto 
targets <- data.frame(
  FileName = c("GSM6364066.txt", "GSM6364067.txt", "GSM6364068.txt", "GSM6364069.txt", "GSM6364070.txt", "GSM6364071.txt", "GSM6364072.txt",
               "GSM6364073.txt", "GSM6364074.txt", "GSM6364075.txt", "GSM6364076.txt", "GSM6364077.txt", "GSM6364078.txt", "GSM6364079.txt",
               "GSM6364080.txt", "GSM6364081.txt", "GSM6364082.txt", "GSM6364083.txt", "GSM6364084.txt", "GSM6364085.txt", "GSM6364086.txt",
               "GSM6364087.txt", "GSM6364088.txt", "GSM6364089.txt", "GSM6364090.txt", "GSM6364091.txt", "GSM6364092.txt", "GSM6364093.txt",
               "GSM6364094.txt","GSM6364095.txt" ),
  Treatment  = rep("A",30),
  GErep = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 2, 3, 1, 2, 3, 1, 2, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2),
  Subject = c(26, 26, 26, 27, 27, 27, 3, 3, 3, 48, 48, 48, 49, 49, 50, 50, 50, 51, 51, 52, 52, 52, 53, 53, 53, 1, 1, 1, 2, 2)
)

#Crear funcion tomado de https://rdrr.io/bioc/AgiMicroRna/src/R/read.agiMicroRna.R#sym-%60read.agiMicroRna%60

`read.agiMicroRna` <-
  function(targets,columns=NULL,other.columns = NULL,annotation=NULL,
           exp.names=NULL,verbose=FALSE){
    
    ## code adapted from read.maimages function from limma package 
    ## it creates the uRNAList class: see AgiMicroRNA-classes.R for 
    ## the definition of this class and associated methods  
    
    filess=as.character(targets$FileName)
    
    if (is.null(filess)) {
      stop("targets frame doesn't contain FileName column")
    }
    
    if (is.null(columns)){
      columns=list(TGS="gTotalGeneSignal",
                   TPS="gTotalProbeSignal",
                   meanS="gMeanSignal",
                   procS="gProcessedSignal")
    }
    
    if (!is.list(columns)){
      stop("columns must be a list")
    }
    
    if (is.null(other.columns)){
      other.columns=list(IsGeneDetected="gIsGeneDetected",
                         IsSaturated="gIsSaturated",
                         IsFeatNonUnifOF="gIsFeatNonUnifOL",
                         IsFeatPopnOL="gIsFeatPopnOL",
                         BGKmd="gBGMedianSignal",
                         BGKus="gBGUsed")
    } 
    
    if (!is.list(other.columns)){
      stop("other.columns must be a list")
    }
    
    if (is.null(exp.names)){
      annotation = c("ControlType", "ProbeName", "SystematicName")
    }
    
    if (is.null(exp.names)){
      exp.names = rownames(targets) 
    }
    
    
    cnames <- names(columns)
    required.col = unique(c(annotation, unlist(columns), unlist(other.columns))) 
    
    
    obj = read.columns(filess[1], required.col, text.to.search="",skip = 9, 
                       sep = "\t",quote="\"", stringsAsFactors = FALSE,flush=TRUE)
    narray = length(filess)
    ngenes = dim(obj)[1]
    
    # ngenes = length(scan(filess[1],skip=10,what="integer",flush=TRUE,quiet=TRUE))
    Y = matrix(NA, nrow = ngenes, ncol = narray)
    colnames(Y) =  exp.names 
    
    # R, B, Rb, Gb 
    Newagi = columns
    for (a in cnames){
      Newagi[[a]] <- Y
    }
    
    # targets 
    Newagi$targets = data.frame(targets$FileName)
    rownames(Newagi$targets) = exp.names
    colnames(Newagi$targets) = "FileName"
    
    # $genes ("ControlType", "ProbeName", "SystematicName") 
    j <- match(annotation, colnames(obj), 0)
    if (any(j > 0)){
      Newagi$genes <- data.frame(obj[, j, drop = FALSE], check.names = FALSE)
    }
    # $other 
    other.columns <- as.character(other.columns)
    j <- match(other.columns, colnames(obj), 0)
    if (any(j > 0)) {
      other.columns <- colnames(obj)[j]
      Newagi$other = list()
      for (j in other.columns) Newagi$other[[j]] <- Y
    }
    
    
    for(n in 1:narray) {
      if(verbose){
        cat("reading file ",n," - ",filess[n],"\n")
      }
      if(n > 1){
        obj = read.columns(filess[n], required.col, text.to.search="",skip = 9, 
                           sep = "\t",quote="\"", stringsAsFactors = FALSE,flush=TRUE)
      }
      
      for (a in cnames){
        Newagi[[a]][, n] <- obj[, columns[[a]]]    # R, G, Rb, Gb
      }
      
      for (j in other.columns) {
        Newagi$other[[j]][, n] <- obj[, j]
      }    
      
    } ## for  
    
    # defined in AgiMicroRNA-classes.R USING classes.R limma FILE 
    new("uRNAList", Newagi)
    
  } ## end function 

# Leer los archivos AFE usando la función definida
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
                       annotation = c( "ControlType", "ProbeName", "SystematicName"),
                       verbose=TRUE)
dim(dd)
print(names(dd))
#plots 
#Guardar plots, cada plot debe ser cerrado para poder obtener el pdf 
png("qcPlots.png", width = 800, height = 600)
pdf("qcPlots.pdf", width = 8, height = 6)
qcPlots(dd,offset=5, 
        MeanSignal=TRUE,
        ProcessedSignal=TRUE,
        TotalProbeSignal=FALSE,
        TotalGeneSignal=FALSE,
        BGMedianSignal=FALSE,
        BGUsed=FALSE,
        targets)

dev.off()

#Construir una matriz de niveles de expresion de los mirnas (sin normalizar )
eset_mirna <- esetMicroRna(dd, targets)
expr_mirna <- exprs(eset_mirna)

annotation <- vroom::vroom("GSM6364088.txt", skip = 9, .name_repair = janitor::make_clean_names)

rownames(expr_mirna) <- annotation$systematic_name


is_mirna <- grep("^hsa", rownames(expr_mirna))


mirna_exprMat <- expr_mirna[is_mirna,]

# Realizar la normalización de datos por quantiles entre microarreglos, guardar plots
# La normalizacion estara escalada a log 2
png("ddNORM.png", width = 800, height = 600)
#pdf("ddNORM.pdf", width = 8, height = 6)

ddNORM=tgsNormalization(dd,"quantile",
                        makePLOTpre=TRUE,
                        makePLOTpost=TRUE,
                        targets=targets,
                        verbose=TRUE)
if (dev.cur() > 1) dev.off()  # Cierra el dispositivo gráfico PNG si está abierto

# Extraer los datos normalizados de la señal procesada
norm_expr <- ddNORM$procS
norm_expr_names <- ddNORM$genes
combined_df <- cbind(norm_expr, norm_expr_names)

# Guardar la matriz como archivo CSV
write.csv(norm_expr, file = "normalized_expression.csv", row.names = TRUE)
write.csv(combined_df, file = "norm_expr_quantil_total.csv")

#Obtener la matriz de normalizacion por quantiles con el nombre de los mirnas 
rownames(norm_expr) <- annotation$systematic_name
is_mirna <- grep("^hsa", rownames(norm_expr))
mirna_matriz <- norm_expr[is_mirna,]
write.csv(mirna_matriz, file = "norm_expr_quantil_mirnas.csv")

#Normalizacion por preprocessCore
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  install.packages("preprocessCore")
}

library(preprocessCore)

mat_normalized <- normalize.quantiles(mirna_exprMat)
rownames(mat_normalized) <- rownames(mirna_exprMat)
colnames(mat_normalized)<- colnames(mirna_exprMat)

#Aplicar mediana para colapsar mirnas y tener un mirna un nivel de expresion
colMedians=function(mat,na.rm=TRUE){
  return(apply(mat,2,median,na.rm=na.rm))}

precolaps <- data.frame(rownames(mat_normalized), mat_normalized)

precolaps[1:10,1:10]

colapsed <- do.call(rbind, lapply(split(precolaps,precolaps[,1]), function(chunk) {colMedians(chunk[,-1])}))

mirna_matrix_colapsed <- as.matrix(colapsed)
write.csv(mirna_matrix_colapsed, file = "matriz_mirnas_colapsada.csv")

#Boxplots (Comparación antes y después de la normalización)
# Abrir un dispositivo gráfico para PDF
pdf("boxplots_normalization.pdf", width = 8, height = 6)
par(mfrow=c(1,2))
boxplot(mirna_exprMat, main = "Antes de la normalización", outline=FALSE, las=2)
boxplot(mat_normalized, main = "Después de la normalización", outline=FALSE, las=2)
dev.off()

#Gráficos de Densidad (Distribuciones de expresión)
plot(density(as.vector(mirna_exprMat)), col = "red", main = "Densidad de microRNAs", xlim = c(0, 16))
lines(density(as.vector(mat_normalized)), col = "blue")
legend("topright", legend = c("Antes", "Después"), col = c("red", "blue"), lty = 1)












