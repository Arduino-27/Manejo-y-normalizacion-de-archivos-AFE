#28 de agosto 2024
#Autor Carmen Trejo
#Normalizacion de archivos AFE 

library(limma)
library(AgiMicroRna)
ls("package:AgiMicroRna")
help(package = "AgiMicroRna")

data_dir_miR <- "~/Documentos/archivos a normalizar/GSE208677_RAW_mirnas"


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
#No se estan guardando las imagenes 
qcPlots(dd,offset=5, 
        MeanSignal=TRUE,
        ProcessedSignal=TRUE,
        TotalProbeSignal=FALSE,
        TotalGeneSignal=FALSE,
        BGMedianSignal=FALSE,
        BGUsed=FALSE,
        targets)
png("qcPlots.png", width = 800, height = 600)
pdf("qcPlots.pdf", width = 8, height = 6)
dev.off()

#Construir una matriz de niveles de expresion de los mirnas (sin normalizar )
eset_mirna <- esetMicroRna(dd, targets)
expr_mirna <- exprs(eset_mirna)

annotation <- vroom::vroom("GSM6364088.txt", skip = 9, .name_repair = janitor::make_clean_names)

rownames(expr_mirna) <- annotation$systematic_name


is_mirna <- grep("^hsa", rownames(expr_mirna))


mirna_exprMat <- expr_mirna[is_mirna,]

# Realizar la normalización de datos de microARN
ddNORM=tgsNormalization(dd,"quantile",
                        makePLOTpre=FALSE,
                        makePLOTpost=TRUE,
                        targets=targets,
                        verbose=TRUE)

# Extraer los datos normalizados de la señal procesada
norm_expr <- ddNORM$procS

# Guardar la matriz como archivo CSV
write.csv(norm_expr, file = "normalized_expression.csv", row.names = TRUE)











