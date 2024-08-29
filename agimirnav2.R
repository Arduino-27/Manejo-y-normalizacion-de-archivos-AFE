#28 de agosto 2024
#Autor Carmen Trejo
#Normalizacion de archivos AFE con agimirna procedimiento chatgpt

library(limma)
library(AgiMicroRna)
ls("package:AgiMicroRna")
help(package = "AgiMicroRna")

data_dir_miR <- "~/Documentos/archivos a normalizar/GSE208677_RAW_mirnas"


#Hacer un dataframe de la lista del nombre de los archivos, con tratamiento, repeticiones y el sujeto 
targets <- data.frame(
  FileName = c("GSM6364066_26_R1.txt", "GSM6364067_26_R2.txt", "GSM6364068_26_R3.txt", "GSM6364069_27_R1.txt", "GSM6364070_27_R2.txt", "GSM6364071_27_R3.txt", "GSM6364072_3_R1.txt",
               "GSM6364073_3_R2.txt", "GSM6364074_3_R3.txt", "GSM6364075_48_R1.txt", "GSM6364076_48_R2.txt", "GSM6364077_48_R3.txt", "GSM6364078_49_R2.txt", "GSM6364079_49_R3.txt",
               "GSM6364080_50_R1.txt", "GSM6364081_50_R2.txt", "GSM6364082_50_R3.txt", "GSM6364083_51_R1.txt", "GSM6364084_51_R2.txt", "GSM6364085_52_R1.txt", "GSM6364086_52_R2.txt",
               "GSM6364087_52_R3.txt", "GSM6364088_53_R1.txt", "GSM6364089_53_R2.txt", "GSM6364090_53_R3.txt", "GSM6364091_C1_R1.txt", "GSM6364092_C1_R2.txt", "GSM6364093_C1_R3.txt",
               "GSM6364094_C2_R1.txt","GSM6364095_C2_R2.txt" ),
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
      annotation = c( "ControlType", "ProbeName","GeneName")
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
    
    # $genes ("ControlType","ProbeName","GeneName") 
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
                       annotation = c( "ControlType", "ProbeName","SystematicName"),
                       verbose=TRUE)
dim(dd)
print(names(dd))

#plots 
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

#Construir una matriz de niveles de expresion de los mirnas (sin normalizar )
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

# Realizar la normalización de datos de microARN
ddNORM <- tgsNormalization(
  ddTGS,                # Datos a normalizar
  "quantile",           # Método de normalización
  makePLOTpre = FALSE,  # No generar gráficos antes de la normalización
  makePLOTpost = FALSE, # No generar gráficos después de la normalización
  targets,        # Información de los objetivos o meta datos
  verbose = TRUE        # Mostrar mensajes detallados
)





