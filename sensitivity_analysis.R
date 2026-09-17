req_packages <- c("AnnotationDbi", "DESeq2", "devtools", "doParallel", 
                  "dynamicTreeCut", "edgeR", "flashClust", "foreach", 
                  "ggdendro", "ggrepel", "ggplot2", "igraph",  
                  "limma", "MODA", "openxlsx", "org.Hs.eg.db", 
                  "plyr", "preprocessCore", "Rcpp", "reticulate", 
                  "RSQLite", "STRINGdb")

# Load required packages
for (i in 1:length(req_packages)) {
  library(req_packages[i], character.only = TRUE)
}

library("matrixStats")
library(EnhancedVolcano)
library(minfi)
library(DOSE)
library(clusterProfiler)

library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

library("RColorBrewer")
library("minfiData")
library("Gviz")
library("stringr")
library("mCSEA")

set.seed(777)

main_dir <- "COVUM/"
input_dir <- "Data/Normalized/"
output_dir <- "DMCs/"

####################################################################
# Input                                                            #
####################################################################

# Input - metadata
load(paste0(main_dir, input_dir, "beta_pd.RData"))
beta_pd <- beta_pd[(beta_pd$Timepoint == 3) | (beta_pd$Timepoint == 12)
                   & (beta_pd$Group == "PCC+"), ]
beta_pd$Sample_Name <- paste0('X', beta_pd$Sentrix_ID, '_', beta_pd$Sentrix_Position)

# Input - preprocessed beta-values
load(paste0(main_dir, input_dir, "beta_enmix_filt_combat.RData"))
betaVal <- beta_combat[,paste0('X', beta_pd$Sentrix_ID, '_', beta_pd$Sentrix_Position)]

# Input - EPIDISH
epidish_cellprop <- as.data.frame(read.csv(paste0(main_dir, input_dir, "beta_epidish_prop.csv")))

#####################################################################
# Fold change: T2 vs T1 - sample removal                            #
#####################################################################

study.ids <- unique(beta_pd$Study.ID)
N <- length(study.ids)

for (ix in 1:N) {
  beta_pd.incomplete <- beta_pd[beta_pd$Study.ID != study.ids[ix], ]
  beta_combat.incomplete <- beta_combat[,paste0('X', beta_pd.incomplete$Sentrix_ID, '_', beta_pd.incomplete$Sentrix_Position)]
  
  # Covariates
  m <- match(beta_pd.incomplete$Sample_Name, epidish_cellprop$X)
  m.cov <- epidish_cellprop[m, 2:8]
  m.cov$timepoint <- factor(beta_pd.incomplete$Timepoint)
  individual <- beta_pd.incomplete$Study.ID
  gender <- factor(beta_pd.incomplete$Sex)
  
  design <- model.matrix(~0+timepoint+gender+m.cov$B+m.cov$NK+m.cov$CD4T+m.cov$CD8T+m.cov$Mono+m.cov$Neutro, data = m.cov)
  colnames(design) <- c("t1", "t2", "gender", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro")
  
  corfit <- duplicateCorrelation(beta_combat.incomplete,design,block=individual)
  corfit$consensus
  
  fit <- lmFit(beta_combat.incomplete,design,block=individual,correlation=corfit$consensus)
  
  contMatrix <- makeContrasts(t2-t1,levels=design)
  fit2 <- contrasts.fit(fit, contMatrix)
  fit2 <- eBayes(fit2)
  summary(decideTests(fit2))
  DMCs <- topTable(fit2, n = Inf, coef=1)
  
  save(DMCs, file = paste0(main_dir, output_dir, "DMCs_longitudinal_", ix ,".RData"))
}

#####################################################################
# Robustness                                                        #
#####################################################################

load("/COVUM/DMCs/DMCs_longitudinal.RData")
A <- rownames(subset(DMCs, DMCs$P.Value < 0.01))
C <- rownames(subset(DMCs, DMCs$P.Value >= 0.01))

sensitivity <- c()
specificity <- c()
precision <- c()
npv <- c()
accuracy <- c()

for (ix in 1:17){
  load(paste0(main_dir, output_dir, "DMCs_longitudinal_", ix ,".RData"))
  B <- rownames(subset(DMCs, DMCs$P.Value < 0.01))
  
  TP <- length(intersect(A, B))
  TN <- length(C) - length(intersect(C, B))
  FP <- length(setdiff(B, A))
  FN <- length(setdiff(A, B))
  
  sensitivity[ix] <- TP/(TP+FN)
  specificity[ix] <- TN/(TN+FP)
  precision[ix] <- TP/(TP+FP)
  npv[ix] <- TN/(TN+FN)
  accuracy[ix] <- (TP+TN)/(TP+FP+TN+FN)
}

median(sensitivity)
median(specificity)
median(precision)
median(npv)
median(accuracy)