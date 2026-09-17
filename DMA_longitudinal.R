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

# Covariates
m <- match(beta_pd$Sample_Name, epidish_cellprop$X)
m.cov <- epidish_cellprop[m, 2:8]
m.cov$timepoint <- factor(beta_pd.postcovid$Timepoint)
m.cov$age <- beta_pd$Age
sample.group <- factor(beta_pd$Group)
gender <- factor(beta_pd$Sex)


#####################################################################
# Fold change: T2 vs. T1                                            #
#####################################################################

design <- model.matrix(~0+timepoint+gender+m.cov$B+m.cov$NK+m.cov$CD4T+m.cov$CD8T+m.cov$Mono+m.cov$Neutro, data = m.cov)
colnames(design) <- c("t1", "t2", "gender", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

corfit <- duplicateCorrelation(beta_combat,design,block=individual)
corfit$consensus

fit <- lmFit(beta_combat,design,block=individual,correlation=corfit$consensus)

contMatrix <- makeContrasts(t2-t1,levels=design)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

DMCs.longitudinal <- topTable(fit2, n = Inf, coef=1)

save(DMCs.longitudinal, file = paste0(main_dir, output_dir, "DMCs_longitudinal.RData"))


# Regional Differential Methylation (DMRcate)

design <- model.matrix(~timepoint+gender+m.cov$B+m.cov$NK+m.cov$CD4T+m.cov$CD8T+m.cov$Mono+m.cov$Neutro, data = m.cov)
corfit <- duplicateCorrelation(beta_combat,design,block=individual)

myAnnotation <- cpg.annotate(object = beta_combat,
                             datatype = "array",
                             what = "Beta",
                             analysis.type = "differential",
                             design = design,
                             coef = 2,
                             arraytype = "EPIC",
                             block = individual,
                             correlation=corfit$consensus.correlation)

DMRs.longitudinal <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = 0.05)

save(DMRs.longitudinal, file = paste0(main_dir, output_dir, "DMRs_longitudinal.RData"))
