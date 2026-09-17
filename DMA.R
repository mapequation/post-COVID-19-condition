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
beta_pd <- beta_pd[beta_pd$Timepoint == 3, ]
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
# Fold change: PCC+ vs. PCC-                                        #
#####################################################################

design <- model.matrix(~0+sample.group+gender+m.cov$age+m.cov$B+m.cov$NK+m.cov$CD4T+m.cov$CD8T+m.cov$Mono+m.cov$Neutro, data = m.cov)
colnames(design) <- c("control", "pcc", "gender", "age", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro")

fit<-lmFit(betaVal,design)
contrast.matrix <- makeContrasts(pcc-control,
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
summary(decideTests(fit2))

DMCs <- topTable(fit2, n = Inf, coef=1)

save(DMCs, file = paste0(main_dir, output_dir, "DMCs_pcc_vs_control.RData"))


# Regional Differential Methylation (DMRcate)

myAnnotation <- cpg.annotate(object = beta_combat,
                             datatype = "array",
                             what = "Beta",
                             analysis.type = "differential",
                             design = design,
                             contrasts = TRUE,
                             cont.matrix = contrast.matrix,
                             coef = "pcc - control",
                             arraytype = "EPIC")

DMRs <- dmrcate(myAnnotation, lambda=1000, C=2, pcutoff = 0.05)

save(DMRs, file = paste0(main_dir, output_dir, "DMRs_pcc_vs_control.RData"))
