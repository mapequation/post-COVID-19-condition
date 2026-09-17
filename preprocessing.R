#####################################################################
# Preprocessing of methylation datasets from IDAT                   #
#####################################################################

# 1. Formatting of phenotypic metadata
# 3. Normalization with ENmix
# 3. Filtering and QC with ChAMP
# 4. Detection of batch effects using SVD
# 5. Batch effect correction with ComBat
# 6. Cell type analysis

library(minfi)
library(bnstruct)
library(sva)
library(ChAMP)
library(FlowSorted.Blood.EPIC)
library(EpiDISH)
library(readxl)
library(wateRmelon) #BMIQ
library(ENmix)

set.seed(777)

main_dir <- "COVUM/"
input_dir <- "Data/Input/"
output_dir <- "Data/Normalized/"

#####################################################################
# 1. Metadata                                                       #
#####################################################################

load(paste0(main_dir, input_dir, "beta_pd.RData"))
beta_pd <- beta_pd[(beta_pd$Timepoint == 3) | (beta_pd$Timepoint == 12), ]

#####################################################################
# 2. Normalization with ENmix                                       #
#####################################################################

load(paste0(main_dir, input_dir, "raw_rgset.RData"))

RGSet <- RGSet[, paste0(beta_pd$Sentrix_ID, '_', beta_pd$Sentrix_Position)]

mdat<-preprocessENmix(RGSet, bgParaEst="oob", dyeCorr="RELIC",nCores=4)
beta_enmix<-rcp(mdat)

# Output: beta_enmix.RData
save(beta_enmix, file = paste0(main_dir, output_dir, "beta_enmix.RData"))

#####################################################################
# 3. Filtering and QC with ChAMP                                    #
#####################################################################

colnames(beta_enmix) <- paste0("X", colnames(beta_enmix))

# Filtering
out_filt <- champ.filter(beta           = as.matrix(beta_enmix), 
                         pd             = beta_pd,
                         autoimpute     = FALSE,
                         detP           = NULL,
                         filterDetP     = FALSE,
                         ProbeCutoff    = 0,
                         SampleCutoff   = 0.1,
                         detPcut        = 0.01,
                         filterBeads    = FALSE,
                         filterNoCG     = TRUE,
                         filterSNPs     = TRUE,
                         filterMultiHit = TRUE,
                         filterXY       = TRUE,
                         fixOutlier     = TRUE, 
                         arraytype      = "EPIC")

beta_filt <- out_filt$beta

# Output: beta_enmix_filt.RData
save(beta_filt, file = paste0(main_dir, output_dir, "beta_enmix_filt.RData"))

#####################################################################
# 4. Detection of batch effects using SVD                           #
#####################################################################

# Select PD file covariates to analyze
beta_pd_svd <- beta_pd[, c("Study.ID", "Sentrix_ID", "Sentrix_Position", "Group",
                           "Sex", "Age", "BMI")]
beta_pd_svd[, ] <- lapply(beta_pd_svd[, ], as.character)

# Calculate first 10 SVD components (champ.SVD)
source(paste0(main_dir, "calculate_svd.R"))

pdf(file = paste0(main_dir, output_dir, "svd_enmix_filt_beta.pdf"), 
    width = 10, height = 8)
calculate_svd(counts = beta_filt, 
              metadata = beta_pd_svd, 
              components = 10, 
              Rplot = TRUE)
dev.off()

#####################################################################
# 5. Batch effect correction with ComBat                            #
#####################################################################

beta_filt <- beta_filt[,paste0('X', beta_pd$Sentrix_ID, '_', beta_pd$Sentrix_Position)]

batchname <- c("Sentrix_ID", "Sentrix_Position")
beta_corr <- beta_filt
for (batchvar in batchname){
  beta_corr <- ComBat(
    dat=beta_corr,
    batch=beta_pd_svd[, batchvar],
    mod = NULL,
  )
}
beta_combat <- beta_corr

# Output: beta_enmix_filt_combat.RData
save(beta_combat, file = paste0(main_dir, output_dir, "beta_enmix_filt_combat.RData"))

pdf(file = paste0(main_dir, output_dir, "svd_enmix_filt_combat_beta.pdf"), 
    width = 10, height = 8)
calculate_svd(counts = beta_combat, 
              metadata = beta_pd_svd, 
              components = 10, 
              Rplot = TRUE)
dev.off()

#####################################################################
# 6. Cell type analysis                                             #
#####################################################################

# EpiDISH : estimated cell proportions for Epi, Fib, Fat and IC, with robust partial correction
data(centEpiFibIC.m)
epidish_cellprop <- epidish(beta.m = beta_combat, 
                            ref.m = centDHSbloodDMC.m, 
                            method = "RPC")$estF

pdf(file = paste0(main_dir, output_dir, "boxplot_epidish.pdf"), 
    width = 8, height = 8)
boxplot(epidish_cellprop)
dev.off()

write.csv(epidish_cellprop, file = paste0(main_dir, output_dir, "/beta_epidish_prop.csv"),
          quote = FALSE, row.names = TRUE)
