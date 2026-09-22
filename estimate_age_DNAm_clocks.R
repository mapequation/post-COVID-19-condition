#####################################################################
# Age Prediction from Fitted Linear Models                          #
#####################################################################

# Re-trained coefficients (WB)

# Author: David Martínez-Enguita (2024)

pack_R <- c("dplyr", "glmnet", "tictoc", "MLmetrics")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

set.seed(777)

# Set directories
main_dir <- "/home/davma27/Documents/COVUM_Project/"
probe_dir <- "Data/DNAm_clock_cpgs/"
model_dir <- "Models/Regression_benchmark/"
output_dir <- "/home/davma27/Documents/COVUM_Project/Results/"

# Set the maximum and minimum age for minmax scaling
maxage <- 114
minage <- 0

# Auxiliary functions
Fage <- function(x, adult.age = 20) {
  #' function to go from age in years to Horvath’s transformed age
  #' INPUT:
  #'      x - the age to be transformed
  #'      adult.age - the age after which age should be log transformed
  #' OUTPUT:
  #'      y - the transformed age
  x <- (x + 1)/(1 + adult.age)
  y <- ifelse(x <= 1, log(x), x - 1)
  return(y)
}

invFage <- function(x, adult.age = 20) {
  #' Function that inverts the transformation in Fage function.
  #' INPUT:
  #'      x - the age to be transformed
  #'      adult.age - the age after which age should be log transformed
  #' OUTPUT:
  #'      y - the transformed age
  ifelse(x < 0, (1 + adult.age)*exp(x) - 1, (1 + adult.age)*x + adult.age)
}

# Load probe coefficients for DNAm age clocks
horvath_probes <- read.csv(paste0(main_dir, probe_dir, "/horvath_probe_coefs.csv"))
horvath_skin_probes <- read.csv(paste0(main_dir, probe_dir, "/horvath_skinblood_probe_coefs.csv"))
hannum_probes <- read.csv(paste0(main_dir, probe_dir, "/hannum_probe_coefs.csv"))
phenoage_probes <- read.csv(paste0(main_dir, probe_dir, "/phenoage_probe_coefs.csv"))
zhang_probes <- read.csv(paste0(main_dir, probe_dir, "/zhang_elnet_probe_coefs.csv"))

# Load fitted models
fitted_horvath <- readRDS(paste0(main_dir, model_dir, "/fitted_Horvath.RDS"))
fitted_hannum <- readRDS(paste0(main_dir, model_dir, "/fitted_Hannum.RDS"))
fitted_horvath_skin <- readRDS(paste0(main_dir, model_dir, "/fitted_Horvath_skin.RDS"))
fitted_phenoage <- readRDS(paste0(main_dir, model_dir, "/fitted_PhenoAge.RDS"))
fitted_zhang <- readRDS(paste0(main_dir, model_dir, "/fitted_Zhang.RDS"))

# Load data
sel_test <- beta_combat2
sel_pd <- beta_pd

# Predict age from fitted models
pred_out_df <- as.data.frame(matrix(data = NA, nrow = ncol(sel_test), ncol = 14))
colnames(pred_out_df) <- c("Horvath", "Hannum", "Horvath_skin", "PhenoAge", "Zhang",
                           "NCAE_Age", "Median_prediction", 
                           "True_Age", "Sample_Name", "Sample_Group")
pred_out_df[, (ncol(pred_out_df)-2)] <- sel_pd$Age
pred_out_df[, (ncol(pred_out_df)-1)] <- sel_pd$Sample_Name
pred_out_df[, ncol(pred_out_df)] <- sel_pd$Group_tag

fitted_list <- list(fitted_horvath, fitted_hannum, fitted_horvath_skin,
                    fitted_phenoage, fitted_zhang)
reg_probes <- list(horvath_probes, hannum_probes, horvath_skin_probes,
                   phenoage_probes, zhang_probes)

for (i in 1:length(fitted_list)) {
  fitted_model <- fitted_list[[i]]
  sel_probes <- reg_probes[[i]]
  
  # Identify intercept
  if (sel_probes[1, 1] == "intercept") {
    intercept <- sel_probes[1, ]
    sel_probes <- sel_probes[-1, ]
  } else {
    intercept <- FALSE
  }
  
  # Subset by age clock probe set
  sel_test_probes <- sel_test[match(sel_probes$probe, rownames(sel_test)), ]
  sel_test_probes <- sel_test_probes[complete.cases(sel_test_probes), ]
  sel_test_probes <- sel_test_probes[rownames(sel_test_probes) %in% fitted_model$beta@Dimnames[[1]], ]
  sel_test_probes <- sel_test_probes[match(rownames(sel_test_probes), fitted_model$beta@Dimnames[[1]]), ]

  # Perform predictions
  for (j in 1:ncol(sel_test)) {
    x <- as.matrix(sel_test_probes[, j])
    pred_age <- as.data.frame(predict(fitted_model, newx = t(x)))
    
    if (nrow(sel_test_probes) == 328 | nrow(sel_test_probes) == 380 ) {
       pred_age <- sapply(pred_age, invFage)
       }
     pred_out_df[j, i] <- pred_age
  }
}

# Store output summary
write.csv(pred_out_df, file = paste0(output_dir, "pred_age_", sel_dataset, ".csv"),
          row.names = FALSE)
          
