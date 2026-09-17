library(ggplot2)
library(ggrastr)
library(bacon)

qq_inflation_gg <- function(p, title="EWAS QQ") {
  p <- p[is.finite(p) & !is.na(p) & p > 0 & p <= 1]
  n <- length(p)
  exp <- -log10(((1:n) - 0.5) / n)
  obs <- -log10(sort(p))
  df <- data.frame(exp, obs)
  
  chisq <- qchisq(1 - p, df = 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1)

  ggplot(df, aes(exp, obs)) +
    geom_point_rast(size = 0.8, alpha = 0.8, color = "#bf812d") +  # brown markers
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(
      x = expression("Expected " * -log[10](p)),
      y = expression("Observed " * -log[10](p)),
      title = sprintf("%s (λ = %.2f)", title, lambda)
    ) +
    theme_minimal(base_size = 12) +           # minimal white theme
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size=12)
    )  
}


# Input
main_dir <- "COVUM/"
input_dir <- "DMCs/"
load(paste0(main_dir, input_dir, "DMCs.RData"))

# QQ-plot
qq_inflation_gg(DMCs$P.Value)

# Run bacon on t-statistics from limma
bc <- bacon(teststatistics = DMCs$t)

# Get corrected p-values
DMCs$P_bacon <- 2 * pnorm(-abs(tstat(bc)))
save(DMCs, file="COVUM/DMCs/DMCs_bacon.RData")

