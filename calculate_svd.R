calculate_svd <- function (counts, metadata, components = 10, Rplot = TRUE) {
  GenPlot <- function(thdens.o, estdens.o, evalues.v) {
    minx <- min(min(thdens.o$lambda), min(evalues.v))
    maxx <- max(max(thdens.o$lambda), max(evalues.v))
    miny <- min(min(thdens.o$dens), min(estdens.o$y))
    maxy <- max(max(thdens.o$dens), max(estdens.o$y))
  }
  
  EstDimRMTv2 <- function(data.m) {
    M <- data.m
    for (c in 1:ncol(M)) M[, c] <- (data.m[, c] - mean(data.m[, c]))/sqrt(var(data.m[, c]))
    sigma2 <- var(as.vector(M))
    Q <- nrow(data.m)/ncol(data.m)
    thdens.o <- thdens(Q, sigma2, ncol(data.m))
    C <- 1/nrow(M) * t(M) %*% M
    eigen.o <- eigen(C, symmetric = TRUE)
    estdens.o <- density(eigen.o$values, from = min(eigen.o$values), 
                         to = max(eigen.o$values), cut = 0)
    GenPlot(thdens.o, estdens.o, eigen.o$values)
    intdim <- length(which(eigen.o$values > thdens.o$max))
    return(list(cor = C, dim = intdim, estdens = estdens.o, 
                thdens = thdens.o))
  }
  
  thdens <- function(Q, sigma2, ns) {
    lambdaMAX <- sigma2 * (1 + 1/Q + 2 * sqrt(1/Q))
    lambdaMIN <- sigma2 * (1 + 1/Q - 2 * sqrt(1/Q))
    delta <- lambdaMAX - lambdaMIN
    roundN <- 3
    step <- round(delta/ns, roundN)
    while (step == 0) {
      roundN <- roundN + 1
      step <- round(delta/ns, roundN)
    }
    lambda.v <- seq(lambdaMIN, lambdaMAX, by = step)
    dens.v <- vector()
    ii <- 1
    for (i in lambda.v) {
      dens.v[ii] <- (Q/(2 * pi * sigma2)) * sqrt((lambdaMAX - i) * (i - lambdaMIN))/i
      ii <- ii + 1
    }
    return(list(min = lambdaMIN, max = lambdaMAX, step = step, 
                lambda = lambda.v, dens = dens.v))
  }
  
  drawheatmap <- function(svdPV.heat) {
    svdPV.heat <- svdPV.heat[, !colnames(svdPV.heat) == "Explained_variance"]
    myPalette <- c("darkred", "red", "orange", "pink", "white")
    breaks.v <- c(-10000, -10, -5, -2, log10(0.05), 0)
    image(x = 1:nrow(svdPV.heat), 
          y = 1:ncol(svdPV.heat), 
          z = log10(svdPV.heat), 
          col = myPalette, breaks = breaks.v, xlab = "", ylab = "", 
          axes = FALSE, main = "Singular Value Decomposition Analysis (SVD)")
    axis(1, at = 1:nrow(svdPV.heat), labels = paste0("C", 1:nrow(svdPV.heat)), las = 2)
    suppressWarnings(axis(2, at = 1:ncol(svdPV.heat), labels = colnames(svdPV.heat), las = 2))
    legend(x = -(topPCA/2.5), y = ncol(metadata), legend = c(expression("p < 1 x" ~ 10^{-10}), 
                                                             expression("p < 1 x" ~ 10^{-5}), 
                                                             "p < 0.01", "p < 0.05", "p > 0.05"), 
           fill = c("darkred", "red", "orange", "pink", "white"), par("usr")[2], 
           par("usr")[4], xpd = NA, border = "black", seg.len = 0)
  }
  
  splot <- function(x = svd.o, y = rmt.o) {
    scp <- list(u = x$u[1:nrow(x$u), 1:y$dim], v = x$v[1:y$dim, 1:y$dim], d = x$d[1:y$dim])
  }
  
  if (class(counts)[1] == "data.frame") {
    counts <- as.matrix(counts) 
  }
  if (class(metadata)[1] == "matrix") {
    metadata <- as.data.frame(metadata)
  }
  PhenoTypes.lv <- metadata
  
  if (!is.null(rownames(metadata))) {
    rownames(PhenoTypes.lv) <- rownames(metadata)
  }
  tmp.m <- counts - rowMeans(counts)
  rmt.o <- EstDimRMTv2(tmp.m)
  svd.o <- svd(tmp.m)
  if (rmt.o$dim > 20) {
    topPCA <- components
  } else topPCA <- components
  svdPV.m <- matrix(nrow = topPCA, ncol = ncol(PhenoTypes.lv))
  colnames(svdPV.m) <- colnames(PhenoTypes.lv)
  rownames(svdPV.m) <- paste0("C", 1:nrow(svdPV.m))
  for (c in 1:topPCA) for (f in 1:ncol(PhenoTypes.lv)) if (class(PhenoTypes.lv[, f]) != "numeric") 
    svdPV.m[c, f] <- kruskal.test(svd.o$v[, c] ~ as.factor(PhenoTypes.lv[[f]]))$p.value
  else svdPV.m[c, f] <- summary(lm(svd.o$v[, c] ~ PhenoTypes.lv[[f]]))$coeff[2, 4]
  
  max_comp <- nrow(metadata)
  scp <- list(u = svd.o$u[1:nrow(svd.o$u), 1:max_comp], v = svd.o$v[1:max_comp, 1:max_comp], d = svd.o$d[1:max_comp])
  pve <- scp$d^2/sum(scp$d^2) * 100
  svdPV.m <- cbind(svdPV.m, pve[1:components])
  colnames(svdPV.m)[ncol(svdPV.m)] <- "Explained_variance"
  
  if (Rplot) {
    par(mar = c(5, 15, 2, 1))
    splot(svd.o, rmt.o)
    drawheatmap(svdPV.heat = svdPV.m)
  }
  return(svdPV.m)
}
