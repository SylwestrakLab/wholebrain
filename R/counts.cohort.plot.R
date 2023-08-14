source(file.path(getwd(),"counts.cohort.plot.R"))
counts.cohort.plot<-function (dataset, device = TRUE, region.lab = "Input region:", xlab= 'Cell count') 
{
  if(ncol(dataset)<3){
    counts<-dataset
  }else{
    genotypeLabes<- matrix(rep(genotype,1),3)
    counts <- table(as.character(dataset$acronym), dataset$right.hemisphere)
    z <- length(dataset$right.hemisphere)
    genotypeLabels<- matrix(rep(genotype,1),z)
    counts <- table(as.character(dataset$acronym), genotypeLabels)
  }
  
  maxCounts<-apply(counts, 1, FUN=max)
  idxKeepRow<- which(maxCounts>0)
  counts<-counts[idxKeepRow,]
  hemisphere.to.sort <- which.max(colSums(counts*is.finite(counts),na.rm=TRUE ))
  #counts <- counts[order(counts[, hemisphere.to.sort], decreasing = TRUE), ]
  counts <- counts[order(counts[, color], decreasing = TRUE), ]
  
  counts <- log10(counts)
  
  if (device) {
    #quartz(width = 7.036585, height = 0.2099039 * nrow(counts))
    quartz(width = 17.036585, height = 20.2099039)
  }
  quartz(width = 7.036585, height = 10.2099039)
  layout(matrix(c(1, 1, 1, 2, 2, 2, 2), nrow = 1))
  par(mar = c(4, 2, 4, 2))
  par(mar=rep(1,4))
  plot(rep(2.5, nrow(counts)), nrow(counts):1, col = 0, 
       axes = F, ylim = c(0.5, nrow(counts) + 0.5), ylab = "", 
       xlab = "", xlim = c(1, 5))
  mtext(region.lab, 3, cex = 0.9)
  for (i in 1:nrow(counts)) {
    regioncolor <- color.from.acronym(row.names(counts)[i])
    regioncolor <- adjustcolor(regioncolor, alpha.f = 0.1)
    y.lab <- (nrow(counts) + 1) - i
    polygon(c(1, 5, 5, 1), c(y.lab - 0.5, y.lab - 0.5, 
                             y.lab + 0.5, y.lab + 0.5), col = regioncolor, 
            border = FALSE)
    text(3, y.lab, name.from.acronym(row.names(counts)[i]), 
         cex = 0.9)
  }
  par(xpd = TRUE)
  legend(3, -0.75, c("Left", "Right"), pch = c(21), pt.bg = c("white", 
                                                              gray(0.2)), title = "Hemisphere:", bg = "white", 
         horiz = TRUE, cex = 1.3, xjust = 0.5)
  par(xpd = FALSE)
  par(mar = c(4, 4, 4, 6))
  zeros <- min(is.finite(counts)) - 1
  counts[!is.finite(counts)] <- zeros
  x.range <- ceiling(range(counts[is.finite(counts)]))
  plot(apply(counts, 1, max), nrow(counts):1 - 0.125, pch = 21, 
       bg = "white", ylim = c(0.5, nrow(counts) + 0.5), 
       xlim = x.range, xlab = "", axes = F, ylab = "", col = 0)
  for (i in 1:nrow(counts)) {
    regioncolor <- color.from.acronym(row.names(counts)[i])
    regioncolor <- adjustcolor(regioncolor, alpha.f = 0.05)
    y.lab <- (nrow(counts) + 1) - i
    polygon(c(x.range, rev(x.range)), c(y.lab - 0.5, 
                                        y.lab - 0.5, y.lab + 0.5, y.lab + 0.5), col = regioncolor, 
            border = FALSE)
  }
  log.range <- 10^seq(x.range[1], x.range[2])
  axis(1, at = seq(x.range[1], x.range[2]), las = 1, labels = c(0, 
                                                                log.range[-1]))
  axis(3, at = seq(x.range[1], x.range[2]), las = 1, labels = c(0, 
                                                                log.range[-1]))
  log.range <- unlist(lapply(1:(length(log.range) - 1), 
                             function(x) {
                               seq(log.range[x], log.range[x + 1], by = log.range[x])
                             }))
  axis(1, at = log10(log.range), labels = FALSE)
  axis(3, at = log10(log.range), labels = FALSE)
  axis(2, at = nrow(counts):1, labels = row.names(counts), 
       las = 1)
  axis(4, at = nrow(counts):1, labels = row.names(counts), 
       las = 1)
  abline(h = 1:nrow(counts), lty = 2, col = "gray")
  abline(v = log10(log.range), col = "lightblue")
  lapply(1:nrow(counts), function(x) {
    lines(counts[x, ], rep(nrow(counts) - x + 1, 2), 
          lwd = 2)
  })
  points(counts[, 2], nrow(counts):1, 
         pch = 21, bg = gray(0.2), cex = 1.2)
  points(counts[, 1], nrow(counts):1, 
         pch = 21, bg = "white", cex = 1.2)
  box()
  par(xpd = TRUE)
  #polygon(c(x.range[1] + 0.25, x.range[1] + 0.5, x.range[1] + 
  #            0.5, x.range[1] + 0.25), c(-2, -2, nrow(counts) + 
  #                                         3, nrow(counts) + 3), col = "white", border = "white")
  par(xpd = FALSE)
  #abline(v = c(x.range[1] + 0.25, x.range[1] + 0.5))
  mtext(xlab, 3, 2.2, cex = 0.8)
  mtext(xlab, 1, 2.2, cex = 0.8)
}
  