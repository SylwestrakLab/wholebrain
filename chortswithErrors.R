 
device = TRUE
region.lab = "Input region:" 
                               
xlab= 'Pct. of Cells'
countThreshold = 3
sorted = FALSE 
cellColors=Colors
genotypesToRun=g

  counts<-cohort
#  if (device) {
 #   #quartz(width = 7.036585, height = 0.2099039 * nrow(counts))
  #  quartz(width = 17.036585, height = 20.2099039)
  #}
  #maxCounts<-apply(counts, 1, FUN=max)
  #idxKeepRow<- which(maxCounts>3)
  #counts<-counts[idxKeepRow,]
  
  
  quartz(width = 7.036585, height = 10.2099039)
  layout(matrix(c(1, 1, 1, 2, 2, 2, 2), nrow = 1))
  par(mar = c(4, 2, 4, 2))
  par(mar=rep(1,4))
  plot(rep(2.5, nrow(counts)), nrow(counts):1, col = 0, 
       axes = F, ylim = c(0.5, nrow(counts) + 0.5), ylab = "", 
       xlab = "", xlim = c(1, 5))
  mtext(region.lab, 3, cex = 0.9)
  
  
  
  
  #Draw Labels of Brain Areas
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
  
  #Make a legend
  par(xpd = TRUE)
  
  legend(3, 1.5, genotypeNames[genotypesToRun], pch = c(21), pt.bg = cellColors[genotypesToRun], 
         title = "Genotype:", bg = "white", 
         horiz = TRUE, cex = .75, xjust = 0.5)
  par(xpd = FALSE)
  par(mar = c(4, 4, 4, 6))
  #zeros <- min(is.finite(counts)) - 1
  #counts[!is.finite(counts)] <- zeros
  
  
  #x.range <- ceiling(range(counts[,0:(length(mice)-1)]))
  x.range<-c(0, 25)
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
  axis(1, at = seq(x.range[1], x.range[2],5), las = 1, labels = seq(x.range[1], x.range[2],5))
  axis(3, at = seq(x.range[1], x.range[2],5), las = 1, labels = seq(x.range[1], x.range[2],5))
  #log.range <- unlist(lapply(1:(length(log.range) - 1), 
  #                          function(x) {
  #                             seq(log.range[x], log.range[x + 1], by = log.range[x])
  #                           }))
  #axis(1, at = log10(log.range), labels = FALSE)
  #axis(3, at = log10(log.range), labels = FALSE)
  axis(2, at = nrow(counts):1, labels = row.names(counts), 
       las = 1)
  axis(4, at = nrow(counts):1, labels = row.names(counts), 
       las = 1)
  abline(h = 1:nrow(counts), lty = 2, col = "gray")
  abline(v = seq(x.range[1], x.range[2],5), col = "lightblue")
  for (cc in genotypesToRun) {
    idx<-counts[, cc]>-1
    points(counts[which(idx), cc], nrow(counts)+1-which(idx), 
           pch = 21, bg = cellColors[cc], cex = 1.2)
  }
  box()
  par(xpd = TRUE)
  polygon(c(x.range[1] + 0.25, x.range[1] + 0.5, x.range[1] + 
              0.5, x.range[1] + 0.25), c(-2, -2, nrow(counts) + 
                                           3, nrow(counts) + 3), col = "white", border = "white")
  par(xpd = FALSE)
  abline(v = c(x.range[1] + 0.25, x.range[1] + 0.5))
  mtext(xlab, 3, 2.2, cex = 0.8)
  mtext(xlab, 1, 2.2, cex = 0.8)
}