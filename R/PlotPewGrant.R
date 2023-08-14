counts.compare.plot<-function (cohort, device = TRUE, region.lab = "Input region:", xlab= 'Pct. of Cells', countThreshold = 10, sorted = TRUE) 
{
  counts<-cohort
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
  legend(3, 1.5, c("Calb1", "Tac1","Tac1-PVT","Th"), pch = c(21), pt.bg = c("orange","blue","white","green"), 
         title = "Genotype:", bg = "white", 
         horiz = TRUE, cex = .75, xjust = 0.5)
  par(xpd = FALSE)
  par(mar = c(4, 4, 4, 6))
  zeros <- min(is.finite(counts)) - 1
  counts[!is.finite(counts)] <- zeros
  x.range <- ceiling(range(counts[,1:4]))
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
  log.range <- unlist(lapply(1:(length(log.range) - 1), 
                             function(x) {
                               seq(log.range[x], log.range[x + 1], by = log.range[x])
                             }))
  #axis(1, at = log10(log.range), labels = FALSE)
  #axis(3, at = log10(log.range), labels = FALSE)
  axis(2, at = nrow(counts):1, labels = row.names(counts), 
       las = 1)
  axis(4, at = nrow(counts):1, labels = row.names(counts), 
       las = 1)
  abline(h = 1:nrow(counts), lty = 2, col = "gray")
  abline(v = seq(x.range[1], x.range[2],5), col = "lightblue")
  points(counts[, 1], nrow(counts):1, 
         pch = 21, bg = "orange", cex = 1.2)
  points(counts[, 2], nrow(counts):1, 
         pch = 21, bg = "blue", cex = 1.2)
  points(counts[, 3], nrow(counts):1, 
         pch = 21, bg = 'white', cex = 1.2)
  points(counts[, 4], nrow(counts):1, 
         pch = 21, bg = "green", cex = 1.2)
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

regions <- unique(Genotype_datasets$acronym)


d = matrix(0, nrow = length(regions), ncol = length(genotypes))
for(i in 1:length(regions)){
  area<-subset(Genotype_datasets, acronym == regions[i])
  for(c in 1:length(genotypes)){
    d[i, c]<-nrow(subset(area, genotype == genotypes[c]))
  }
}

colnames(d) <- genotypes
rownames(d) <- regions
cohort <- as.table(d)
row.names.remove <- c("LH", "MH", "PVT","fiber tracts")
cohort<-cohort[!(row.names(cohort) %in% row.names.remove), ]
cohort <- cohort[order(cohort[, 3:4], decreasing = TRUE), ]

cohort <-cohort[,4]
maxCounts<-apply(cohort), 1, FUN=max)
idxKeepRow<- which(cohort[,0]>5)
cohort2<-cohort[idxKeepRow,]


counts.compare.plot(cohort, device = TRUE, region.lab = "Input region:", xlab= 'Pct. of Cells', countThreshold = 10, sorted = TRUE) 
filename <- file.path('/Users/emily/Dropbox/Grants and Fellowships/Pew',paste('CompareCohorts.pdf',sep = "", collapse=NULL))
quartz.save(filename,type="pdf")



