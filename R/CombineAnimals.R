#Combine for each genotype
library(wholebrain)
setwd(cwd)
setwd("/Users/esylwestrak/Dropbox (University of Oregon)/UO-Sylwestrak Lab/Rabies_registration/wholebrain/")
genotype <- 'Tac1-LHB'
#set path for data combined by animal
combindedAnimalFolder<- file.path(getwd(),'Composite Data',genotype,'CompositeData')
mouseFiles<-list.files(combindedAnimalFolder, full.names=TRUE, pattern = "\\.Rdata$")

#Start with the first section
i = 1
firstMouse = mouseFiles[i]
if (length(firstMouse)==1) {
  load(firstMouse[1])
  Genotype_datasets<-datasets
}


#Use Rbind to add addictional sections
for(i in 1:length(mouseFiles)){
    load(mouseFiles[i])
    Genotype_datasets <-rbind(Genotype_datasets, datasets)
}

##Plot the combined data and save out a combined .Rdata file
counts.plot(Genotype_datasets)
filename <- file.path(getwd(),'Composite Data',genotype,paste(genotype, '_composite_plot.pdf',sep = "", collapse=NULL))
quartz.save(filename,type="pdf")
save(Genotype_datasets, file=file.path(getwd(), 'Composite Data',genotype,paste(genotype, '_compositeData.Rdata',sep = "", collapse=NULL)))



#Plot all areas so to compare across genotypes
counts.plot(Genotype_datasets, countThreshold=3)
filename <- file.path(getwd(),'Composite Data',genotype,paste(genotype, '_composite_plot.pdf',sep = "", collapse=NULL))
quartz.save(filename,type="pdf")

genotypeFolders<-list.dirs(file.path(getwd(),'Composite Data'), full.names = TRUE, recursive = FALSE)
genotypeFolders=genotypeFolders[c(1,2,4,5, 7,8,10)]
genotypes <-list.dirs(file.path(getwd(),'Composite Data'), full.names = FALSE, recursive = FALSE)
genotypes = genotypes[c(1,2,4,5, 7,8,10)]
i<-1
fname = list.files(genotypeFolders[i], full.names=TRUE, pattern = "\\.Rdata$")
load(fname)
cohortData <-Genotype_datasets
cohortData$genotype = genotypes[i]
cohortData <- cohortData[,!names(cohortData) %in% c("scaleFactor", "somaX", "somaY")]

i<-2
for(i in 2:length(genotypeFolders)){
  fname = list.files(genotypeFolders[i], full.names=TRUE, pattern = "\\.Rdata$")
  load(fname)
  Genotype_datasets <- Genotype_datasets[,!names(Genotype_datasets) %in% c("scaleFactor", "somaX", "somaY")]
  Genotype_datasets$genotype = genotypes[i]
  cohortData<-rbind(cohortData, Genotype_datasets)
  
}

regions <- unique(cohortData$acronym)


d = matrix(0, nrow = length(regions), ncol = length(genotypes))
for(i in 1:length(regions)){
  area<-subset(cohortData, acronym == regions[i])
  for(c in 1:length(genotypes)){
  d[i, c]<-nrow(subset(area, genotype == genotypes[c]))
  #d[i, c+1]<- area$color[1]
  }
}

colnames(d) <- genotypes

rownames(d) <- regions
cohort <- as.table(d)
row.names.remove <- c("LH", "MH","MPT","TH", "PVT","fiber tracts","TH","PF","PRC","V3","IMD","grey")
cohort<-cohort[!(row.names(cohort) %in% row.names.remove), ]

total <-colSums(cohort[,c(1,2,3,4,5,6,7)])

for(i in 1:length(genotypes)-2){
  cohort[,i]<-cohort[,i]/total[i]*100
}
cohort[,6]<-cohort[,6]/total[5]*100
cohort[,7]<-cohort[,7]/total[7]*100
#cohort <- cbind(cohort, Total = rowSums(cohort))
#cohort <- cohort[order(cohort[, 3], decreasing = TRUE), ]

#maxCounts<-apply(cohort, 1, FUN=max)

for(i in 1:nrow(cohort)){
  maxCounts[i]<-max(cohort[i,],na.rm = TRUE)
}

idxKeepRow<- which(maxCounts>3)
cohort<-cohort[idxKeepRow,]



Colors<-c("orange","red","yellow","purple","blue","white","green")
genotypeNames<-c("Calb1", "ChAT","Gad2","Oprm","Tac1-LHb", "Tac1-noG","Th")

#Colors<-c("orange","blue","black")
#genotypeNames<-c("Gad2", "Tac1LHb","Tac1-noG")

#Colors<-c("blue","white")
#genotypeNames<-c("Tac1-LHb","Tac1-noG")


#d <- d[order(d[, length(genotypes)+1], increaseing = TRUE), ]
#d <- d[,-c(length(genotypes)+1)]
#cohort<-d

g<-c(1,2,3,4,5,6,7)
counts.compare.plot<-function (cohort, device = TRUE, region.lab = "Input region:", 
                               xlab= 'Pct. of Cells', countThreshold = 3, sorted = FALSE, cellColors=Colors, genotypesToRun=g)
{
  counts<-cohort
  if (device) {
    #quartz(width = 7.036585, height = 0.2099039 * nrow(counts))
    quartz(width = 17.036585, height = 20.2099039)
  }
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
counts.compare.plot(cohort, device = TRUE, region.lab = "Input region:", xlab= 'Pct. of Cells', countThreshold = 1, sorted = TRUE) 
filename <- file.path(getwd(),'Composite Data',paste('CompareCohorts.pdf',sep = "", collapse=NULL))
quartz.save(filename,type="pdf")

