dot.plot<-function (dataset, device = TRUE, region.lab = "Input region:", xlab= 'Cell count') 
{
  if(ncol(dataset)<3){
    counts<-dataset
  }else{
    counts <- table(as.character(dataset$acronym), dataset$right.hemisphere)
  }
  hemisphere.to.sort <- which.max(colSums(counts*is.finite(counts),na.rm=TRUE ))
  counts <- counts[order(counts[, hemisphere.to.sort], decreasing = TRUE), ]
  counts <- log10(counts)
  
  if (device) {
    #quartz(width = 7.036585, height = 0.2099039 * nrow(counts))
    quartz(width = 7.036585, height = 10.2099039)
  }
  layout(matrix(c(1, 1, 1, 2, 2, 2, 2), nrow = 1))
  par(mar = c(4, 0, 4, 0))
  #par(mar=rep(2,4))
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
  polygon(c(x.range[1] + 0.25, x.range[1] + 0.5, x.range[1] + 
              0.5, x.range[1] + 0.25), c(-2, -2, nrow(counts) + 
                                           3, nrow(counts) + 3), col = "white", border = "white")
  par(xpd = FALSE)
  abline(v = c(x.range[1] + 0.25, x.range[1] + 0.5))
  mtext(xlab, 3, 2.2, cex = 0.8)
  mtext(xlab, 1, 2.2, cex = 0.8)
  
}


plot.outlines<-function(regi, plot=FALSE){
  scale.factor<-regi$transformationgrid$width/dim(regi$transformationgrid$mx)[2]
  if(plot)
    plot(c(regi$atlas$outlines[[1]]$xrT*scale.factor, regi$atlas$outlines[[1]]$xlT*scale.factor), c(regi$atlas$outlines[[1]]$yrT*scale.factor, regi$atlas$outlines[[1]]$ylT*scale.factor), type='l', lwd=2, xlim=c(regi$transformationgrid$width, 0), ylim=c(regi$transformationgrid$height, 0), axes=FALSE, asp=1, xlab='', ylab='')
  
  for(i in 1:regi$atlas$numRegions){
    if(regi$atlas$col[i]=='#cccccc'){
      bg<-'#cccccc'
      lwd<-2
    }else{
      bg<-gray(0.95)
      lwd<-1
    }
    polygon(regi$atlas$outlines[[i]]$xlT*scale.factor, regi$atlas$outlines[[i]]$ylT*scale.factor, col=bg, lwd=lwd)
    polygon(regi$atlas$outlines[[i]]$xrT*scale.factor, regi$atlas$outlines[[i]]$yrT*scale.factor, col=bg, lwd=lwd)
  }
}


#' Cortical topview plot
#'
#' Plots a schematic brain from the top showing only parts that are easily accessible with cranial windows without removing parts of the temporalis muscle.
#' @param dataset a dataset frame obtained by get.cell.ids() or by inspect.registration().
#' @param dv.cut numeric, value dorso-ventral plane where everything above dv.cut will be plotted. Default is -2.8 mm.
#' @param show.sections boolean, if lines should be drawn for each unique AP section. Default is TRUE.
#' @param labels if labels should be plotted with acronyms for each cortical region.
#' @param col color code or name.
#' @param pch plotting ‘character’, i.e., symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols.
#' @param pt.bg background (fill) color for the open plot symbols given by pch = 21:25.
#' @examples
#' cortical.topview(dataset, right.hemisphere = FALSE, labels = FALSE)
cortical.topview<-function(dataset, dv.cut = -2.8, show.sections = TRUE, labels = TRUE,  col = NULL, pch = 16, cex = 0.5, pt.bg = NULL){
  
  brain.color<-rgb(0.25, 0.25, 0.25, 0.1)
  
  par(pty='s')
  plot(0,0, col=0, ylab='Anterior-Posterior (mm)', xlab='Medial-Lateral (mm)', xlim=c(-5,5), axes=F, ylim=c(-5,3.7))
  
  
  abline(v=c(-5:5), col='lightblue')
  abline(h=c(-5:5), col='lightblue')
  
  
  polygon(topview$outline[,1]/1.95, topview$outline[,2]/1.95, lwd=2, col= brain.color )
  polygon(-(topview$outline[,1]/1.95), topview$outline[,2]/1.95, lwd=2, col= brain.color)
  
  
  axis(1, at=seq(-4,4, by=0.1), labels=F, col='gray', tck=-0.018)
  axis(2, at=seq(-5,5, by=0.1), labels=F, col='gray', tck=-0.018)
  axis(1, at=-4:4, labels=F)
  axis(2, at=c(-5:5), labels=F)
  
  axis(1, at=seq(-4,4, by=2))
  axis(2, at=seq(-4,4, by=2))
  
  abline(h=0)
  points(0,0,pch=23, bg='white', cex=1.3, lwd=2)  
  axis(4, at=0, labels='Bregma')
  abline(h=-3.8)
  points(0,-3.8,pch=23, bg='lightblue', cex=1.3, lwd=2)   
  axis(4, at=-3.8, labels='lambda')
  
  for(i in unique(topview$regions[,3])){
    polygon(as.numeric(topview$regions[which(topview$regions[,3]==i),1])/1.92, as.numeric(topview$regions[which(topview$regions[,3]==i),2])/1.98, lty=3, border='black')
    polygon(-as.numeric(topview$regions[which(topview$regions[,3]==i),1])/1.92, as.numeric(topview$regions[which(topview$regions[,3]==i),2])/1.98, border='black', lty=3)
  }
  
  isocortex <- get.acronym.child(get.acronym.child(get.acronym.child("Isocortex")))
  isocortex <- append(isocortex, get.acronym.child(get.acronym.child("Isocortex")))
  isocortex <- na.omit(isocortex)
  isocortex <- dataset$acronym %in% isocortex
  isocortex <- (isocortex & (dataset$DV > dv.cut) )
  
  if(!is.null(col)){
    color<-col
  }else{
    color<-dataset$color[isocortex]
  }
  
  if(show.sections)
    abline(h = unique(dataset$AP), col="gray")
  
  points(dataset$ML[isocortex], dataset$AP[isocortex], pch = pch, cex = cex, col = color, pt.bg = pt.bg)  
  
}

#' Cortical sideview plot
#'
#' Plots a schematic brain from the side with cortical regions outlined.
#' @param dataset a dataset frame obtained by get.cell.ids() or by inspect.registration().
#' @param right.hemisphere boolean, if true then right side of the hemisphere will be plotted, if false left hemisphere will be plotted.
#' @param ml.cut numeric, value medio-lateral plane where everything more medial than ml.cut will be excluded. Default is 2.0 mm. 
#' @param show.sections boolean, if lines should be drawn for each unique AP section. Default is TRUE.
#' @param labels if labels should be plotted with acronyms for each cortical region.
#' @param col color code or name.
#' @param pch plotting ‘character’, i.e., symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols.
#' @param pt.bg background (fill) color for the open plot symbols given by pch = 21:25.
#' @examples
#' cortical.sideview(dataset, right.hemisphere = FALSE, labels = FALSE)
cortical.sideview<-function(dataset, right.hemisphere = TRUE, ml.cut = 2.0, show.sections = TRUE, labels = TRUE,  col = NULL, pch = 16, pt.bg = NULL, cex = 1){
  hemisphere<-c('Left hemisphere', 'Right hemisphere')
  xlim<-rbind(c(5,-7), c(-7, 5))
  plot(sideview[,1:2], xlim=xlim[right.hemisphere + 1,], asp=1, ylim=c(-7,0), col=0, axes=F, ylab='Dorso-ventral [mm]', xlab='Anterior-posterior [mm]')
  mtext(hemisphere[right.hemisphere + 1], 3, 1, font=2)
  graylight<-0.9
  lapply(5:-8, function(x)abline(v=x, col='lightblue'))
  lapply(1:-8, function(x)abline(h=x, col='lightblue'))
  
  polygon(sideview[which(sideview$acronym=='brainstem'),1:2], col=gray(graylight))
  polygon(sideview[which(sideview$acronym=='cerebellum'),1:2], col=gray(graylight))
  polygon(sideview[which(sideview$acronym=='telencephalon'),1:2], col=gray(graylight))
  
  cortical<-unique(sideview$acronym[which(nchar(sideview$acronym)<=5)] )
  
  lapply(cortical, function(x){polygon(sideview[which(sideview$acronym==x),1:2], lty=3, col=gray(graylight+0.05))} )
  
  
  isocortex <- get.acronym.child(get.acronym.child(get.acronym.child("Isocortex")))
  isocortex <- append(isocortex, get.acronym.child(get.acronym.child("Isocortex")))
  isocortex <- na.omit(isocortex)
  isocortex <- dataset$acronym %in% isocortex
  isocortex <- (isocortex & (abs(dataset$ML) > ml.cut) )
  
  if(!is.null(col)){
    color<-col
  }else{
    color<-dataset$color[isocortex & dataset$right.hemisphere == right.hemisphere]
  }
  
  if(show.sections)
    abline(v = unique(dataset$AP), col="gray")
  
  points(dataset$AP[isocortex & dataset$right.hemisphere == right.hemisphere], dataset$DV[isocortex & dataset$right.hemisphere == right.hemisphere], pch = pch, col = color, pt.bg = pt.bg, cex = cex)
  
  if(labels)
    lapply(cortical, function(x){coord<-apply(sideview[which(sideview$acronym==x),1:2], 2, median);text(coord[1], coord[2], x, cex=0.85, col=rgb(0.3,0,0.3))} )
  
  axis(2, at=c(0:-7), las=1)
  axis(1, at=c(5:-7))
}


#' Cortical plot
#'
#' Plots schematic brain(s) of cells in cortex from different views.
#' @param dataset a dataset frame obtained by get.cell.ids() or by inspect.registration().
#' @param type character vector, what type of plots should be made. Default is c('left', 'right', 'top').
#' @param show.sections boolean, if lines should be drawn for each unique AP section. Default is TRUE.
#' @param ml.cut numeric, value medio-lateral plane where everything more medial than ml.cut will be excluded. Default is 2.0 mm. 
#' @param dv.cut numeric, value dorso-ventral plane where everything above dv.cut will be plotted. Default is -2.8 mm.
#' @param labels if labels should be plotted with acronyms for each cortical region.
#' @param col color code or name.
#' @param pch plotting ‘character’, i.e., symbol to use. This can either be a single character or an integer code for one of a set of graphics symbols.
#' @param pt.bg background (fill) color for the open plot symbols given by pch = 21:25.
#' @examples
#' cortical.plot(dataset, type = c('left', 'right', 'top') )
cortical.plot<-function(dataset, type = c('left', 'right', 'top'), show.sections = TRUE, labels = TRUE, ml.cut = 2.0, dv.cut = -2.8, col = NULL, cex = 0.5, pch = 16, pt.bg = NULL){
  if(all(c('left', 'right')%in%type)  & (!'top'%in%type)){
    quartz(width=8.330275*2, height= 5.201835)
    par(mar=c(4,4,2,1), mfrow=c(1,2))
    cortical.sideview(dataset, right.hemisphere = FALSE, show.sections = show.sections, ml.cut = ml.cut, labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
    cortical.sideview(dataset, right.hemisphere = TRUE, show.sections = show.sections,ml.cut = ml.cut, labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
  }
  
  #plot only both hemispheres
  if(all(c('left', 'right', 'top')%in%type)){
    quartz(width=14.691892, height= 4.864865)
    layout(matrix(c(1, 1, 1, 1, 2, 2, 2,2,2,3, 3, 3, 3,3), 1, 14, byrow = TRUE))
    par(mar=c(4,4,4,2), yaxs="r")
    cortical.topview(dot, show.sections = FALSE)
    par(mar=c(4,2,3,0), yaxs="i" )
    cortical.sideview(dataset, right.hemisphere = FALSE, show.sections = FALSE)
    par(mar=c(4,2,3,4), yaxs="i" )
    cortical.sideview(dataset, right.hemisphere = TRUE, show.sections = FALSE)
  }
  #plot  top
  if(all(type == 'top')){
    cortical.topview(dataset, show.sections = show.sections, dv.cut = dv.cut, labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
  }
  
  #plot only left or right
  if(xor('left'%in%type, 'right'%in%type)){
    cortical.sideview(dataset, right.hemisphere = (type == 'right'), show.sections = show.sections, ml.cut = ml.cut, labels = labels, col = col,  pch = pch, pt.bg = pt.bg)
  }
  
}