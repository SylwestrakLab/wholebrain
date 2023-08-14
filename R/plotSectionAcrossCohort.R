# This script will load the .Rdata files containing the dataset and regi variables
# it will find all sections that match a particular AP loction and plot the 
# corresponding soma positions onto the common coordinates.  
allGenotypes = c('Th','Tac1','ChAT')
APcoord <- .9


for(c in 1:3){
  genotype = allGenotypes[c]
  genotype = 'Th'
  genotypeFolder <- file.path(getwd(),'Composite Data',genotype,'CompositeData')

  #Get Mouse IDs of analyzed files
  mice=0
  fname = list.files(genotypeFolder, full.names=FALSE, pattern = "\\.Rdata$")
  for(i in 1:length(fname)){
    parts <- strsplit(fname[i], "_")[[1]]
    mice[i]<-parts[1]
  }
  
  quartz()
  
  # Set q to zero so that the atlas image is only plotted once
  q=0
  #for each mouse in the list
  for(i in 1:length(mice)){
    print(i)
    mouse = mice[i]
    animalFolder<-file.path(getwd(),'BHdata', mouse)
    outputFolders<-list.files(animalFolder, full.names=TRUE)
    
    #Loop through folders to check which sections satisfy the AP position
    for(j in 1:length(outputFolders)){
      dataFiles<-list.files(outputFolders[j], full.names=TRUE, pattern = "\\.Rdata$")
      if (length(dataFiles>0)){
        load(dataFiles[1])
        if (dataset$AP[1]==APcoord){
          print(dataFiles[1])
          
          # Get forward warp
          registration<-get.forward.warpRCPP(regi)
          scale.factor<-mean(dim(registration$transformationgrid$mx)/c(registration$transformationgrid$height,registration$transformationgrid$width) )
           
          #if this is the first section found, make the atlas image
          if (q==0){
            xMax<-max(c(registration$transformationgrid$mx,registration$transformationgrid$mxF),na.rm=TRUE)*(1/scale.factor)
            xMin<-min(c(registration$transformationgrid$mx,registration$transformationgrid$mxF),na.rm=TRUE)*(1/scale.factor)
            yMax<-max(c(registration$transformationgrid$my,registration$transformationgrid$myF),na.rm=TRUE)*(1/scale.factor)
            yMin<-min(c(registration$transformationgrid$my,registration$transformationgrid$myF),na.rm=TRUE)*(1/scale.factor)
            plot(c(xMin, xMax), c(yMin, yMax), ylim=c(yMax,yMin), xlim=c(xMin, xMax), asp=1, axes=F, xlab='', ylab='', col=0, main=paste('Bregma:',registration$coordinate,'mm'),font.main = 1
            )
            polygon(c(0,rep(registration$transformationgrid$width,2),0),c(0, 0,rep(registration$transformationgrid$height,2)))
            numPaths<-registration$atlas$numRegions
            outlines<-registration$atlas$outlines
            
            mtext('Dorso-ventral (mm)',side=2,line=1.5)
            mtext('Medio-lateral (mm)',side=1,line=-1.5)
            
            lapply(1:numPaths, function(x){polygon(outlines[[x]]$xr/scale.factor,outlines[[x]]$yr/scale.factor, border='gray', col=as.character(registration$atlas$col[x]) )})
            lapply(1:numPaths, function(x){polygon(outlines[[x]]$xl/scale.factor,outlines[[x]]$yl/scale.factor, border='gray', col=as.character(registration$atlas$col[x]) )})
            
            hei<-dim(registration$transformationgrid$mx)[1]
            wid<-dim(registration$transformationgrid$mx)[2]
            q=1
          }
          
          #Get the soma locations
          index<-round(scale.factor*cbind(dataset$y, dataset$x))
          somaX<-registration$transformationgrid$mxF[index]/scale.factor
          somaY<-registration$transformationgrid$myF[index]/scale.factor
          
          #Plot them over the atlas image
          circle.color<-rep('', nrow(dataset))
          circle.color[dataset$id>0]<-'black'
          circle.color[dataset$id==0]<-'red'
          cex=.5
          points(somaX,somaY,pch=21,bg= dataset$color, col= circle.color, cex=cex)
          
          axis(1, at=stereotactic.coordinates(seq(-4,4,by=0.1),NA,registration, inverse=TRUE)$x, line=-4, labels=FALSE, tck=-0.01, col.ticks='lightblue')
          axis(1, at=stereotactic.coordinates(seq(-4,4,by=0.5),NA,registration, inverse=TRUE)$x, line=-4, labels=FALSE, tck=-0.02, col.ticks='coral')
          axis(1, at=stereotactic.coordinates(c(-4:4),c(0:-6),registration, inverse=TRUE)$x, line=-4, labels=c(-4:4))
          
        }
      }
    }
  }
  
  filename <- file.path(getwd(),'Composite Data','Figures', paste(genotype, '_', APcoord,'.pdf',sep = ""))
  quartz.save(filename,type="pdf")
}
  
  
