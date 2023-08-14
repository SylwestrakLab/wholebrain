#Load Raster package to flip image
library(raster)


#Load the prevous dataset (to get prior cells removed)
load(.....Rdata)



#Retroactively add in an field to "seg" to keep track of false positives
idxKeep = as.numeric(row.names(dataset))
somas = logical(length = length(seg$soma$x))
somas[idxKeep]=TRUE
seg$truePos = somas

#Also keep track of what had previously been a false negative
somas = logical(length = length(seg$soma$x))
seg$falseNeg = somas

#Get current number of soma in seg
orig<- length(seg$soma$x)

#Load Image and Flip
b<-raster::raster(imagepath)
b<-flip(b,'y')

#Plot Previous points, excluding false positives
quartz(width = 20, height = 14)
plot(b,col = gray.colors(maxValue(b)),legend=FALSE)
points(seg$soma$x[seg$truePos], seg$soma$y[seg$truePos],cex=0.1)

#Select Zoom Area and then Select a new Point
zoom(b, ext=drawExtent(), maxpixels=100000, layer=1, new=TRUE, useRaster=TRUE, col = gray.colors(maxValue(b)))
legends_coord <- locator(1)

#Add the point to the existing seg data and indicate that it was a false negative, but now a known positive.
seg$soma$x[orig+1]<-legends_coord$x
seg$soma$y[orig+1]<-legends_coord$y
seg$truePos[orig+1]<-TRUE
seg$falseNeg[orig+1]<-TRUE

#Seems necessary to add values for area and intensity too.  
seg$soma$intensity[orig+1]<-20000
seg$soma$area[orig+1]<-5

#Rerun the registration
regi<-registration(imagepath, coordinate= dataset$AP[1], filter=seg$filter, correspondance = regi)
idxRemove = which(!seg$truePos, arr.ind = FALSE, useNames = TRUE)
dataset<-inspect.registration(regi, seg, forward.warps = TRUE,removePoints=idxRemove)



