# Initialize the library packages
library(wholebrain)
#library(raster)
#library(zoom)
#library(sp)
#library(rgdal)

# initial registration and remove false positives -------------------------
setwd("C:/Users/bholcom5/Dropbox (University of Oregon)/UO-Sylwestrak Lab/Rabies_registration/wholebrain")


# Define what section you want to analyze
mouse = "m1107"
sectionNumber = '33'
Genotype= 'Tac1-LHb' #Tac1, TH, ChAT, CCK
folder<-file.path(getwd(),'images', mouse)
filePattern <- paste0(mouse,"_", sectionNumber, "c1.tif")
imagefile <- list.files(path = folder, pattern = filePattern, full.names = TRUE)


# Create output folder
subfolder = paste(mouse,sectionNumber, sep = "-", collapse=NULL)
outputfolder<- file.path(getwd(),'BHData',mouse,paste('output_',subfolder,sep = "", collapse=NULL))

if (dir.exists(outputfolder)){
}else{
  dir.create(outputfolder,recursive = TRUE,showWarnings = TRUE)
}


#Load previous filter, or generate new one
if (exists("seg"))  {
myfilter<- seg$filter
}else{
myfilter<-list(alim = c(1, 1000),
               threshold.range = c(30000, 65636),
               eccentricity = 1000,#
               Max = 65536, Min = 0,
               brain.threshold = 5791,
               resize = 0.04,
               blur = 3,
               downsample = 0.75)
}

# Determine Section Position
fromBregma <- -.46
#segment
seg<-segment(imagefile, filter= myfilter)

#Registration
windows()
regi<-registration(imagefile,coordinate= fromBregma, filter=seg$filter)

#manual curration and rerun registartion 
regi<-change.corrpoints(regi,1:20)

regi<-remove.corrpoints(regi, c(21:32))

regi<- add.corrpoints(regi,30)

regi<-registration(imagefile[35], coordinate= fromBregma, filter=seg$filter, correspondance=regi)

#save out segmentation and registration (?)
#setwd(cwd)

# Plot Cells on Atlas and create dataset variable
dataset<-inspect.registration(regi, seg, forward.warps = TRUE)

# Exclude false positives
#(1) label false positives (will not work if there are no cells selected in segmentation)
dataset<-inspect.registration(regi, seg, forward.warps = TRUE,labelPoints=TRUE, cex=.01)
#(2) make a vector with points to remove
idxRemove<- c(1:136,138:148,150:153,155:157,160:400)
idxRemove<- c(dataset$x>10000)

#Remove points from dataset
dataset<-inspect.registration(regi, seg, forward.warps = TRUE,removePoints=idxRemove,labelPoints = TRUE, cex=0.1)

dataset<-inspect.registration(regi, seg, forward.warps = TRUE,removePoints=idxRemove)
dataset$animal=mouse

SNR <- get.pixel.intensity(path.to.your.image, dataset$x, dataset$y, type = "SNR", roi = 9, background = 40)
hist(SNR$intensity)

#dataset[dataset$id==0]
dataset<-dataset[dataset$acronym %in% c('fiber tracts', 'MH','PVT','LH'),]

# reassign location for cell in dataset
dataset$acronym[dataset$acronym =='V3']<-'MPT'
dataset$name[dataset$name =='third ventricle']<-'Medial pretectal area'

# Plotting
bargraph(dataset)
dot.plot(dataset)

#save plots and seg/regi objects
save(seg, regi, dataset,Genotype,idxRemove, file= file.path(outputfolder, paste(subfolder, '.Rdata',sep = "", collapse=NULL)))
#setwd(cwd)

# move distorted images
#base="C:/Users/bholcom5/Dropbox (University of Oregon)/UO-Sylwestrak Lab/Rabies_registration/wholebrain/images"
#imgtemp=imagename[15]
#imgtemp= file.path(base,paste('output_',substr(imgtemp,1,nchar(imgtemp)-4),'/',sep = "", collapse=NULL))
#if (dir.exists(outputfolder)){
#  file.copy(from = imgtemp,to=outputfolder)
#}else{
#}

# add false negatives -----------------------------------------------------

### Add false negatives after initial save

mouse = "m1155"
mouse2="m1080_bright"
sectionNumber = '2'
subfolder = paste(mouse,sectionNumber, sep = "-", collapse=NULL)
subfolder2 = paste(mouse,sectionNumber, sep = "_", collapse=NULL)

# Get images
folder<-file.path(getwd(), mouse)
folder2<-file.path(getwd(), subfolder2)
imagefile<-list.files(folder, full.names=TRUE)
imagename<-list.files(folder)

#Retroactively =add in an field to "seg" to keep track of false positives
idxKeep = as.numeric(row.names(dataset))
somas = logical(length = length(seg$soma$x))
somas[idxKeep]=TRUE
seg$truePos = somas

#Also keep track of what had previously been a false negative
somas = logical(length = length(seg$soma$x))
seg$falseNeg = somas

#Get current number of soma in seg
origWidx<- length(dataset$x)
origWOidx<- length(seg$soma$x)

#Load Image and Flip
b<-raster::raster(imagefile[35])
b<-flip(b,'y')
quartz(width = 20, height = 15)
plot(b,col = gray.colors(maxValue(b)),legend=FALSE)
points(seg$soma$x[seg$truePos], seg$soma$y[seg$truePos],cex=0.1)

#Select Zoom Area and then Select a new Point
options(device="quartz")

#green/yellow
zoom(b, ext=drawExtent(), maxpixels=100000, layer=1, new=TRUE, useRaster=TRUE)
points(seg$soma$x[seg$truePos], seg$soma$y[seg$truePos],cex=0.5)
points(add_coord$x,add_coord$y,cex=0.1)

#grayscale
zoom(b, ext=drawExtent(), maxpixels=100000, layer=1, new=TRUE, useRaster=TRUE, col = gray.colors(maxValue(b)))
points(seg$soma$x[seg$truePos], seg$soma$y[seg$truePos],cex=0.5)
#points(dataset$x, dataset$y,cex=0.3)

### add coordiantes and update plot to prevent accidental duplicates/ overlap
add_coord<-c()
pointnum=1
for (i in 1:pointnum){
legends_coord <- locator(1)
points(legends_coord$x,legends_coord$y,cex=0.5)
add_coord$x<-append(add_coord$x,legends_coord$x, after = length(add_coord$x))
add_coord$y<-append(add_coord$y,legends_coord$y, after = length(add_coord$y))
}

x<-unlist(add_coord$x)
y<-unlist(add_coord$y)

#Seems necessary to add values for area and intensity too.  
#area<-rep(4,length(add_coord$x))
#intensity<-rep(30000,length(add_coord$x))
area<-rep(4,length(x))
intensity<-rep(30000,length(x))

#indicate that it was a false negative, but now a known positive.
truePos<-rep(TRUE,length(add_coord$x))
falseNeg<-rep(TRUE,length(add_coord$x))
seg$truePos<-append(seg$truePos,truePos, after = length(seg$truePos))
seg$falseNeg<-append(seg$falseNeg,falseNeg, after = length(seg$falseNeg))


#Add the point to the existing seg data and indicate that it was a false negative, but now a known positive.
#seg$soma$x<-append(dataset$x,x, after = length(origWidx))
#seg$soma$y<-append(dataset$y,y, after = length(origWidx))
#seg$soma$intensity<-append(dataset$intensity,intensity, after = length(origWidx))
#seg$soma$area<-append(dataset$area,area, after = length(origWidx))

#Add the point to the existing seg data and indicate that it was a false negative, but now a known positive. 
seg$soma$x<-append(dataset$x,x, after = length(origWOidx))
seg$soma$y<-append(dataset$y,y, after = length(origWOidx))
seg$soma$intensity<-append(dataset$intensity,intensity, after = length(origWOidx))
seg$soma$area<-append(dataset$area,area, after = length(origWOidx))

#Rerun the registration
#regi<-registration(imagefile[13], coordinate= dataset$AP[1], filter=seg$filter, correspondance = regi)
fromBregma=fromBregma
regi<-registration(imagefile[64], coordinate= fromBregma, filter=seg$filter, correspondance = regi)

#idxRemove = which(!seg$truePos, arr.ind = FALSE, useNames = TRUE)
#dataset<-inspect.registration(regi, seg, forward.warps = TRUE,removePoints=idxRemove)
dataset<-inspect.registration(regi, seg, forward.warps = TRUE)

save(seg, regi, dataset,Genotype, idxRemove, file= file.path(outputfolder, paste(subfolder, '.Rdata',sep = "", collapse=NULL)))
rm(add_coord,legends_coord,dataset)
setwd(cwd)


