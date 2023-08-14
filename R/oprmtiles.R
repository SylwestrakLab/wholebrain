#### incomplete. Meant to iterate through image tiles and add points that can be later applied to larger image in wholebrain

library(wholebrain)
library(raster)

### Define animal you want to analyze and pull sections
mouse = "oprm"
sectionNumber = '1'
subfolder = paste(mouse,sectionNumber, sep = "-", collapse=NULL)
subfolder2 = paste(mouse,sectionNumber, sep = "_", collapse=NULL)

folder<-file.path(getwd(), mouse)
folder2<-file.path(getwd(), subfolder2)
imagefile<-list.files(folder, full.names=TRUE)
imagename<-list.files(folder)

options(device="quartz")

######## Count cells ####

# initialize cell count for whole image
all_coord<-c()

#larger process for each image
sections=length(imagefile)
for (i in 1:sections){
b<-raster::raster(imagefile[i])
b<-flip(b,'y')
quartz(width = 20, height = 14)
plot(b,col = gray.colors(maxValue(b)),legend=FALSE)
cells<-readline(prompt = 'Cells Present?:')

#if cells are present if not skip tile
if (cells != 'yes'){
  next
}else{
overestimate<-readline(prompt = 'how many cells ?:')
overestimate<- as.numeric(overestimate)

# click cells in sections
add_coord<-c()
for (i in 1:overestimate){
legends_coord <- locator(1)
points(legends_coord$x,legends_coord$y,cex=0.5)
add_coord$x<-append(add_coord$x,legends_coord$x, after = length(add_coord$x))
add_coord$y<-append(add_coord$y,legends_coord$y, after = length(add_coord$y))
}

#add cells to larger dataset
all_coord$x<-append(all_coord$x,add_coord$x, after = length(add_coord$x))
all_coord$y<-append(all_coord$y,add_coord$y, after = length(add_coord$y))

# 
}

}

####### apply count to wholbrain pipeline ######

# Define what section you want to analyze
mouse = "m1107"
mouse2="m1080_bright"
sectionNumber = '7'
subfolder = paste(mouse,sectionNumber, sep = "-", collapse=NULL)
subfolder2 = paste(mouse,sectionNumber, sep = "_", collapse=NULL)


# Get images
folder<-file.path(getwd(), mouse)
folder2<-file.path(getwd(), subfolder2)
imagefile<-list.files(folder, full.names=TRUE)
imagename<-list.files(folder)


myfilter<-list(alim = c(1, 1000),
               threshold.range = c(17678, 65636),
               eccentricity = 1000,
               Max = 65536, Min = 0,
               brain.threshold = 5791,
               resize = 0.04,
               blur = 3,
               downsample = 0.75)

# Determine Section Position
fromBregma <-  2.4
#segment
seg<-segment(imagefile[1], filter= seg$filter)
#seg<-segment(imagefile[1], downsample=0.75)


#Registration
quartz()
regi<-registration(imagefile[1], coordinate= fromBregma, filter=seg$filter)

#manual curration and rerun registartion 
regi<-change.corrpoints(regi, 1:32)
regi<-change.corrpoints(regi, 38:39)
regi<-change.corrpoints(regi, 15:19)

regi<-remove.corrpoints(regi, 20)
regi<-registration(imagefile[1], coordinate= fromBregma, filter=seg$filter, correspondance=regi)

regi<- add.corrpoints(regi,25)
regi<-registration(imagefile[9], coordinate= fromBregma, filter=seg$filter, correspondance=regi)

regi<-add.corrpoints(regi, 10)
regi<-registration(imagefile[47], coordinate= fromBregma, filter=seg$filter, correspondance=regi)

#save out segmentation and registration (?)

# Plot Cells on Atlas and create dataset variable
dataset<-inspect.registration(regi, seg, forward.warps = TRUE)

# Exclude false positives
#(1) label false positives
dataset<-inspect.registration(regi, seg, forward.warps = TRUE,labelPoints=TRUE, cex=0.001)
#(2) make a vector with points to remove
idxRemove<- c(177,179:181)
idxRemove<- c(dataset$x>10000)

#Remove points from dataset
dataset<-inspect.registration(regi, seg, forward.warps = TRUE,removePoints=idxRemove,labelPoints = TRUE, cex=0.0001)

dataset<-inspect.registration(regi, seg, forward.warps = TRUE,removePoints=idxRemove)


SNR <- get.pixel.intensity(path.to.your.image, dataset$x, dataset$y, type = "SNR", roi = 9, background = 40)
hist(SNR$intensity)

#dataset[dataset$id==0]
dataset<-dataset[dataset$acronym %in% c('fiber tracts', 'MH','PVT','LH'),]
# Plotting
bargraph(dataset)
dot.plot(dataset)

#save plots and seg/regi objects
#filename <- file.path('BHdata',paste('output_', subfolder,sep = "", collapse=NULL),paste(subfolder, 'Reg.pdf',sep = "", collapse=NULL))
#quartz.save(file==filename,type="pdf")


#save(seg,regi,dataset,file=file.path('BHdata',mouse,paste('output',subfolder,sep = "_", collapse=NULL),paste(subfolder, '.Rdata',sep = "", collapse=NULL)))
save(seg, regi, dataset,file=file.path('BHdata',mouse,paste(subfolder2,sep = "", collapse=NULL),paste(subfolder, '.Rdata',sep = "", collapse=NULL)))
