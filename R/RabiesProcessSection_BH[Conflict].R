# Initialize the library packages
library(wholebrain)

# Define what section you want to analyze
mouse = "m1107"
#mouse2="m1085_bright"
sectionNumber = '36'
subfolder = paste(mouse, sectionNumber, sep = "-", collapse=NULL)
#subfolder = paste(mouse2, sectionNumber, sep = "-", collapse=NULL)


# Get images
#folder<-file.path(getwd(), mouse,subfolder)
folder<-file.path(getwd(), mouse)
imagefile<-list.files(folder, full.names=TRUE)
imagename<-list.files(folder)


myfilter<-list(alim = c(1, 1000),
               threshold.range = c(17678, 65636),
               eccentricity = 1000,
               Max = 65536, Min = 0,
               brain.threshold = 3000,
               resize = 0.05,
               blur = 12,
               downsample = 0.75)

# Determine Section Position
fromBregma <-     -2.1
#segment
seg<-segment(imagefile[1], filter= seg$filter)
#seg<-segment(imagefile[1], downsample=0.75)

#Registration
quartz()
regi<-registration(imagefile[1], coordinate= fromBregma, filter=seg$filter)

#manual curration and rerun registartion 
regi<-change.corrpoints(regi, 1:32)
regi<-change.corrpoints(regi, 1:3)
regi<-change.corrpoints(regi, 2:4)

regi<-registration(imagefile[1], coordinate= fromBregma, filter=seg$filter, correspondance=regi)

regi<-remove.corrpoints(regi, 44:45)
regi<-registration(imagefile[5], coordinate= fromBregma, filter=seg$filter, correspondance=regi)

regi<- add.corrpoints(regi,10)
regi<-registration(imagefile[2], coordinate= fromBregma, filter=seg$filter, correspondance=regi)

regi<-add.corrpoints(regi, 2)
regi<-registration(imagefile[1], coordinate= fromBregma, filter=seg$filter, correspondance=regi)

#save out segmentation and registration (?)

# Plot Cells on Atlas and create dataset variable
dataset<-inspect.registration(regi, seg, forward.warps = TRUE)

SNR <- get.pixel.intensity(path.to.your.image, dataset$x, dataset$y, type = "SNR", roi = 9, background = 40)
hist(SNR$intensity)

#dataset[dataset$id==0]

# Plotting
bargraph(dataset)
dot.plot(dataset)

#save plots and seg/regi objects
filename <- file.path('BHdata',paste('output_', subfolder,sep = "", collapse=NULL),paste(subfolder, 'Reg.pdf',sep = "", collapse=NULL))
quartz.save(file==filename,type="pdf")

save(seg, regi, dataset, file=file.path('BHdata',mouse,paste('output_', subfolder,sep = "", collapse=NULL),paste(subfolder, '.Rdata',sep = "", collapse=NULL)))

save(seg, regi, dataset,idxRemove, file=file.path('BHdata',mouse,paste(subfolder,sep = "_", collapse=NULL),paste(subfolder, '.Rdata',sep = "_", collapse=NULL)))
