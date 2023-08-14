# Initialize the library packages
if (!require(devtools)) {
  install.packages("devtools")
}
devtools::install_github("tractatus/wholebrain")

# Define what section you want to analyze
mouse = "m1107"
sectionNumber = '34'
subfolder = paste(mouse, sectionNumber, sep = "-", collapse=NULL)

# Get images
folder<-file.path(getwd(), mouse,subfolder)
imagefile<-list.files(folder, full.names=TRUE)
imagename<-list.files(folder)

# Determine Section Position
fromBregma <- 1.3


# Segment out neurons and the brain outline from autofluorescence.
seg<-wholebrain::segment(imagefile[1],filter = defaultFilter)

#Registration
quartz()
regi<-registration(imagefile[1], coordinate= fromBregma, filter=seg$filter)

# Add Additional Points
regi<-add.corrpoints(regi, 5)

# Add Additional Points
regi<-change.corrpoints(regi, 30:32)

#Update Registration with new correspondance
regi<-registration(imagefile, coordinate= fromBregma, filter=seg$filter, correspondance = regi)

# Plot Cells on Atlas and create dataset variable
dataset<-inspect.registration(regi, seg, forward.warps = TRUE)

# Save out plots
filename <- file.path(mouse,paste('output_', subfolder,sep = "", collapse=NULL),paste(subfolder, 'Reg.pdf',sep = "", collapse=NULL))
quartz.save(filename,type="pdf")

# Save our R file
save(seg, regi, dataset, file=file.path(mouse,paste('output_', subfolder,sep = "", collapse=NULL),paste(subfolder, '.Rdata',sep = "", collapse=NULL)))
save(seg, regi, dataset, file=file.path(mouse,"all","data",paste(subfolder, '.Rdata',sep = "", collapse=NULL)))

# Plot the combined data set
counts.plot(dataset)
filename <- file.path(mouse,paste('output_', subfolder,sep = "", collapse=NULL),paste(subfolder, '_composite_plot.pdf',sep = "", collapse=NULL))
quartz.save(filename,type="pdf")


