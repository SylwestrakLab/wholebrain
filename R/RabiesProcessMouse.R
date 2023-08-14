#set animal folder path
animal<-file.path(getwd(), mouse,"all")
images<-get.images(animal)

#check number of images in folder
length(images)


#set spacing between periods in millimeters
smp<--2.2


#set images to manually assign position in the atlas.
inspected<-c(1,2)
images[inspected]


#assign brain coordinates for inspected sections.
smp.coord<-c(.14,-2.04 ) #brain coordinates based on bregma from openbrainmap.org


#assign brain coordinates for all sections in this brain.
coord<-map.to.atlas(image.number=inspected, 
                    coordinate=smp.coord, 
                    sampling.period=smp, 
                    number.of.sections=length(images))
coord<-rev(coord)
myfilter<-structure(list(alim = c(22, 118), threshold.range = c(26912, 35020), eccentricity = 900L, Max = 50000, Min = 2000, brain.threshold = 2226,
                         resize = 0.04, blur = 4, downsample = 0.25), .Names = c("alim",
                                                                                 "threshold.range", "eccentricity", "Max", "Min", "brain.threshold",
                                                                                 "resize", "blur", "downsample"))




# Do first section
i = 1
seg<-segment(images[i], filter = myfilter, display=FALSE)
regi<-registration(images[i], coordinate=coord[i], filter = myfilter, display=FALSE)
dev.copy(pdf, paste0(tools::file_path_sans_ext(basename(images[i])), '.pdf'))
dataset<-inspect.registration(regi, seg, soma = TRUE, forward.warps = TRUE, batch.mode = TRUE)
dev.off()
save(file=paste0(tools::file_path_sans_ext(basename(images[i])), '.RData'), seg, regi, dataset)
#use datasets to rbind all dataset to.
datasets<-dataset


#loop through the rest and append them
for(i in seq_along(images)[-1]){
  seg<-segment(images[i], display=FALSE, filter = myfilter)
  regi<-registration(images[i], coordinate=coord[i], display=FALSE, filter = myfilter)
  dev.copy(pdf, paste0(tools::file_path_sans_ext(basename(images[i])), '.pdf'))
  dataset<-inspect.registration(regi, seg, soma = TRUE, forward.warps = TRUE, batch.mode = TRUE)
  dev.off()
  save(file=paste0(tools::file_path_sans_ext(basename(images[i])), '.RData'), seg, regi, dataset)
  datasets<-rbind(datasets, dataset)
}




