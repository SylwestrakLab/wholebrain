# This script will load the .Rdata files containing the dataset 
# and determine what AP locations have been analyzed for all cohorts listed
allGenotypes = c('Th','Tac1','ChAT')
dataLocation <- data.frame("AP" = 1.0, "Genotype" = 'Th', "Mouse" = 'm812', "Filename" = 'fname')

for(c in 1:3){
  genotype = allGenotypes[c]
  genotypeFolder <- file.path(getwd(),'Composite Data',genotype,'CompositeData')
  
  #Get Mouse IDs of analyzed files
  mice=0
  fname = list.files(genotypeFolder, full.names=FALSE, pattern = "\\.Rdata$")
  for(i in 1:length(fname)){
    parts <- strsplit(fname[i], "_")[[1]]
    mice[i]<-parts[1]
  }
  
  # Set q to zero so that the atlas image is only plotted once
  q=1
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
        x <- data.frame("AP" = dataset$AP[1], "Genotype" = genotype, "Mouse" = mouse, "Filename" = dataFiles[1])
        dataLocation = rbind(dataLocation,x)
        q=q+1
      
      }
    }
  }
}


filename <- file.path(getwd(),'Composite Data','Figures', 'analyzedSections.Rdata')
save(sections, file=filename)
