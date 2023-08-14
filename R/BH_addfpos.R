###-------------------
#Shiny is a package for app design. Unsure if I'm able to save coordinates 
library(shiny)

ui <- basicPage(
  plotOutput("plot1", click = "plot_click"),
  verbatimTextOutput("info")
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    plot(mtcars$wt, mtcars$mpg)
  })
  
  output$info <- renderText({
    paste0("x=", input$plot_click$x, "\ny=", input$plot_click$y)
  })

}
###------ 
#x<-c()
#xlen<-length(x)
#append(x,input$plot_click$x,after = xlen)
#y<-c()
#ylen<-length(y)
#append(y,input$plot_click$y,after = ylen)


shinyApp(ui, server)

mydata = reactiveValues(x_values = c(), y_values = c())

observeEvent(input$myclick, {
  
  mydata$x_values = c(mydata$x_values, input$myclick$x)
  mydata$y_values = c(mydata$y_values, input$myclick$y)
  
})

plot(runif(100)) 
legends_coord <- locator(1) 
print(legends_coord) 
legend(x= legends_coord[1], y= legends_coord[2], legend= "First Legend") 

### -------------------------
library(rgdal)
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("EBImage")

# Read Image
library(EBImage)
Image1 <- readImage("~/Desktop/pic1.jpg")
Image2 <- readImage("~/Desktop/pic2.jpg")


# Pull mouse number 
mouse = "m1119"
#mouse2="m1080_bright"
sectionNumber = '55'
subfolder = paste(mouse, sectionNumber, sep = "_", collapse=NULL)
subfolder2 = paste(mouse,sectionNumber, sep = "-", collapse=NULL)


# Get images
folder<-file.path(getwd(), mouse,subfolder)
#folder<-file.path(getwd(), subfolder2)
imagefile<-list.files(folder, full.names=TRUE)
imagename<-list.files(folder)

#show image
quartz()
imshow(imagefile,resize = 0.9)
grid.locator(unit = "native")

#pull coordinates(from locator() not grid.locator)
legends_coord <- locator(1) 
print(legends_coord) 
legend(x= legends_coord[1], y= legends_coord[2], legend= "First Legend") 

b <- brick(imagefile)
quartz()
plotRGB(b)
grid.locator(unit = "native")
###-----

#### This segment is meant to allow the user to add false positives just after the segmentation step
#### From there the user can register selected objects with original objects

# Load Libraries 
library(sp)
library(rgdal)
library(raster)
library(zoom)

#Pull mouse number
mouse = "m1119"
#mouse2="m1080_bright"
sectionNumber = '55'
subfolder = paste(mouse, sectionNumber, sep = "_", collapse=NULL)
subfolder2 = paste(mouse,sectionNumber, sep = "-", collapse=NULL)

# Get images
folder<-file.path(getwd(), mouse,subfolder)
#folder<-file.path(getwd(), subfolder2)
imagefile<-list.files(folder, full.names=TRUE)
imagename<-list.files(folder)

#Determine prexisting cells before adding false poitives
orig<- length(seg$soma$x)

#Conver to raster and plot
b<-raster::raster(imagefile)
quartz()
plot(b)
zoom::zm()
#test1<-grid.locator(unit = "native")
legends_coord <- locator(2)
#turn into numeric instead of list
x<-unlist(legends_coord$x)
y<-unlist(legends_coord$y)

#Create vectors of arbitrary values for area and intensity so 
area<-rep(4,length(legends_coord$x))
intensity<-rep(30000,length(legends_coord$x))
  
#appen numerics
seg$soma$x<-append(seg$soma$x,legends_coord$x, after = length(orig))
seg$soma$y<-append(seg$soma$y,legends_coord$y, after = length(orig))
seg$soma$intensity<-append(seg$soma$intenstity,intensity, after = length(orig))
seg$soma$area<-append(seg$soma$area,area, after = length(orig))









save(seg, regi, dataset,idxRemove,orig,file=file.path('BHdata',mouse,paste('output_',subfolder2,sep = "", collapse=NULL),paste(subfolder, '.Rdata',sep = "", collapse=NULL)))
