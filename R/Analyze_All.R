#Analyze All

#load my functions
setwd("C:/Users/Brian/Documents/R")
source("runGAMS.R")
# source("plotGAMS.R")

#go to data directory
setwd("C:/Data")
miceDirs<-readLines("mice.txt",warn=FALSE,skipNul=TRUE)
numMice<-length(miceDirs)

for(m in 1:numMice)
  {
  print(c("Reading data for: ",miceDirs[m]),quote=FALSE)
  setwd(miceDirs[m])
  
  # read the dates.txt file
  dates<-readLines("dates.txt",warn=FALSE,skipNul=TRUE)
  numDates<-length(dates)
  
  for(d in 1:numDates)
    {
    print(dates[d])
    setwd(paste(miceDirs[m],"\\",dates[d],sep=""))
    inFile<-"infile - Copy.txt"
    
    runGAMS(inFile)
    # plotGAMS()
    
    }
}
