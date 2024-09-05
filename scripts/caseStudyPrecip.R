#load libraries
library(readr) #read and download files
library(tidyr) #data wrangling
library(ggplot2) #plotting

#download precipitation data (location: Aachen, Germany) and read it
##link leads to an zip-archive
downLink <- "https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/daily/kl/historical/tageswerte_KL_00003_18910101_20110331_hist.zip"
download.file(downLink, "data.zip")
fileName <- "produkt_klima_tag_18910101_20110331_00003.txt"
dataRaw <- unz("data.zip", fileName) %>% read_delim(delim=";", trim_ws = T)
##delete the downloaded file
file.remove("data.zip")
if(!file.exists("data.zip")){
  print("File succesfully deleted!")
}else{
  print("File could not be deleted!")
}
##only need precipitation height in mm (RSK) and date of measuring (MESS_DATUM)
data <- dataRaw %>% select(MESS_DATUM, RSK) %>% 
  rename(day = "MESS_DATUM")
##
