
### set path
setwd("~/Documents/FST/results")

# packages
library(rlang)
library(stringi)
library(stringr)
library(dplyr)
library(tidyr)

# list all folder inside
list_folders <- list.dirs(path = "~/Documents/FST/results_arp2", full.names = TRUE, recursive = TRUE)
list_folders = list_folders[-1]

# go inside each folder:
for (fol in list_folders){
  # most files
  list_files <- list.files(path = fol, pattern = NULL, all.files = FALSE,
             full.names = FALSE, recursive = FALSE,
             ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
  
  # read 'xml' table
  data <- read.table(paste(fol, '/', list_files[which(endsWith(list_files, '.xml')==TRUE)], sep=""), header=F, sep="\n", quote='"', fill=TRUE, stringsAsFactors = FALSE) 
  
  if (length(grep("pairwiseDifferenceMatrix", data$V1))>0){
    # select the interesting matrix in XML file
    pfm <- grep("pairwiseDifferenceMatrix", data$V1)
    AvPaDi <- data.frame(V1 = data$V1[pfm[1]:pfm[2]])
    
    # convert from factor to character
    AvPaDi$V1 <- as.character(AvPaDi$V1)
    
    # delete white spaces at the beginning of each row
    AvPaDi$V1 <- trimws(AvPaDi$V1, "l")
    
    # find columns index
    column_name <- unlist(strsplit(as.character(AvPaDi[2,]), '        '))
    column_name <- prepend(column_name, "0", before = 1)
    
    # delete useless rows
    AvPaDi <- data.frame(V1 = AvPaDi[c(3:(nrow(AvPaDi)-1)),])
    
    # split column into multiple columns
    AvPaDi %>% separate(V1, column_name, "   ") -> AvPaDi_2
    
    # delete first useless column
    AvPaDi_2[,1] <- NULL
    
    ################
    
    # select the interesting matrix in XML file (population name)
    pn <- grep("pairDistPopLabels", data$V1)
    PopName <- data.frame(V1 = data$V1[pn[1]:pn[2]])
    
    # convert from factor to character
    PopName$V1 <- as.character(PopName$V1)
    
    # delete white spaces at the beginning of each row
    PopName$V1 <- trimws(PopName$V1, "l")
    
    # delete useless rows
    PopName <- data.frame(V1 = PopName[c(4:(nrow(PopName)-1)),])
    
    # split column into multiple columns
    PopName %>% separate(V1, c("Label", "PopName"), ":\t") -> PopName_2
    
    # rename columns AvPaDi
    colnames(AvPaDi_2) <- PopName_2$PopName
    rownames(AvPaDi_2) <- PopName_2$PopName
    
    file_name <- list_files[which(endsWith(list_files, '.xml')==TRUE)]
    file_name = substr(file_name, 1, nchar(file_name)-4)
    
    write.csv(AvPaDi_2, paste("/home/mathilde/Documents/FST/final_tables_arp2/", file_name, ".csv", sep=""), row.names = TRUE)
  } else { print(paste('error file ', fol, sep="")) }
  
}

