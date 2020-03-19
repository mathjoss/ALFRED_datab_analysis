

### set path
setwd("~/Documents/FST")

library(stringr)

# whatever name you want
database.name <- "ALFRED"

# name of the empty folder where you want your output files to be stored
name_folder <- "all_files_arp2"

#--------------------------------------------
# SPLIT BIG FILE INTO MULTIPLE FILES...
#--------------------------------------------

#data_full <- read.table("Microhap_alleleF_198.txt", header=F, sep="\t", quote='"', fill=TRUE, stringsAsFactors = FALSE) 
data_full <- read.table("other_files_from_alfred/ALL.csv", header=F, sep=",", quote='"', fill=TRUE, stringsAsFactors = FALSE) 
indexes <- which(startsWith(as.character(data_full$V1), '--')==TRUE)

list_filenames <- c()

for (i in 1:(length(indexes)-1)){
  data <- data_full[c((indexes[i]+1):(indexes[i+1]-1)),]
  data <- data[, colSums(data != "") != 0]
  name.file <- as.character(data[1,1])
  
  if (name.file %in% list_filenames==FALSE){
    list_filenames <- c(list_filenames, name.file) 
    
    data <- data.frame(data[-1,])
    colnames(data) <- data.frame(data[1,])
    data <- data[-1, ] 
    freq.tabl <- data[,3:ncol(data)]
    nb.hapl <- ncol(freq.tabl)
    
    filen <- gsub(" ", "", name.file, fixed = TRUE)
    filen <- gsub("|", "", filen, fixed = TRUE)
    
    outputfile.name <- paste(name_folder, "/", filen, ".arp", sep="")
    
    #--------------------------------------------
    # WRITE PROFILE 
    #--------------------------------------------
    
    line= paste('\n[Profile]\n\n\tTitle="first trial "\n\n\t#Reference:\n\t#\t', database.name ,' DATABASE\n', sep="")
    write(line,file=outputfile.name,append=TRUE)
    
    line='\tNbSamples=6\n\tDataType=FREQUENCY\n\tGenotypicData=0\n\tLocusSeparator=" "\n\tMissingData="?"\n\tFrequency= REL\n'
    write(line,file=outputfile.name,append=TRUE)
    
    
    
    #--------------------------------------------
    # WRITE DATA (commented) 
    #--------------------------------------------
    
    line='[Data]\n\n#Commented microhap definition\n\n#\t[[HaplotypeDefinition]]\n'
    write(line,file=outputfile.name,append=TRUE)
    
    line='#\t\tHaplListName="nameblabla"\n#\t\tHaplList={'
    write(line,file=outputfile.name,append=TRUE)
    
    for (el in 1:nb.hapl){
      line=paste('#\t\t\t', as.character(el), ' Yap- ',colnames(freq.tabl)[el] , sep="")
      write(line,file=outputfile.name,append=TRUE)
    }
    line='#\t\t}\n\n'
    write(line,file=outputfile.name,append=TRUE)
    
    
    
    #--------------------------------------------
    # WRITE SAMPLES 
    #--------------------------------------------
    
    line='[[Samples]]\n'
    write(line,file=outputfile.name,append=TRUE)
    
    for (pop in 1:nrow(data)){
      line=paste('\tSampleName="', as.character(data$popName[pop]), '"', sep="")
      write(line,file=outputfile.name,append=TRUE)  
      
      line=paste('\tSampleSize=', as.character(data$chrom[pop]), sep="")
      write(line,file=outputfile.name,append=TRUE)  
      
      line='\tSampleData= {'
      write(line,file=outputfile.name,append=TRUE) 
      
      for (col in 1:ncol(freq.tabl)){
        line=paste('\t\t', as.character(col), '\t', as.character(freq.tabl[pop, col], sep=""))
        write(line,file=outputfile.name,append=TRUE) 
      }
      line='\t}\n' 
      write(line,file=outputfile.name,append=TRUE) 
    }
  }
}
