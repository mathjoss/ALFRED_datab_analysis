

### set path
setwd("~/Documents/FST")

# Loading DEP and packages required for data handling
library("DEP")
library("dplyr")
library("tidyr")
library("purrr")
library("ggplot2")
library("SummarizedExperiment")
library(stringr)

database.name <- "ALFRED"

#--------------------------------------------
# SPLIT BIG FILE INTO MULTIPLE FILES...
#--------------------------------------------

data_full <- read.table("Microhap_alleleF_198_2.txt", header=F, sep="\t", quote='"', fill=TRUE, stringsAsFactors = FALSE) 
indexes <- which(startsWith(as.character(data_full$V1), '--')==TRUE)

# experimental design
exp.design <- data.frame()

# split datafull in many subdataframes

for (i in 1:(length(indexes)-1)){
  # subset dataframe
  data <- data_full[c((indexes[i]+1):(indexes[i+1]-1)),]
  
  # delete empty columns and empty rows
  data <- data[, colSums(data != "") != 0]
  if (data[1,1] ==""){ data <- data[-1,]}
  
  # take condition name
  name.file <- as.character(data[1,1])
  filen <- gsub(" ", "", name.file, fixed = TRUE)
  filen <- gsub("|", "", filen, fixed = TRUE)
  
  # delete useless first row
  data <- data[-1,]
  
  # add column names by pasting condition+previous column name
  data[1,3:ncol(data)] <- paste(data[1,3:ncol(data)], filen, sep="")
  colnames(data) <- data[1,]
  
  # delete useless first row to get the good dataframe
  data <- data[-1,]
  
  #### do not include this haplotype if there is less than 30 populations. 
  # Why? Because after, we will need to make missing data imputations, and the algorithms do not want more than 80% of missing data.
  if (nrow(data)>30){
    
    # complete experimental design (in order to handle missing values later )
    labels <- colnames(data[3:ncol(data)])
    condition <- rep(filen, length(labels))
    replicate <- 1:length(labels)
    mini.df <- data.frame(labels= labels, condition=condition, replicate=replicate)
    exp.design <- rbind(exp.design, mini.df)
    
    # delete what's inside the parenthesis in population name (to make things easier)
    data$popName <- gsub("\\s*\\([^\\)]+\\)","",as.character(data$popName))
    
    # change columns format to numeric
    for (col in 2:ncol(data)){data[,col] <- as.numeric(data[,col])}
    
    # merge population duplicates:
    # (1) mean duplicates rows for frequency value
    data2 <- aggregate(data, by=list(popName=data$popName), FUN=mean)
    data2 <- data2[,-2]
    # (2) sum duplicates for number of chromosome values
    data3 <- aggregate(data[,c(1,2)]$chrom, by=list(popName=data[,c(1,2)]$popName), FUN=sum)
    # (3) combine the two
    data2$chrom <- data3$x
    
    # increment the big dataset with those subdatafrales
    if (i==1){
      big.df <- data2
    } else {
      # merge dataframes by population name
      big.df <- merge(big.df, data2, by='popName', all.x=TRUE, all.y=TRUE)
      
      # sum chromosomes numbers (previous and new)
      big.df$chrom <- rowSums(big.df[, c("chrom.x", "chrom.y")], na.rm=TRUE)
      
      # delete useless columns
      big.df$chrom.x <- NULL
      big.df$chrom.y <- NULL
    }
  }
}

#--------------------------------------------
# DATA IMPUTATION (with bioconductor)
#--------------------------------------------



# change small things on dataframe to make it compatible (maybe not everything was necessary but it does not hurt)
big.df2 <- big.df
big.df2$ID <- 1:nrow(big.df2)
exp.design$label <- as.character(exp.design$label)
exp.design$condition <- as.character(exp.design$condition)
finaltable <- make_unique(big.df2, "popName", "ID", delim = ";")

# Apparently, the algorithm considers values=0 as NA......... So change the values=0 to value =0.00001... So stupid
finaltable[finaltable==0] <- 0.0000001

# used SummarizedExperiment function
se2 <- make_se(finaltable,  2:(ncol(finaltable)-3), exp.design)

# No filtering
no_filter <- se2

# Filter for proteins that are quantified in all replicates of at least one condition
condition_filter <- filter_proteins(se2, "condition", thr = 0)

# Filter for proteins that have no missing values
complete_cases <- filter_proteins(se2, "complete")

# Filter for proteins that are quantified in at least 2/3 of the samples.
frac_filtered <- filter_proteins(se2, "fraction", min = 0.66)

# Function to extract number of proteins
number_prots <- function(se2) {
  names <- rownames(get(se2))
  data.frame(Dataset = se2,
             bg_proteins = sum(grepl("bg", names)),
             DE_proteins = sum(grepl("DE", names)))
}

# Number of bg and DE proteins still included
objects <- c("no_filter", 
             "condition_filter",
             "complete_cases",
             "frac_filtered")

map_df(objects, number_prots)

# Scale and variance stabilize
no_filter <- normalize_vsn(se2)
#condition_filter <- normalize_vsn(condition_filter)
#complete_cases <- normalize_vsn(complete_cases)
#frac_filtered <- normalize_vsn(frac_filtered)

# Mean versus Sd plot
meanSdPlot(no_filter)

# Plot a heatmap of proteins with missing values
plot_missval(no_filter)

# Plot intensity distributions and cumulative fraction of proteins 
# with and without missing values
plot_detect(no_filter)


### Data imputation

# -------- (1) "Simple" imputation methods -----------

# No imputation
no_imputation <- no_filter

# Impute missing data using random draws from a 
# Gaussian distribution centered around a minimal value (for MNAR)
MinProb_imputation <- impute(no_filter, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a 
# manually defined left-shifted Gaussian distribution (for MNAR)
manual_imputation <- impute(no_filter, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
knn_imputation <- impute(no_filter, fun = "knn", rowmax = 0.9)

# Plot intensity distributions before and after imputation
plot_imputation(no_filter, MinProb_imputation, 
                manual_imputation, knn_imputation)

# ----------- (2) Mixed imputation on populations (rows) -----------------

# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(no_filter) %>%
  group_by(name, condition) %>%
  summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(no_filter) %in% proteins_MNAR

# Perform a mixed imputation
mixed_imputation <- impute(
  no_filter, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn", # imputation function for MAR
  mnar = "zero") # imputation function for MNAR


noimp <- as.data.frame(assay(no_imputation))
knnimp <- as.data.frame(assay(knn_imputation))
mpimp <- as.data.frame(assay(MinProb_imputation))
nanimp <- as.data.frame(assay(manual_imputation))
miximp <- as.data.frame(assay(mixed_imputation))

# write csv
write.csv(knnimp, "/home/mathilde/Documents/FST/knnimp.csv", row.names = TRUE)
write.csv(miximp, "/home/mathilde/Documents/FST/miximp.csv", row.names = TRUE)
write.csv(mpimp, "/home/mathilde/Documents/FST/mpimp.csv", row.names = TRUE)

#--------------------------------------------
# WRITE OUTPUT FILE
#--------------------------------------------

 
print("All tables are gathered inside one dataframe.")
print("Let's convert it into a .arp file...")
freq.tabl <- miximp[,2:(ncol(miximp)-1)]
nb.hapl <- ncol(freq.tabl)

  
outputfile.name <- paste("miximp.arp", sep="")

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

for (pop in 1:nrow(big.df)){
  
  line=paste('\tSampleName="', as.character(big.df$popName[pop]), '"', sep="")
  write(line,file=outputfile.name,append=TRUE)  
  
  line=paste('\tSampleSize=', as.character(big.df$chrom[pop]), sep="")
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

