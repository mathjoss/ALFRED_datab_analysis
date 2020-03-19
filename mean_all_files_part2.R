
### set path
setwd("~/Documents/FST/final_tables_arp2")

list_files <- list.files(path = "~/Documents/FST/final_tables_all", pattern = NULL, all.files = FALSE,
                         full.names = FALSE, recursive = FALSE,
                         ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

vec_pop <- c()
# retrieve maximum number of population
for (file in list_files){
  data <- read.table(paste("~/Documents/FST/final_tables_all/", file, sep=""), header=1, check.names=FALSE, row.names=1, sep=",", quote='"', fill=TRUE, stringsAsFactors = FALSE) 
  vec_pop <- c(vec_pop, colnames(data))
}
vec_pop <- unique(vec_pop)


# initialize list that will contain all matrices
multiarray = list();
i=1 # and index for incrementing list

for (file in list_files){
  # read table and convert it to matrix
  data <- read.table(paste("~/Documents/FST/final_tables_all/", file, sep=""), header=1, check.names=FALSE, row.names=1, sep=",", quote='"', fill=TRUE, stringsAsFactors = FALSE) 
  data <- as.matrix(data)

  # find which population are not in data
  col_names <- colnames(data)
  new <- vec_pop[which(!vec_pop %in% col_names)]
  
  if (length(new)!=0){
    # create new matrix with NA row/column for population not there
    empty_row <- matrix(NA, nrow=length(new), ncol=ncol(data))
    d <- list(data, empty_row)
    mm <- do.call(rbind,d)
    empty_col <- matrix(NA, nrow=nrow(mm), ncol=length(new))
    e <- list(mm, empty_col)
    mm <- do.call(cbind,e)

    # rename empty row and column
    colnames(mm)[(ncol(data)+1):ncol(mm)] <- new
    rownames(mm)[(nrow(data)+1):nrow(mm)] <- new
  } else {mm <- data}
  
  # reorder matrix
  mm <- mm[order(as.character(rownames(mm))), order(as.character(colnames(mm)))]

  # store new matrice in list
  multiarray[[i]] = mm;
  i = i+1
  
}

# mean list of matrices
mat_total <- do.call(cbind, multiarray)
mat_total <- array(mat_total, dim=c(dim(multiarray[[1]]), length(multiarray)))
mat_total <- apply(mat_total, c(1, 2), mean, na.rm = TRUE)

# rename column and row name
colnames(mat_total) <- colnames(mm)
rownames(mat_total) <- rownames(mm)



# ---- merge same population in ROWS ---------------

# convert to dataframe
mat_total <- data.frame(mat_total)

# get vector with rownames excluding what's inside the parenthesis
nams <- as.character(rownames(mat_total))
population_names <- gsub("\\s*\\([^\\)]+\\)","", nams)

# add a new row in data with population name vector
mat_total$popName <- population_names

# merge population duplicates based on this new row:
mat_total2 <- aggregate(mat_total, by=list(popName=mat_total$popName), FUN=mean, na.rm=TRUE)

# remove empty columns generated
mat_total2 <- mat_total2[!sapply(mat_total2, function (x) all(is.na(x) | x == ""))]

# rename rows and delete useless row
rownames(mat_total2) <- as.character(mat_total2$popName)
mat_total2$popName <- NULL

# ---- merge same population in COLUMNS ------------

# transpose matrix
mat_total3 = data.frame(t(mat_total2))

# same process...
mat_total3$popName <- population_names
mat_total4 <- aggregate(mat_total3, by=list(popName=mat_total3$popName), FUN=mean, na.rm=TRUE)
mat_total4 <- mat_total4[!sapply(mat_total4, function (x) all(is.na(x) | x == ""))]
mat_total4$popName <- NULL
rownames(mat_total4) <- as.character(colnames(mat_total4))

# ----  Write final CSV :-) youhou!!  ------------

write.csv(mat_total4, "/home/mathilde/Documents/FST/mat_mean.csv", row.names = TRUE)
