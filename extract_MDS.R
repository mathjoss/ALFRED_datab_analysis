
library(ape)
library(ggplot2)
# -------------- DATASET -----------------

### set path
setwd("~/Documents/FST")

# read final matrix
data <- read.table("mat_mean.csv", header=1, check.names=FALSE, row.names=1, sep=",", quote='"', fill=TRUE, stringsAsFactors = FALSE) 
#data <- read.table("qlc_data_dist.csv", header=1, check.names=FALSE, row.names=1, sep=",", quote='"', fill=TRUE, stringsAsFactors = FALSE) 

# convert dataframe to matrix
data <- as.matrix(data)

# rename columns???
col <- gsub("\\s*\\([^\\)]+\\)","",as.character(colnames(data)))
colnames(data) <- col
rownames(data) <- col

# delete upper part of table
data[upper.tri(data)] <- NA

# replace missing data
d1 <- as.dist(additive(data)) # additive procedure
d2 <- as.dist(ultrametric(data)) # ultrametric procedure             

# convert to dist
d <- as.dist(data)

# --------------- Multidimensional scaling MDS ----------------

# select table
tab_to_analyze <- d

# plot dimension/rsquared
vect_rsquared <- c()
for (dim in 1:50){
  mds <- cmdscale(tab_to_analyze, eig=TRUE, k=dim)
  mdscoords = mds$points #take coords from cmdscale solution
  mdsdists <- dist(mdscoords, diag=TRUE, upper=TRUE)
  rsquared <- (cor(c(tab_to_analyze), c(mdsdists))) * (cor(c(tab_to_analyze), c(mdsdists)))
  vect_rsquared <- c(vect_rsquared, rsquared)
}

vec_gof <- c()
for (dim in 1:50){
  mds <- cmdscale(tab_to_analyze, eig=TRUE, k=dim)
  vec_gof <- c(vec_gof, mds$GOF[1])
}

vec_gof2 <- c()
for (dim in 1:50){
  mds <- cmdscale(tab_to_analyze, eig=TRUE, k=dim)
  vec_gof2 <- c(vec_gof2, mds$GOF[2])
}

# plot variance
df <- data.frame(dim=c(1:50), rsquared = vect_rsquared)
ggplot(df, aes(x=dim, y=rsquared)) + geom_point()

df2 <- data.frame(dim=c(1:50), gof = vec_gof)
ggplot(df, aes(x=dim, y=vec_gof)) + geom_point()

df3 <- data.frame(dim=c(1:50), gof = vec_gof2)
ggplot(df, aes(x=dim, y=vec_gof2)) + geom_point()


# Perform MDS on the data using cmdscale; k determines the number of dimensions
mds <- cmdscale(tab_to_analyze, eig=TRUE, k=4)

# Copy coordinates from MDS solution into a new matrix
mdscoords = mds$points #take coords from cmdscale solution

#Compute a new set of distances from the MDS coordinates delivered by cmdscale
mdsdists <- dist(mdscoords, diag=TRUE, upper=TRUE)

#String the two distance matrices out into two columns and calculate the correlation between these two columns
rsquared <- (cor(c(tab_to_analyze), c(mdsdists))) * (cor(c(tab_to_analyze), c(mdsdists)))
rsquared 
mds$GOF
# r squared is quite low... and so is goodness of fit


# % of variance explained by the MDS axes
round(mds$eig*100/sum(mds$eig),1)

# plot to understand variance??
spear <- round(cor(tab_to_analyze, dist(mds$points), method = "spearman"),3)
plot(tab_to_analyze, dist(mds$points), main = "Shepard diagram of MDS", 
     xlab = "True Bray-Curtis distance", ylab = "Distance in the reduced space")
mtext(line = 0.1, text = paste0("Spearman correlation = ", spear), cex = 0.7)


# plot results
plot(mdscoords[,1], mdscoords[,2], type="n", xlab="", ylab="")
text(jitter(mdscoords[,1]), jitter(mdscoords[,2]),
     rownames(mdscoords), cex=0.8)
abline(h=0,v=0,col="gray75")

# plot tree0
#plot(hclust(dist(1-tab_to_analyze), method="single"))
rownames(mdscoords) <- col

# write table
write.csv(mdscoords, "/home/mathilde/Documents/FST/mdscoords_additive_4.csv", row.names = TRUE)

