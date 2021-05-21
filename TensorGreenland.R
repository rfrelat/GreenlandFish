# Tensor Decomposition of Greenland data
# Script to reproduce the results presented in the manuscript 
# "Deep demersal fish communities respond to surface climate fluctuations" 
# by Margrete Emblemsvag, Karl Michael Werner, Ismael Núñez-Riboni, 
# Romain Frelat, Helle Torp Christensen, Heino O. Fock and Raul Primicerio.

# The script run the Tensor decomposition on the log-scaled abundance 
# and the clustering of the species based on their spatio-temporal dynamics. 
# Author: R. Frelat, last update: 21.05.2021


## A. Loading and transforming dataset ----------
# Load package
library(ade4)
library(PTAk)
library(RColorBrewer)
library(raster)
library(maps)
library(mapdata)
#source("tools.R")

load("TensorGreenland.Rdata")
dim(logtensor) # 55 species, 18 years, 7 depth layers


#1. Scaling the abundances per species ----------
logscale<-array(0,dim=dim(logtensor))
#Loop scanning each species
for (i in 1:dim(logtensor)[1]){
  #Calculating the mean and sd of the log CPUE for species i
  ma<-mean(logtensor[i,,], na.rm=TRUE)
  sa<-sd(logtensor[i,,], na.rm=TRUE)
  #Saving the anomaly in the array
  logscale[i,,]<-(logtensor[i,,]-ma)/sa
}

#Copy the labels to the new array
dimnames(logscale)<-dimnames(logtensor)


#2. Run a PTA -----------------------------------

pta<-PTA3(logscale, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(pta,testvar = 0)

# Create the scree plot
out <- !substr(pta[[3]]$vsnam, 1, 1) == "*"
gct<-(pta[[3]]$pct*pta[[3]]$ssX/pta[[3]]$ssX[1])[out]
barplot(sort(gct, decreasing = TRUE), xlab="PT",
        ylab="Percentage of variance")

# Choosing the number of PT
# PT that are kept in the analysis
keep <- c(1, 6, 7, 11) 
labkeep <- paste0("PT",1:4, " (", round((100 * (pta[[3]]$d[keep])^2)/pta[[3]]$ssX[1],1), "%)")
nkeep <- length(keep)

# Figure S1: Interpretation of PT ---------------
# Selection of the color palette
colpal<-brewer.pal(9,"BrBG")

# Create the heatmaps
layout(matrix(c(1,3,2,4,5,5), ncol=3), width=c(3,3,1))
#par(mfrow=c(2,2))
par(mar = c(3,0,3,1))
for (i in seq_along(keep)){
  temp <- pta[[3]]$v[keep[i],] %o% pta[[2]]$v[keep[i],]
  dimnames(temp) <- list(dimnames(logtensor)[[3]], dimnames(logtensor)[[2]])
  myHeatmap(temp, pal="BrBG", title=labkeep[i], colscale = FALSE, 
            mary = -3,cex.x = 1, cex.y = 1, rm.empty = TRUE,
            tck.x = -0.02, tck.y = -0.02, padj.x = -0.9, hadj.y = 0.7) 
}
par(mar=c(5,3,5,0))
plot(0, xlim=c(0,4), ylim=c(0.1,8.9), type="n", 
     xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
for (i in 1:9){
  rect(0,i-1,1,i,col = colpal[i])
}
axis(2,at=c(0.5,8.5), labels = c("low", "high"), las=1, cex=0.4)
mtext("Anomaly", side = 2, line=2, xpd=NA, cex=0.8)

#2. Clustering ----------------------------------
#Get the score of species on the 4 PT
coo<-t(pta[[1]]$v[c(keep),])

#Compute the distance between species
dist1=dist(coo, method = "euclidean")

#Build a tree with Ward linkage
den=hclust(dist1,method = "ward.D2")

#Choose the number of clusters
nclust<-6

# Figure S2: Plot the dendogram -----------------
par(mfrow=c(2,1), mar=c(1,3,1,1))
plot(den, hang=-1, ax = T, ann=F, xlab="", sub="",labels = FALSE)

#Create the clusters
clust_3D <- as.factor(cutree(den, k=nclust))
levels(clust_3D) <- c("d1275", "d675-", "d375-","d675+",
                      "d225", "d375+")
clust_3D <- factor(clust_3D, levels = levels(clust_3D)[c(5,3,6,2,4,1)])
#Visualize the cutting
rect.hclust(den, k=nclust, border=rainbow(nclust)[as.numeric(unique(clust_3D[den$order]))])
text(tapply(seq_along(den$order), clust_3D[den$order], median),
     par()$usr[3]-0.05, levels(clust_3D), xpd=NA, cex=0.8)

par(mar=c(4,4,1,1))
# Select the number of clusters based on Within sum of square
res <- sapply(1:20, hc_wss, hc = den, x = coo)


plot(seq_along(res), res, type = "b", pch = 19,
     col=ifelse(seq_along(res)==nclust, "red", "black"),
     xlab="# cluster", ylab="Within sum of square")

# Table S1: Species score and cluster -----------------

# Merge species score and cluster in a table
tabS1 <- data.frame("Species"= gsub("_", " ", dimnames(logtensor)[[1]]),
                  "Cluster"= as.character(clust_3D), 
                  #"Biogeography"=biogeography,
                  "PT1"=round(coo[,1],5),
                  "PT2"=round(coo[,2],5),
                  "PT3"=round(coo[,3],5),
                  "PT4"=round(coo[,4],5))

# Set the order of the rows
ordR <- order(clust_3D, dimnames(logtensor)[[1]])

tabS1 <- tabS1[ordR,]
head(tabS1) # or View(tabS1)


# Figure S3: Species score on PT4 ---------------
par(mfrow=c(1,1), mar=c(4,1,1,1))
# get order on PT4
oPT4 <- order(-coo[,4])
barp <- barplot(-coo[oPT4,4], horiz=TRUE, 
        xlab="PT4 species score", 
        col=rainbow(nclust)[clust_3D[oPT4]])
text(ifelse(coo[oPT4,4]>0, 0.15, -0.15), 
     barp, cex = 0.6,
     labels = gsub("_", " ", dimnames(logtensor)[[1]][oPT4]))


#3. Interpretation of the clusters -----------------

#reconstruct dynamics based on the selected tensor
Xapp <- REBUILD(pta, nTens = keep, testvar = 1)

#Compute the mean spatio-temporal dynamics per cluster
meanAbo<-array(0, dim=c(nclust, dim(logtensor)[c(2,3)]))
for (c in 1:nclust){
  meanAbo[c,,]<-apply(Xapp[as.numeric(clust_3D)==c,,],c(2,3),mean)
}
dimnames(meanAbo)<-list(1:nclust, dimnames(logtensor)[[2]], dimnames(logtensor)[[3]]) 

#Set graphical parameters
ab<-max(abs(quantile(meanAbo,probs = c(0.05, 0.95))))
colcut<-seq(-ab, ab, length.out = 8)
colcut <- sort(c(min(meanAbo) -0.01, colcut, max(meanAbo) +0.01))
colpal<-brewer.pal(9,"BrBG")
nye <- dim(logtensor)[2]
nar <- dim(logtensor)[3]

# Figure 1 --------------------------------------
layout(matrix(c(1:((nclust*2)+2)), ncol=nclust+1), heights = c(1.1, 1, 1))
par(oma=c(5,3,0,0))
for (c in 1:nclust){
  temp <- meanAbo[c,,]
  dimnames(temp)[[1]][-seq(6,nye,10)] <- ""
  if (c!= 1){
    dimnames(temp)[[2]]<- rep("", nar)
  } else {
    dimnames(temp)[[2]] <- dimnames(logtensor)[[3]]
  }
  myHeatmap(t(temp), colscale = FALSE, pal = colpal, 
            title=levels(clust_3D)[c], mary = -5, breaks = colcut, 
            rm.empty = TRUE, tck.x = -0.07, tck.y = -0.05, 
            cex.y = 1, cex.x = 1, padj.x = -0.8, hadj.y=0.5)
  axis(1, at = 1:nye,tck = -0.03, labels = FALSE)
  axis(2,tck = -0.05, labels = FALSE)
  if (c==1){
    mtext("Depth", side = 2, line = 2, xpd=NA, cex=0.8)
  }
  
  par(mar=c(2,0,0,1))
  plot(dimnames(logtensor)[[2]],apply(meanAbo[c,,],1,mean), type = "l", 
       ylim=c(-0.6, 0.5), xaxt="n", yaxt="n")
  axis(1, at = dimnames(logtensor)[[2]], tck = -0.03, labels=FALSE)
  axis(1, at = seq(1998, 2016,10), tck = -0.05, cex=0.9, padj = -0.7)
  if (c==1){
    axis(2, tck = -0.05, las=1, cex=0.9, hadj = 0.7)
  } else {
    axis(2, tck = -0.05, cex=0.9, hadj = 0.5, label=FALSE)
  }
  abline(h=0, lty=3, col="grey")
  if (c==1){
    mtext("Abundance anomaly", side = 2, line = 2, xpd=NA, cex=0.8)
  }
}
par(mar=c(1,3,1,1))
plot(0, xlim=c(0,1), ylim=c(0.1,8.9), type="n", 
     xaxt="n", yaxt="n", bty="n", xlab="", ylab="")
for (i in 1:9){
  rect(0,i-1,1,i,col = colpal[i])
}
axis(2,at=c(0.5,8.5), labels = c("low", "high"), las=1, cex=0.4)
mtext("Anomaly", side = 2, line=2, xpd=NA, cex=0.8)

# Figure 2: Biogeographic affiliation per cluster ---------
#percentage of boreal/artic species per clusters
perc <- t(t(table(biogeography, clust_3D))/as.numeric(table(clust_3D)))

pal <- RColorBrewer::brewer.pal(nrow(perc), "PuBu")

layout(matrix(1:2, ncol=2), width=c(2,1))
par(mar=c(4,4,1,0))
barp <- barplot(perc*100, col=pal, #hadj=0.6, tck=-0.02,
                xlab="Clusters", ylab="% of species", axisnames=FALSE) #,las=2, cex.names = 0.8)
axis(1,barp, levels(clust_3D), 
     lwd=0, padj = -1, cex.axis=0.9)

par(mar=c(4,0,1,0))
plot.new()
legend("center", rev(row.names(perc)), bty="n",
       fill=rev(pal))

# Figure 3: Tensor time series -----------------
year <- dimnames(logtensor)[[2]]
#tensor time series = temporal score on PT4
tstensor <- pta[[2]]$v[keep[4],]
par(mar=c(4,4,1,1), mfrow=c(1,1))
plot(year, -tstensor, type="l", ylab="PT Score")

# Figure 4: Spatial correlation ---------
#Parameters of the plots
#area plotted
ylim <- c(56, 70)
xlim <- c(-50, -10) 
colCountry <- "grey70"
colCoast <- "grey40"
#color palette
pal <- rev(brewer.pal(11, "RdBu"))
#breaks
brk <- seq(-1, 1, length.out = 12)
#aspect ratio to balance long vs lat
aspR <- 1/cos(mean(ylim)*pi/180)

# Map
layout(matrix(c(1:5,5), ncol=3), width=c(2,2,1))
par(mar=c(3,3,1,1), las=1)

# Sea level air temperature
image(slat, xlim=xlim, ylim=ylim, 
      col=pal, breaks=brk, asp=aspR,
      xaxt="n", yaxt="n")
map("world", add=TRUE, fill=TRUE, 
    col=colCountry, border =colCoast)
box()
axis(1, seq(-50, -10, by=10), labels = paste0(seq(50, 10, by=-10), "°W"))
axis(2, seq(55, 70, by=5), labels = paste0(seq(55, 70, by=5), "°N"), 
     las=1)
title("Sea level air temperature")

# Sea ice concentration
image(ice, xlim=xlim, ylim=ylim, 
      col=pal, breaks=brk, asp=aspR,
      xaxt="n", yaxt="n")
map("world", add=TRUE, fill=TRUE, 
    col=colCountry, border =colCoast)
box()
axis(1, seq(-50, -10, by=10), labels = paste0(seq(50, 10, by=-10), "°W"))
axis(2, seq(55, 70, by=5), labels = paste0(seq(55, 70, by=5), "°N"), 
     las=1)
title("Sea ice concentration")

# Sea surface temperature
image(sst, xlim=xlim, ylim=ylim, 
      col=pal, breaks=brk, asp=aspR,
      xaxt="n", yaxt="n")
map("world", add=TRUE, fill=TRUE, 
    col=colCountry, border =colCoast)
box()
axis(1, seq(-50, -10, by=10), labels = paste0(seq(50, 10, by=-10), "°W"))
axis(2, seq(55, 70, by=5), labels = paste0(seq(55, 70, by=5), "°N"), 
     las=1)
title("Sea surface temperature")

# Sea surface salinity
image(sss, xlim=xlim, ylim=ylim, 
      col=pal, breaks=brk, asp=aspR,
      xaxt="n", yaxt="n")
map("world", add=TRUE, fill=TRUE, 
    col=colCountry, border =colCoast)
box()
axis(1, seq(-50, -10, by=10), labels = paste0(seq(50, 10, by=-10), "°W"))
axis(2, seq(55, 70, by=5), labels = paste0(seq(55, 70, by=5), "°N"), 
     las=1)
title("Sea surface salinity")

par(mar=c(10,4,10,1))
plot(0, xlim=c(0.1,1), ylim=c(0.1,length(pal)-0.1), type="n", 
     bty="n", xlab="", ylab="", xaxt="n",yaxt="n")
for (i in 1:length(pal)){
  rect(0,i-1,1,i,col = pal[i])
}
showAx <- round(seq(0,length(pal), length.out=length(brk)))
axis(2, at=showAx, labels = round(brk,2)[showAx+1])
mtext("Correlation (no units)", 3, line = 0.5, adj = 0.6)



