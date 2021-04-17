#clear the memory
rm(list=ls(all=TRUE))

source('src.R')
library(maxLik)

#load the coordinates data
map.data=read.csv(file="Coordinates.csv")

# load the methylation data of offspring, mother, and father; 
# note that the real data in the paper is not publicly available
# we simulated a data set with K=4 clusters
# the first column is family ID, which is deleted when loading
child=read.csv(file="offspring.csv")[,-1]
mother=read.csv(file="mother.csv")[,-1]
father=read.csv(file="father.csv")[,-1]

set.seed(12345) #set starting seed to make the results the same 

# Approach 1: Partition blocks when neighboring CpGs are highly correlated

# First draw scree plot

K_vec=3:6
BIC.value.approach1=GBClust.screeplot(mother=mother,father=father,child=child,K_vec,max.iter=100,block.method='distance',by_dataset='child',map.data=map.data,neighbor.corr.cutoff=0.7)

plot(K_vec,BIC.value.approach1,type='l',xlab='',ylab='',main=paste('BIC scree plot'),lty=1,col=1)
points(K_vec,BIC.value.approach1,pch=1,col=1)
title(ylab="BIC", cex.lab=1.2)
title(xlab="K", cex.lab=1.2)

res.approach1=GBClust(mother=mother,father=father,child=child,K=4,max.iter=100,block.method='distance',by_dataset='child',map.data=map.data,neighbor.corr.cutoff=0.7)

head(res.approach1$cluster.assignment)
res.approach1$Coef

#Approach 2: Partition blocks based on observed correlation through k-means clustering

K_vec=3:6
BIC.value.approach2=GBClust.screeplot(mother=mother,father=father,child=child,K_vec,max.iter=100,block.method='correlation',by_dataset='child',map.data,neighbor.corr.cutoff=0.7,kmean.num=2500,kmean.nstart=50,kmean.iter.max=50,kmean.block.cutoff=0.7)

plot(K_vec,BIC.value.approach2,type='l',xlab='',ylab='',main=paste('BIC scree plot'),lty=1,col=1)
points(K_vec,BIC.value.approach2,pch=1,col=1)
title(ylab="BIC", cex.lab=1.2)
title(xlab="K", cex.lab=1.2)

res.approach2=GBClust(mother=mother,father=father,child=child,K=4,max.iter=100,block.method='correlation',by_dataset='child',map.data,neighbor.corr.cutoff=0.7,kmean.num=2500,kmean.nstart=50,kmean.iter.max=50,kmean.block.cutoff=0.7)

head(res.approach2$cluster.assignment)
res.approach2$Coef

