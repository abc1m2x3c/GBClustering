# GBClustering
This repository contains R codes and examples for 'Identifying intergenerational patterns of correlated methylation sites'

The necessary source code is src.R. Please download it to your working directory of R or Rstudio.


## A simple example 
Below is an example of R codes to manually input and analyze a DNA methylation data set. You can also find the integrated codes in 'Example.R' in this repository. The DNA data presented in the paper is not publicly available. To demonstrate the use of the program, we simulated the mother, father, and offspring's DNA methylation of 5000 CpG sites from 50 families stored in 'mother.csv', 'father.csv', and 'offspring.csv', sepretately. The true number of clusters equals 4.

We also generated a dataset 'Coordinates.csv' that stores the CpG names, Chromosome numbers, and Chromosomal Coordinates of the 5000 CpGs. Please also download these files to your working directory if you want to test this example.

### Step 1. Clean the memory, install (if has not) and load the required R packages, and source the code in this repository  
```
packages <- c("maxLik")
install.packages(setdiff(packages, rownames(installed.packages())))  

rm(list=ls(all=TRUE))
library(maxLik)
source("src.R")
```
### Step 2. Input the data
```
#load the coordinates data
map.data=read.csv(file="Coordinates.csv")

#load the DNA methylation data. The first column is family ID, which is deleted when loading
child=read.csv(file="offspring.csv")[,-1]
mother=read.csv(file="mother.csv")[,-1]
father=read.csv(file="father.csv")[,-1]

#set the starting seed to make results replicable
set.seed(12345)
```

### Step 3. Data analysis. 
Approach 1: Partition blocks when neighboring CpGs are highly correlated.
First, we need to draw the scree plot to determine the correct number of clusters.
```
K_vec=3:6
#max.iter: the maximum number allowed for stochastic EM algorithm
#K_vec: a vector of K candidates you want to check
#block.method: block-partition method
#by_dataset: Based on which dataset (mother, father, or child) to partition the blocks
#map.data: coodinates information
#neighbor.corr.cutoff: the correlation threshold of neighboring CpGs to form a block
BIC.value.approach1=GBClust.screeplot(mother=mother,father=father,child=child,K_vec,max.iter=100,block.method='distance',by_dataset='child',map.data=map.data,neighbor.corr.cutoff=0.7)

#scree plot
plot(K_vec,BIC.value.approach1,type='l',xlab='',ylab='',main=paste('BIC scree plot'),lty=1,col=1)
points(K_vec,BIC.value.approach1,pch=1,col=1)
title(ylab="BIC", cex.lab=1.2)
title(xlab="K", cex.lab=1.2)
```
In this example, we check the BIC values of K=3,4,5,6. The output is shown below.

![Optional Text](https://github.com/abc1m2x3c/GBClustering/blob/master/Approach1.png)

One can see that the turning point of BIC values is K=4, which is the true K when we simulated the data. Next, we use the following code to identify the patterns of CpG sites.
```
res.approach1=GBClust(mother=mother,father=father,child=child,K=4,max.iter=100,block.method='distance',by_dataset='child',map.data=map.data,neighbor.corr.cutoff=0.7)
```
The cluster assignments can be assessed using the following command.
```
res.approach1$cluster.assignment
```
Below is the first 6 lines:
```
         CpG Cluster
1 cg00012362       3
2 cg00019366       2
3 cg00019759       3
4 cg00026919       4
5 cg00030627       2
6 cg00034130       3
```
The coefficients of each cluster can be obtained by
```
res.approach1$Coef
```
The result is
```
     cluster      gamma0       gamma1       gamma2
[1,]       1 -0.26281500 -0.155635289  0.097613429
[2,]       2 -0.01754879 -0.221697419  0.001473189
[3,]       3  0.14745187  0.006859664 -0.005148698
[4,]       4  0.39935156 -0.008855642  0.119962693
```
One can observe that the cluster 2 is significantly maternally domintaed since the absolute value of gamma1 (0.2217) is much larger than that of gamma2 (0.0015).

Approach 2: Partition blocks based on observed correlation through k-means clustering. The procedure is similar. First, scree plot:
```
K_vec=3:6
#max.iter: the maximum number allowed for stochastic EM algorithm
#K_vec: a vector of K candidates you want to check
#block.method: block-partition method
#by_dataset: Based on which dataset (mother, father, or child) to partition the blocks
#kmean.num: number of kmean clusters. Recommend J/2=5000/2=2500
#kmean.nstart: Parameter for 'kmeans' function. Recommend 50.
#kmean.iter.max: Parameter for 'kmeans' function. Recommend 50.
#kmean.block.cutoff: The threshold of the average pairwise correlation aming CpGs in the same kmean cluster to determine whether the kmean cluster is a block. Recomend 0.7.
BIC.value.approach2=GBClust.screeplot(mother=mother,father=father,child=child,K_vec,max.iter=100,block.method='correlation',by_dataset='child',kmean.num=2500,kmean.nstart=50,kmean.iter.max=50,kmean.block.cutoff=0.7)

#scree plot
plot(K_vec,BIC.value.approach2,type='l',xlab='',ylab='',main=paste('BIC scree plot'),lty=1,col=1)
points(K_vec,BIC.value.approach2,pch=1,col=1)
title(ylab="BIC", cex.lab=1.2)
title(xlab="K", cex.lab=1.2)
```

The output is shown below.

![Optional Text](https://github.com/abc1m2x3c/GBClustering/blob/master/Approach2.png)

Based on the scree plot, we choose K=4. We use the following codes to identify the patterns of CpG sites.
```
res.approach2=GBClust(mother=mother,father=father,child=child,K=4,max.iter=100,block.method='correlation',by_dataset='child',kmean.num=2500,kmean.nstart=50,kmean.iter.max=50,kmean.block.cutoff=0.7)
```
The cluster assignments can be assessed using the following command.
```
res.approach2$cluster.assignment
```
Below is the first 6 lines:
```
         CpG Cluster
1 cg00012362       3
2 cg00019366       2
3 cg00019759       3
4 cg00026919       3
5 cg00030627       2
6 cg00034130       3
```
The coefficients of each cluster can be obtained by
```
res.approach2$Coef
```
The result is
```
     cluster      gamma0       gamma1      gamma2
[1,]       1 -0.26307116 -0.155183857 0.102134392
[2,]       2 -0.02011784 -0.226067167 0.022639365
[3,]       3  0.14480687  0.003063664 0.004519321
[4,]       4  0.39875758 -0.002559992 0.120563266
```
