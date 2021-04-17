#data generation
data.generate<-function(I,J,L,B,K,gamma_O,p_cluster,alpha_M,alpha_F,beta_M,beta_F,beta_O,beta_M_B,beta_F_B,beta_O_B){
  cluster_assign=rep(0,J)
  #initialize alpha_O
  alpha_O=rep(0,J)
  Z1=matrix(0,I,J)
  Z2=matrix(0,I,J)
  y=matrix(0,I,J)
  
  for (b in 1:B){
    alpha_M_sub=alpha_M[(1+(b-1)*L):(b*L)]
    beta_M_sub=beta_M[(1+(b-1)*L):(b*L)]
    alpha_F_sub=alpha_F[(1+(b-1)*L):(b*L)]
    beta_F_sub=beta_F[(1+(b-1)*L):(b*L)]
    beta_O_sub=beta_O[(1+(b-1)*L):(b*L)]
    
    cluster_assign_sub=sample(x=1:K,size=L,replace=TRUE,prob=p_cluster)
    temp_A=gamma_O[cluster_assign_sub,]
    temp_B=cbind(rep(1,L),log(alpha_M_sub)-log(beta_M_sub),log(alpha_F_sub)-log(beta_F_sub))
    alpha_O_sub=exp(rowSums(temp_A*temp_B)+log(beta_O_sub))
    rm(temp_A,temp_B)
    X_M_sub=matrix(rgamma(I*L,alpha_M_sub,rep(1,L)),nrow=I,ncol=L,byrow=TRUE)
    # the generation of Y_M_sub is a bit tricky, since s*gamma(alpha,beta)=gamma(s*alpha,beta),
    # we can
    Y_M_sub=matrix(rgamma(I,beta_M_B[b],rep(1,I)),nrow=I,ncol=L,byrow=FALSE) 
    #check:
    #colMeans(X_M_sub)
    #alpha_M_sub
    Z1_sub=X_M_sub/(X_M_sub+Y_M_sub)
    #colMeans(Z1)
    #alpha_M/(beta_M+alpha_M)
    
    X_F_sub=matrix(rgamma(I*L,alpha_F_sub,rep(1,L)),nrow=I,ncol=L,byrow=TRUE)
    Y_F_sub=matrix(rgamma(I,beta_F_B[b],rep(1,I)),nrow=I,ncol=L,byrow=FALSE)
    Z2_sub=X_F_sub/(X_F_sub+Y_F_sub)
    
    X_O_sub=matrix(rgamma(I*L,alpha_O_sub,rep(1,L)),nrow=I,ncol=L,byrow=TRUE)
    Y_O_sub=matrix(rgamma(I,beta_O_B[b],rep(1,I)),nrow=I,ncol=L,byrow=FALSE)
    y_sub=X_O_sub/(X_O_sub+Y_O_sub)
    
    cluster_assign[((b-1)*L+1):(b*L)]=cluster_assign_sub
    
    alpha_O[((b-1)*L+1):(b*L)]=alpha_O_sub
    
    Z1[,((b-1)*L+1):(b*L)]=Z1_sub
    Z2[,((b-1)*L+1):(b*L)]=Z2_sub
    y[,((b-1)*L+1):(b*L)]=y_sub
  }
  res=list(cluster_assign,alpha_O,Z1,Z2,y)
  return (res)
}




#data correlation
#use very large I to find out the correlation in neighbering CpG sites
#specifically, the correlation between X_1/(X_1+Y) and X_2/(X_2+Y), 
#where X_1 ~ gamma(alpha_1,1), X_2 ~ gamma(alpha_2,1), Y ~ gamma(beta,1)
data.correlation<-function(I,alpha_1,alpha_2,beta){
  X1=rgamma(I,alpha_1,1)
  X2=rgamma(I,alpha_2,1)
  Y=rgamma(I,beta,1)
  res=cor(X1/(X1+Y),X2/(X2+Y))
  return (res)
}

colCors = function(x, y) {
  sqr = function(x) x*x
  if(!is.matrix(x)||!is.matrix(y)||any(dim(x)!=dim(y)))
    stop("Please supply two matrices of equal size.")
  x   = sweep(x, 2, colMeans(x))
  y   = sweep(y, 2, colMeans(y))
  cor = colSums(x*y) /  sqrt(colSums(sqr(x))*colSums(sqr(y)))
  return(cor)
}



fit.D.realdata<-function(cluster.data,K,initial,start.loc.B,end.loc.B,mle.est,max.iter){
  
  B=length(start.loc.B)
  
  #constant variables
  y  = cluster.data[[1]]
  Z1 = cluster.data[[2]]
  Z2 = cluster.data[[3]]
  
  I=nrow(y)
  J=ncol(y)
  
  
  mle.alpha.M = mle.est[[1]]
  mle.beta.M  = mle.est[[2]]
  mle.alpha.F = mle.est[[3]]
  mle.beta.F  = mle.est[[4]]
  mle.alpha.O = mle.est[[5]]
  mle.beta.O = mle.est[[6]]
  
  mle.M = log(mle.alpha.M)-log(mle.beta.M)
  mle.F = log(mle.alpha.F)-log(mle.beta.F)
  mle.O = log(mle.alpha.O)-log(mle.beta.O)
  
  ini.gamma.O=initial[[1]]
  ini.Pi=initial[[2]]
  ini.S=initial[[3]]
  ini.alpha.O=mle.beta.O*exp(rowSums(cbind(rep(1,J),mle.M,mle.F)*ini.gamma.O[ini.S,]))
  
  gamma.O.old =ini.gamma.O
  Pi.old      =ini.Pi
  S.old       =ini.S
  alpha.O.old =ini.alpha.O
  
  gamma.O.new =ini.gamma.O
  Pi.new      =ini.Pi
  S.new       =ini.S
  alpha.O.new =ini.alpha.O
  
  llk.old=-Inf
  record.llk=rep(-Inf,max.iter)
  
  #some pre-calculation to speed up
  log_y=log(y) #I*L
  log_y_colsum=colSums(log_y)
  log_1minusy=log(1-y) #I*L
  log_1minusy_colsum=colSums(log_1minusy)
  logit_y=log(y/(1-y))
  logit_y_colsum=colSums(logit_y)
  log_1pluslogit=matrix(0,nrow=I,ncol=B)
  
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    if (length(index)>1){
      log_1pluslogit[,b]=log(rowSums(y[,index]/(1-y[,index]))+rep(1,I))
    }else{
      log_1pluslogit[,b]=log(y[,index]/(1-y[,index])+rep(1,I))
    }
  }
  log_1pluslogit_colsum=colSums(log_1pluslogit)
  cond=1
  iter=0
  flag=0 #20200312 add
  count.restart=0 ##20200312 add
  while (cond){
    #####20200312 begins
    # if S.new all equals one status, find a new start, this happens when ini.S gives a large intercept of certain cluster
    if (flag){  
      start_cluster_assign=matrix(sample(x=1:K,size=J,replace=TRUE,prob=rep(1/K,K))) #starting point
      start_gamma_O=matrix(rnorm(3*K,0,1),nrow=K,ncol=3)#gamma_O #
      initial=list(start_gamma_O,rep(1/K,K),start_cluster_assign)
      ini.gamma.O=initial[[1]]
      ini.Pi=initial[[2]]
      ini.S=initial[[3]]
      ini.alpha.O=mle.beta.O*exp(rowSums(cbind(rep(1,J),mle.M,mle.F)*ini.gamma.O[ini.S,]))
      gamma.O.new =ini.gamma.O
      Pi.new      =ini.Pi
      S.new       =ini.S
      alpha.O.new =ini.alpha.O
      llk.old=-Inf
      record.llk=rep(-Inf,max.iter)
      count.restart=count.restart+1
      flag=0
      iter=0 #reset everyting
    }
    #####20200312 end
    iter=iter+1
    #E-step
    alpha.O.new=mle.beta.O*exp(rowSums(cbind(rep(1,J),mle.M,mle.F)*gamma.O.new[S.new,])) # gai le
    
    Loop=1
    record=matrix(0,nrow=Loop,ncol=J)
    for (loop in 1:Loop){
      for (b in 1:B){
        index=start.loc.B[b]:end.loc.B[b]
        for (l in 1:(end.loc.B[b]-start.loc.B[b]+1)){
          #to avoid repeated computation, only calculate the difference based on llk.
          loc=start.loc.B[b]+l-1 #location
          temp.weight=rep(0,K)
          for (k in 1:K){
            alpha.O.new[loc]=mle.beta.O[loc]*exp(sum(c(1,mle.M[loc],mle.F[loc])*gamma.O.new[k,]))
            temp.weight[k]=
              I*lgamma(sum(alpha.O.new[index])+mle.beta.O[index][1])+
              sum(log_y_colsum[index]*(alpha.O.new[index]-1))-
              sum(log_1minusy_colsum[index]*(alpha.O.new[index]+1))-
              I*lgamma(mle.beta.O[index][1])-I*sum(lgamma(alpha.O.new[index]))-
              (sum(alpha.O.new[index])+mle.beta.O[index][1])*log_1pluslogit_colsum[b]+
              log(Pi.new[k])
          }
          temp.weight=temp.weight-max(temp.weight) #in case too large
          temp.exp.weight=exp(temp.weight)
          S.new[loc]=sample(x=1:K,size=1,replace=TRUE,prob=temp.exp.weight/sum(temp.exp.weight))
          alpha.O.new[loc]=mle.beta.O[loc]*exp(sum(c(1,mle.M[loc],mle.F[loc])*gamma.O.new[S.new[loc],])) #update alpha_O
        }
      }
      record[loop,]=S.new
    }
    #####20200312 add
    if (length(table(S.new))<K) flag=1
    if (flag==0){
      #####20200312 end
      
      #M-step
      #maximize Pi
      for (k in 1:K){
        Pi.new[k]=sum(S.new==k)/length(S.new)
      }
      if (sum(Pi.new==0)){ # gai le
        Pi.new[which(Pi.new==0)]=0.01
        Pi.new=Pi.new/sum(Pi.new) #if there is one cluster has no assignments, change the probability to 0.01 in case of crash
      } 
      
      #maximize gamma.O
      #ini.gamma.O.Mstep=gamma.O.new
      ini.gamma.O.Mstep=matrix(0,nrow=K,ncol=3)
      for (k in 1:K){
        if (sum(S.new==k)>=10){ #gai le
          ini.gamma.O.Mstep[k,]=lm(mle.O[S.new==k]~mle.M[S.new==k]+mle.F[S.new==k])$coefficients
        }else{ #in case 0
          ini.gamma.O.Mstep[k,]=rnorm(3)#gamma.O.new[k,]
        }
      }
      
      Mstep.gamma<-function(param){ #param is coefficient matrix of gamma.O
        res=0
        alpha.O.Mstep=mle.beta.O*exp(rowSums(cbind(rep(1,J),mle.M,mle.F)*param[S.new,]))  
        for (b in 1:B){
          index=start.loc.B[b]:end.loc.B[b]
          res=res+
            I*lgamma(sum(alpha.O.Mstep[index])+mle.beta.O[index][1])+
            sum(log_y_colsum[index]*(alpha.O.Mstep[index]-1))-
            sum(log_1minusy_colsum[index]*(alpha.O.Mstep[index]+1))-
            I*lgamma(mle.beta.O[index][1])-I*sum(lgamma(alpha.O.Mstep[index]))-
            (sum(alpha.O.Mstep[index])+mle.beta.O[index][1])*log_1pluslogit_colsum[b]
        }
        return (res)
      }
      
      
      
      gamma.O.new=ini.gamma.O.Mstep
      llk.new=Mstep.gamma(gamma.O.new)+sum(log(Pi.new[S.new]))
      
      cat('iteration: ', iter, 'new - old likelihood=',llk.new-llk.old,'\n')
      cond= as.logical(abs(llk.new-llk.old)>1e-7) & as.logical(iter<max.iter)
      # cat('Gibbs sampler trapped or not:',length.record,'\n')
      llk.old=llk.new
      cat('likelihood.new is:',llk.new,'\n')
      record.llk[iter]=llk.new
    }#####20200312 add
  }
  return (list(S.new,gamma.O.new,iter,llk.new,record.llk))
} 


llk.Z1Z2.realdata<-function(cluster.data,mle.est,start.loc.B,end.loc.B){
  B=length(start.loc.B)
  y  = cluster.data[[1]]
  Z1 = cluster.data[[2]]
  Z2 = cluster.data[[3]]
  I=nrow(y)
  J=ncol(y)
  
  #M
  #some pre-calculation to speed up
  log_Z1=log(Z1) #I*L
  log_Z1_colsum=colSums(log_Z1)
  log_1minusZ1=log(1-Z1) #I*L
  log_1minusZ1_colsum=colSums(log_1minusZ1)
  logit_Z1=log(Z1/(1-Z1))
  logit_Z1_colsum=colSums(logit_Z1)
  log_1pluslogit=matrix(0,nrow=I,ncol=B)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    if (length(index)>1){
      log_1pluslogit[,b]=log(rowSums(Z1[,index]/(1-Z1[,index]))+rep(1,I))
    }else{
      log_1pluslogit[,b]=log(Z1[,index]/(1-Z1[,index])+rep(1,I))
    }
  }
  
  
  
  
  
  
  
  log_1pluslogit_colsum=colSums(log_1pluslogit)
  llk.M=0
  mle.alpha.M=mle.est[[1]]
  mle.beta.M=mle.est[[2]]
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    llk.M=llk.M+
      I*lgamma(sum(mle.alpha.M[index])+mle.beta.M[index][1])+
      sum(log_Z1_colsum[index]*(mle.alpha.M[index]-1))-
      sum(log_1minusZ1_colsum[index]*(mle.alpha.M[index]+1))-
      I*lgamma(mle.beta.M[index][1])-I*sum(lgamma(mle.alpha.M[index]))-
      (sum(mle.alpha.M[index])+mle.beta.M[index][1])*log_1pluslogit_colsum[b]
  }
  
  
  #F
  #some pre-calculation to speed up
  log_Z2=log(Z2) #I*L
  log_Z2_colsum=colSums(log_Z2)
  log_1minusZ2=log(1-Z2) #I*L
  log_1minusZ2_colsum=colSums(log_1minusZ2)
  logit_Z2=log(Z2/(1-Z2))
  logit_Z2_colsum=colSums(logit_Z2)
  log_1pluslogit=matrix(0,nrow=I,ncol=B)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    if (length(index)>1){
      log_1pluslogit[,b]=log(rowSums(Z2[,index]/(1-Z2[,index]))+rep(1,I))
    }else{
      log_1pluslogit[,b]=log(Z2[,index]/(1-Z2[,index])+rep(1,I))
    }
  }
  log_1pluslogit_colsum=colSums(log_1pluslogit)
  llk.F=0
  mle.alpha.F=mle.est[[3]]
  mle.beta.F=mle.est[[4]]
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    llk.F=llk.F+
      I*lgamma(sum(mle.alpha.F[index])+mle.beta.F[index][1])+
      sum(log_Z2_colsum[index]*(mle.alpha.F[index]-1))-
      sum(log_1minusZ2_colsum[index]*(mle.alpha.F[index]+1))-
      I*lgamma(mle.beta.F[index][1])-I*sum(lgamma(mle.alpha.F[index]))-
      (sum(mle.alpha.F[index])+mle.beta.F[index][1])*log_1pluslogit_colsum[b]
  }
  return (llk.M+llk.F)
}


# ######################################################
#If the block size is larger than 1, then use MLE
#If the block size is equal to 1, then use method of moment to get the estimation
MLE.nuisance.param.B<-function(cluster.data,start.mle.est,start.loc.B,end.loc.B){
  B=length(start.loc.B)
  y  = cluster.data[[1]]
  Z1 = cluster.data[[2]]
  Z2 = cluster.data[[3]]
  I=nrow(y)
  J=ncol(y)
  B.candidate=which((end.loc.B-start.loc.B+1)>1)
  
  start.alpha.M = start.mle.est[[1]]
  start.beta.M  = start.mle.est[[2]]
  start.alpha.F = start.mle.est[[3]]
  start.beta.F  = start.mle.est[[4]]
  start.alpha.O = start.mle.est[[5]]
  start.beta.O  = start.mle.est[[6]]
  
  alpha_M=start.alpha.M
  beta_M =start.beta.M 
  alpha_F=start.alpha.F
  beta_F =start.beta.F 
  alpha_O=start.alpha.O
  beta_O =start.beta.O 
  
  
  
  #Z1: mother
  #some pre-calculation to speed up
  log_Z1=log(Z1) #I*L
  log_Z1_colsum=colSums(log_Z1)
  log_1minusZ1=log(1-Z1) #I*L
  log_1minusZ1_colsum=colSums(log_1minusZ1)
  logit_Z1=log(Z1/(1-Z1))
  logit_Z1_colsum=colSums(logit_Z1)
  log_1pluslogit=matrix(0,nrow=I,ncol=B)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    if (length(index)>1){
      log_1pluslogit[,b]=log(rowSums(Z1[,index]/(1-Z1[,index]))+rep(1,I))
    }else{
      log_1pluslogit[,b]=log(Z1[,index]/(1-Z1[,index])+rep(1,I))
    }
  }
  log_1pluslogit_colsum=colSums(log_1pluslogit)
  
  for (b in B.candidate){
    index=start.loc.B[b]:end.loc.B[b]
    #object function
    mle<-function(theta){
      mle.alpha.M=theta[1:(length(theta)-1)]
      mle.beta.M=theta[length(theta)]
      llk.M=
        I*lgamma(sum(mle.alpha.M)+mle.beta.M[1])+
        sum(log_Z1_colsum[index]*(mle.alpha.M-1))-
        sum(log_1minusZ1_colsum[index]*(mle.alpha.M+1))-
        I*lgamma(mle.beta.M[1])-I*sum(lgamma(mle.alpha.M))-
        (sum(mle.alpha.M)+mle.beta.M[1])*log_1pluslogit_colsum[b]
      
      return (llk.M)
    }
    
    AA=diag(end.loc.B[b]-start.loc.B[b]+2)
    BB=rep(0,end.loc.B[b]-start.loc.B[b]+2)
    res=maxNM(mle,start=c(start.alpha.M[index],start.beta.F[index[1]]),constraints=list(ineqA=AA, ineqB=BB),control=list(printLevel=0))
    alpha_M[index]=res$estimate[1:(end.loc.B[b]-start.loc.B[b]+1)]
    beta_M[index]=res$estimate[end.loc.B[b]-start.loc.B[b]+2]
  }
  
  ## Father
  log_Z2=log(Z2) #I*L
  log_Z2_colsum=colSums(log_Z2)
  log_1minusZ2=log(1-Z2) #I*L
  log_1minusZ2_colsum=colSums(log_1minusZ2)
  logit_Z2=log(Z2/(1-Z2))
  logit_Z2_colsum=colSums(logit_Z2)
  log_1pluslogit=matrix(0,nrow=I,ncol=B)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    if (length(index)>1){
      log_1pluslogit[,b]=log(rowSums(Z2[,index]/(1-Z2[,index]))+rep(1,I))
    }else{
      log_1pluslogit[,b]=log(Z2[,index]/(1-Z2[,index])+rep(1,I))
    }
  }
  log_1pluslogit_colsum=colSums(log_1pluslogit)
  for (b in B.candidate){
    index=start.loc.B[b]:end.loc.B[b]
    mle<-function(theta){
      mle.alpha.F=theta[1:(length(theta)-1)]
      mle.beta.F=theta[length(theta)]
      llk.F=
        I*lgamma(sum(mle.alpha.F)+mle.beta.F[1])+
        sum(log_Z2_colsum[index]*(mle.alpha.F-1))-
        sum(log_1minusZ2_colsum[index]*(mle.alpha.F+1))-
        I*lgamma(mle.beta.F[1])-I*sum(lgamma(mle.alpha.F))-
        (sum(mle.alpha.F)+mle.beta.F[1])*log_1pluslogit_colsum[b]
      return (llk.F)
    }
    AA=diag(end.loc.B[b]-start.loc.B[b]+2)
    BB=rep(0,end.loc.B[b]-start.loc.B[b]+2)
    res=maxNM(mle,start=c(start.alpha.M[index],start.beta.F[index[1]]),constraints=list(ineqA=AA, ineqB=BB),control=list(printLevel=0))
    alpha_F[index]=res$estimate[1:(end.loc.B[b]-start.loc.B[b]+1)]
    beta_F[index]=res$estimate[end.loc.B[b]-start.loc.B[b]+2]
  }
  
  #y
  #some pre-calculation to speed up
  log_y=log(y) #I*L
  log_y_colsum=colSums(log_y)
  log_1minusy=log(1-y) #I*L
  log_1minusy_colsum=colSums(log_1minusy)
  logit_y=log(y/(1-y))
  logit_y_colsum=colSums(logit_y)
  log_1pluslogit=matrix(0,nrow=I,ncol=B)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    if (length(index)>1){
      log_1pluslogit[,b]=log(rowSums(y[,index]/(1-y[,index]))+rep(1,I))
    }else{
      log_1pluslogit[,b]=log(y[,index]/(1-y[,index])+rep(1,I))
    }
  }
  log_1pluslogit_colsum=colSums(log_1pluslogit)
  for (b in B.candidate){
    index=start.loc.B[b]:end.loc.B[b]
    mle<-function(theta){
      alpha.O.Mstep=theta[1:(length(theta)-1)]
      mle.beta.O=theta[length(theta)]
      llk.O=
        I*lgamma(sum(alpha.O.Mstep)+mle.beta.O[1])+
        sum(log_y_colsum[index]*(alpha.O.Mstep-1))-
        sum(log_1minusy_colsum[index]*(alpha.O.Mstep+1))-
        I*lgamma(mle.beta.O[1])-I*sum(lgamma(alpha.O.Mstep))-
        (sum(alpha.O.Mstep)+mle.beta.O[1])*log_1pluslogit_colsum[b]
      return (llk.O)
    }
    AA=diag(end.loc.B[b]-start.loc.B[b]+2)
    BB=rep(0,end.loc.B[b]-start.loc.B[b]+2)
    res=maxNM(mle,start=c(start.alpha.M[index],start.beta.F[index[1]]),constraints=list(ineqA=AA, ineqB=BB),control=list(printLevel=0))
    alpha_O[index]=res$estimate[1:(end.loc.B[b]-start.loc.B[b]+1)]
    beta_O[index]=res$estimate[end.loc.B[b]-start.loc.B[b]+2]
  }
  res=list(alpha_M,beta_M,alpha_F,beta_F,alpha_O,beta_O)
  return (res)
}

#input
#mo, fa, chd: Mother, father, and child's DNAm data (matrix: number of families * number of CpG sites)
#by:          Based on which dataset's dependence structure to form blocks
#mod:         
#             'distance' 
#                based on the correlation strength of neighbored CpGs to determine the blocks (neighbored correlation>neighbor_corr_cutoff)
#             'correlation' 
#               first apply kmeans method to cluster similar CpGs together, then determine whether the
#               cluster qualifies 

#ouput
#
#
form.block.by_corr<-function(mo,fa,chd,by_dataset='child',num_kmean,kmean.nstart=50,kmean.iter.max=100,block.cutoff=0.7){
  Z2=fa
  Z1=mo
  y=chd
  
  
  if (!all.equal(colnames(mother),colnames(father),colnames(child))) stop("Error: CpG names do not match. Make sure CpG names are the same in mother, father, and child's data")
    
  
  if (by_dataset=='mother') by_data=as.matrix(Z1)
  if (by_dataset=='father') by_data=as.matrix(Z2)
  if (by_dataset=='child')  by_data=as.matrix(y)
  
  M.kmean=kmeans(t(by_data),centers=num_kmean,nstart=kmean.nstart,iter.max=kmean.iter.max)
  Block.ID=M.kmean$cluster 
  Block.table=table(Block.ID)
  
  #reorder the data belonging to the same block together
  y=y[,order(Block.ID)]
  Z1=Z1[,order(Block.ID)]
  Z2=Z2[,order(Block.ID)]
  
  #the beginning and ending index of each block
  end.loc.num_kmean=cumsum(Block.table)
  start.loc.num_kmean=c(1,(cumsum(Block.table)+1)[1:(num_kmean-1)]) #start location for each block
  names(start.loc.num_kmean)=1:num_kmean
  
  #rm(Block.ID)#remove some useless variable
  
  
  
  #filter block of size between 2 and 30
  block.index=which( Block.table>1)
  length(block.index)
  
  #calculate the correlation of 1000 random pairs within a block to measure the correlation strengh of that block
  record_cor=rep(0,length(block.index))
  for (i in 1:length(block.index)){
    temp=rep(0,1000)
    for (j in 1:1000){
      rd.ind=sample(x=start.loc.num_kmean[block.index[i]]:end.loc.num_kmean[block.index[i]],2)
      temp[j]=cor(Z1[,rd.ind[1]],Z1[,rd.ind[2]])
    }
    record_cor[i]=mean(temp)
    # rd.ind=c(start.loc.num_kmean[block.index[i]],end.loc.num_kmean[block.index[i]])
    # record_cor[i]=cor(Z1[,rd.ind[1]],Z1[,rd.ind[2]])
  }
  
  #further filter blocks with correlation stronger than block.cutoff
  block.index=block.index[record_cor>block.cutoff] 
  length(block.index)
  
  
  
  #only treat those in block.index in block, the rest, i.e., individual.index as individuals (L=1)
  #redefine L.vec
  new.L.vec=c()
  for (b in 1:num_kmean){
    if (b %in% block.index){
      new.L.vec=c(new.L.vec,end.loc.num_kmean[b]-start.loc.num_kmean[b]+1)
    }else{
      new.L.vec=c(new.L.vec,rep(1,end.loc.num_kmean[b]-start.loc.num_kmean[b]+1))
    }
  }
  B=length(new.L.vec)
  end.loc.B=cumsum(new.L.vec)
  start.loc.B=c(1,(cumsum(new.L.vec)+1)[1:(B-1)]) #start location for each block
  names(start.loc.B)=NULL
  names(end.loc.B)=NULL
  
  return(list(mother=Z1,father=Z2,child=y,start.index=start.loc.B,end.index=end.loc.B))
  
}





form.block.by_distance<-function(mo,fa,chd,by_dataset='child',map.data,corr_cutoff=0.5){
  Z1=mo
  Z2=fa
  y=chd
  
  
  names(map.data)=c('CpG','CHR','MAPINFO')
  
  if (!all.equal(colnames(mother),colnames(father),colnames(child),map.data$MAPINFO)) stop("Error: CpG names do not match. Make sure CpG names are the same in mother, father, child, and map data")
  
  
  order.mapinfo=order(map.data$MAPINFO)
  
  if (by_dataset=='mother') by_data=as.matrix(Z1[,order.mapinfo])
  if (by_dataset=='father') by_data=as.matrix(Z2[,order.mapinfo])
  if (by_dataset=='child')  by_data=as.matrix(y[,order.mapinfo])
  
  map.data=map.data[order.mapinfo,]
  
  #######################
  ##form blocks
  #######################
  
  
  #!!! here we assume CHR and mapinfo is correct, i.e., if CHR2>CHR1, then mapinfo2>mapinfo1
  ##get start.loc.B and end.loc.B
  start.loc.B=c()
  end.loc.B=c()
  count=0
  for (chr.num in names(table((map.data$CHR)))){
    count=count+1
    tbl=table((map.data$CHR))
    
    which(names(tbl)==chr.num)
    index=(cumsum(tbl)[count]-tbl[count]+1):cumsum(tbl)[count]
    
    
    record_cor=colCors(by_data[,index[-length(index)]],by_data[,index[-1]])
    
    
    consecutive.table=rle(as.integer(record_cor>corr_cutoff))
    block.len=consecutive.table[[1]]
    block.val=consecutive.table[[2]]
    
    temp.loc=1:length(index)
    temp.block.start=rep(0,sum(block.val==1))
    temp.block.end=rep(0,sum(block.val==1))
    for (i in 1:sum(block.val==1)){
      temp.block.start[i]=(cumsum(block.len)[block.val==1]-block.len[block.val==1]+1)[i]
      temp.block.end[i]=cumsum(block.len)[block.val==1][i]+1
    }
    
    temp.start.loc.B=1:length(index)
    temp.end.loc.B=1:length(index)
    if (sum(block.val==1)>0){
      for (i in 1:sum(block.val==1)){
        temp.start.loc.B[temp.block.start[i]:temp.block.end[i]]=temp.block.start[i]
        temp.end.loc.B[temp.block.start[i]:temp.block.end[i]]=temp.block.end[i]
      }
      temp.start.loc.B=unique(temp.start.loc.B)
      temp.end.loc.B=unique(temp.end.loc.B)
    }
    
    
    if (count==1){
      start.loc.B=temp.start.loc.B
      end.loc.B=temp.end.loc.B
    }else{
      start.loc.B=c(start.loc.B,temp.start.loc.B+cumsum(tbl)[count-1])
      end.loc.B=c(end.loc.B,temp.end.loc.B+cumsum(tbl)[count-1])
    }
    
  }
  
  ##end of building blocks
  names(start.loc.B)=NULL
  names(end.loc.B)=NULL
  
  return(list(mother=Z1,father=Z2,child=y,start.index=start.loc.B,end.index=end.loc.B))
  
  
}



GBClust.screeplot<-function(mother,father,child,K_vec=3:6,max.iter=100,block.method='distance',by_dataset='child',map.data,neighbor.corr.cutoff,kmean.num,kmean.nstart=50,kmean.iter.max=50,kmean.block.cutoff=0.7){
  Z1=mother
  Z2=father
  y=child
  I=nrow(y)
  J=ncol(y)

  print('block partitioning...')
  #For block-partition apporach one based on distance
  if (block.method=='distance'){
    
    res=form.block.by_distance(mo=Z1,fa=Z2,chd=y,by_dataset='child',map.data,corr_cutoff=neighbor.corr.cutoff)
  }
  if (block.method=='correlation'){
    res=form.block.by_corr(Z1,Z2,y,by_dataset='child',num_kmean=kmean.num,kmean.nstart=kmean.nstart,kmean.iter.max=kmean.iter.max,block.cutoff=kmean.block.cutoff)
  }
  
  y=res$child
  Z1=res$mother
  Z2=res$father
  
  start.loc.B=res$start.index
  end.loc.B=res$end.index
  B=length(start.loc.B)
  
  print('calculating MLE...')
  ###MLE
  #MOM
  #step 1: calculate the MOM of each CpG, each CpG has a unique estimation of alpha and beta
  start_alpha_M=rep(0,J)
  start_beta_M=rep(0,J)
  start_alpha_F=rep(0,J)
  start_beta_F=rep(0,J)
  start_alpha_O=rep(0,J)
  start_beta_O=rep(0,J)
  
  for (j in 1:J){
    start_alpha_O[j]=mean(y[,j])*(mean(y[,j])*(1-mean(y[,j]))/var(y[,j])-1)
    start_alpha_M[j]=mean(Z1[,j])*(mean(Z1[,j])*(1-mean(Z1[,j]))/var(Z1[,j])-1)
    start_alpha_F[j]=mean(Z2[,j])*(mean(Z2[,j])*(1-mean(Z2[,j]))/var(Z2[,j])-1)
    start_beta_O[j]=(1/mean(y[,j])-1)*start_alpha_O[j]
    start_beta_M[j]=(1/mean(Z1[,j])-1)*start_alpha_M[j]
    start_beta_F[j]=(1/mean(Z2[,j])-1)*start_alpha_F[j]
  }
  
  #step 2: for those CpG in the same block, average their beta's to be the initial value for generalized beta distribution
  MLE_alpha_M = start_alpha_M
  MLE_beta_M=rep(0,J)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    MLE_beta_M[index]=rep(mean(start_beta_M[index]),length(index))
  }
  
  
  MLE_alpha_F = start_alpha_F
  MLE_beta_F=rep(0,J)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    MLE_beta_F[index]=rep(mean(start_beta_F[index]),length(index))
  }
  
  MLE_alpha_O = start_alpha_O
  MLE_beta_O=rep(0,J)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    MLE_beta_O[index]=rep(mean(start_beta_O[index]),length(index))
  }
  start.mle.est=list(MLE_alpha_M,MLE_beta_M,MLE_alpha_F,MLE_beta_F,MLE_alpha_O,MLE_beta_O)
  
  
  cluster.data=list(y,Z1,Z2)
  mle.est=MLE.nuisance.param.B(cluster.data,start.mle.est,start.loc.B,end.loc.B)
  
  print('stochastic EM...')
  
  #calculate the MLE for Z1 and Z2
  llk_Z1Z2=llk.Z1Z2.realdata(cluster.data,mle.est,start.loc.B,end.loc.B)
  llk_Z1Z2
  
  
  
  
  BIC_record=rep(0,length(K_vec))
  
  count_K=0
  for (K in K_vec){
    count_K=count_K+1
    start_cluster_assign=matrix(sample(x=1:K,size=J,replace=TRUE,prob=rep(1/K,K))) #starting point
    start_gamma_O=matrix(rnorm(3*K),nrow=K,ncol=3)#gamma_O #
    
    initial=list(start_gamma_O,rep(1/K,K),start_cluster_assign)
    
    res_D=fit.D.realdata(cluster.data,K,initial,start.loc.B,end.loc.B,mle.est,max.iter)
    BIC_record[count_K]=-2*res_D[[4]]+log(I*J*3)*(J*2+3*B+4*K-1)
    
  }
  
  return (BIC_record)
}
# BIC_record=GBClust.screeplot(Z1,Z2,y,K_vec)
# plot(K_vec,BIC_record,type='l',xlab='',ylab='',main=paste('BIC scree plot'),lty=1,col=1)
# points(K_vec,BIC_record,pch=1,col=1)
# title(ylab="BIC", cex.lab=1.2)
# title(xlab="K", cex.lab=1.2)


GBClust<-function(mother,father,child,K,max.iter=100,block.method='distance',by_dataset='child',map.data,neighbor.corr.cutoff=0.7,kmean.num,kmean.nstart=50,kmean.iter.max=50,kmean.block.cutoff=0.7){
  Z1=mother
  Z2=father
  y=child
  
  I=nrow(y)
  J=ncol(y)
  print('block partitioning...')
  if (block.method=='distance'){
    res=form.block.by_distance(Z1,Z2,y,by_dataset='child',map.data,corr_cutoff=neighbor.corr.cutoff)
  }
  if (block.method=='correlation'){
    res=form.block.by_corr(Z1,Z2,y,by_dataset='child',num_kmean=kmean.num,kmean.nstart=kmean.nstart,kmean.iter.max=kmean.iter.max,block.cutoff=kmean.block.cutoff)
  }
  
  
  y=res$child
  Z1=res$mother
  Z2=res$father
  
  start.loc.B=res$start.index
  end.loc.B=res$end.index
  B=length(start.loc.B)
  
  print('calculating MLE...')
  
  ###MLE
  #MOM
  #step 1: calculate the MOM of each CpG, each CpG has a unique estimation of alpha and beta
  start_alpha_M=rep(0,J)
  start_beta_M=rep(0,J)
  start_alpha_F=rep(0,J)
  start_beta_F=rep(0,J)
  start_alpha_O=rep(0,J)
  start_beta_O=rep(0,J)
  
  for (j in 1:J){
    start_alpha_O[j]=mean(y[,j])*(mean(y[,j])*(1-mean(y[,j]))/var(y[,j])-1)
    start_alpha_M[j]=mean(Z1[,j])*(mean(Z1[,j])*(1-mean(Z1[,j]))/var(Z1[,j])-1)
    start_alpha_F[j]=mean(Z2[,j])*(mean(Z2[,j])*(1-mean(Z2[,j]))/var(Z2[,j])-1)
    start_beta_O[j]=(1/mean(y[,j])-1)*start_alpha_O[j]
    start_beta_M[j]=(1/mean(Z1[,j])-1)*start_alpha_M[j]
    start_beta_F[j]=(1/mean(Z2[,j])-1)*start_alpha_F[j]
  }
  
  #step 2: for those CpG in the same block, average their beta's to be the initial value for generalized beta distribution
  MLE_alpha_M = start_alpha_M
  MLE_beta_M=rep(0,J)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    MLE_beta_M[index]=rep(mean(start_beta_M[index]),length(index))
  }
  
  
  MLE_alpha_F = start_alpha_F
  MLE_beta_F=rep(0,J)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    MLE_beta_F[index]=rep(mean(start_beta_F[index]),length(index))
  }
  
  MLE_alpha_O = start_alpha_O
  MLE_beta_O=rep(0,J)
  for (b in 1:B){
    index=start.loc.B[b]:end.loc.B[b]
    MLE_beta_O[index]=rep(mean(start_beta_O[index]),length(index))
  }
  start.mle.est=list(MLE_alpha_M,MLE_beta_M,MLE_alpha_F,MLE_beta_F,MLE_alpha_O,MLE_beta_O)
  

  cluster.data=list(y,Z1,Z2)
  mle.est=MLE.nuisance.param.B(cluster.data,start.mle.est,start.loc.B,end.loc.B)
  print('stochastic EM...')
  start_cluster_assign=matrix(sample(x=1:K,size=J,replace=TRUE,prob=rep(1/K,K))) #starting point
  start_gamma_O=matrix(rnorm(3*K),nrow=K,ncol=3)#gamma_O #
  
  initial=list(start_gamma_O,rep(1/K,K),start_cluster_assign)
  
  res_D=fit.D.realdata(cluster.data,K,initial,start.loc.B,end.loc.B,mle.est,max.iter)
  
  
  
  temp=res_D[[2]]
  temp_order=order(temp[,1])
  temp_rank=rank(temp[,1])
  #final_CpG_mat[i,]=match(CpG.record[[i]],temp_rank)
  final_CpG_mat=match(res_D[[1]],temp_order)
  final_coef_mat=temp[temp_order,]
  
  cluster.assignment=as.data.frame(cbind(names(y),final_CpG_mat))
  names(cluster.assignment)=c('CpG','Cluster')
  cluster.assignment=cluster.assignment[order(cluster.assignment$CpG),]
  rownames(cluster.assignment)=NULL
  
  Coef=final_coef_mat
  
  return (list(cluster.assignment=cluster.assignment, Coef=Coef))
}
