#functions

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



fit.D.Gibbs<-function(cluster.data,initial,L,max.iter,fast=1){
  tt=proc.time()
  
  #constant variables
  y  = cluster.data[[1]]
  Z1 = cluster.data[[2]]
  Z2 = cluster.data[[3]]
  
  I=nrow(y)
  J=ncol(y)
  B=J/L #assume L divides J 
  
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
  #ini.S.B=matrix(ini.S,nrow=B,ncol=L) #block
  ini.alpha.O=mle.beta.O*exp(rowSums(cbind(rep(1,J),mle.M,mle.F)*ini.gamma.O[ini.S,]))
  #ini.alpha.O.B=matrix(ini.alpha.O,nrow=B,ncol=L) #block
  
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
  if (L!=1){
    for (b in 1:B){
      log_1pluslogit[,b]=log(rowSums(y[,L*(b-1)+1:L]/(1-y[,L*(b-1)+1:L]))+rep(1,I))
    }
  }else{
    for (b in 1:B){
      log_1pluslogit[,b]=log(y[,b]/(1-y[,b])+rep(1,I))
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
        index=L*(b-1)+1:L
        for (l in 1:L){
          #to avoid repeated computation, only calculate the difference based on llk.
          loc=L*(b-1)+l #location
          
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
          index=L*(b-1)+1:L
          res=res+
            I*lgamma(sum(alpha.O.Mstep[index])+mle.beta.O[index][1])+
            sum(log_y_colsum[index]*(alpha.O.Mstep[index]-1))-
            sum(log_1minusy_colsum[index]*(alpha.O.Mstep[index]+1))-
            I*lgamma(mle.beta.O[index][1])-I*sum(lgamma(alpha.O.Mstep[index]))-
            (sum(alpha.O.Mstep[index])+mle.beta.O[index][1])*log_1pluslogit_colsum[b]
        }
        return (res)
      }
      
      # Mstep.gamma.grad<-function(param){ #param is coefficient matrix of gamma.O
      #   res=0
      #   alpha.O.true = mle.beta.O*exp(rowSums(cbind(rep(1,J),mle.M,mle.F)*param[cluster_assign,]))
      #   alpha.O.2dev = rep(list(matrix(nrow=3*K,ncol=3*K)),J) #second derivative
      #   alpha.O.1dev = rep(list(matrix(nrow=3*K,ncol=1)),J)   #first derivative
      #   alpha.O.dev.sq = rep(list(matrix(nrow=3*K,ncol=3*K)),J)
      #   for (j in 1:J){
      #     temp.matA=matrix(0,nrow=K,ncol=K)
      #     temp.matA[cluster_assign[j],cluster_assign[j]]=1
      #     temp.matB=matrix(c(1,mle.M[j],mle.F[j]),nrow=3,ncol=1)%*%matrix(c(1,mle.M[j],mle.F[j]),nrow=1,ncol=3)*alpha.O.true[j]^2/mle.beta.O[j]
      #     temp.matC=matrix(c(1,mle.M[j],mle.F[j]),nrow=3,ncol=1)%*%matrix(c(1,mle.M[j],mle.F[j]),nrow=1,ncol=3)*alpha.O.true[j]^2
      #     alpha.O.2dev[[j]]=kronecker(temp.matA,temp.matB)
      #     alpha.O.dev.sq[[j]]=kronecker(temp.matA,temp.matC)
      #     temp.vec=matrix(0,nrow=K,ncol=1)
      #     temp.vec[cluster_assign[j]]=1
      #     alpha.O.1dev[[j]]=kronecker(temp.vec,alpha.O.true[j]*matrix(c(1,mle.M[j],mle.F[j]),nrow=3,ncol=1))
      #   }
      #   #precalculate certain values to avoid repeated computation:
      #   AA=0
      #   BB=0
      #   CC=0
      #   DD=0
      #   for (b in 1:B){
      #     index=L*(b-1)+1:L
      #     #bs.alpha.O.2dev=matrix(0,nrow=3*K,ncol=3*K)
      #     bs.alpha.O.1dev=matrix(0,nrow=3*K,1) #bs=block sum
      #     #bs.alpha.O.2dev.B=matrix(0,nrow=3*K,ncol=3*K)
      #     #bs.alpha.O.dev.sq.B=matrix(0,nrow=3*K,ncol=3*K) #bs=block sum
      #     
      #     for (loc in index){
      #       #bs.alpha.O.2dev=bs.alpha.O.2dev+alpha.O.2dev[[loc]]
      #       bs.alpha.O.1dev=bs.alpha.O.1dev+alpha.O.1dev[[loc]]
      #       BB=BB+digamma(alpha.O.true[loc])*alpha.O.1dev[[loc]]
      #       CC=CC+logit_y_colsum[loc]*alpha.O.1dev[[loc]]
      #     }
      #     AA=AA+digamma(sum(alpha.O.true[index])+mle.beta.O[index][1])*bs.alpha.O.1dev
      #     DD=DD+log_1pluslogit_colsum[b]*bs.alpha.O.1dev 
      #   }
      #   res=I*(AA-BB)+CC-DD
      #   return (res)
      # }
      
      if (fast==1){
        #newtonraphson.new=maxNM(Mstep.gamma,start=ini.gamma.O.Mstep,control=list(printLevel=0))
        #newtonraphson.new=maxNM(fn=Mstep.gamma,start=ini.gamma.O.Mstep,control=list(printLevel=0))
        #gamma.O.new=newtonraphson.new$estimate
        gamma.O.new=ini.gamma.O.Mstep
        
        #llk.new=newtonraphson.new$maximum+sum(log(Pi.new[S.new]))
        llk.new=Mstep.gamma(gamma.O.new)+sum(log(Pi.new[S.new]))
      }else{
        newtonraphson.new=maxNM(Mstep.gamma,start=ini.gamma.O.Mstep,control=list(printLevel=0))
        gamma.O.new=newtonraphson.new$estimate
        llk.new=newtonraphson.new$maximum+sum(log(Pi.new[S.new]))
        
      }
      cat('iteration: ', iter, 'new - old likelihood=',llk.new-llk.old,'\n')
      cond= as.logical(abs(llk.new-llk.old)>1e-7) & as.logical(iter<max.iter)
      # cat('Gibbs sampler trapped or not:',length.record,'\n')
      llk.old=llk.new
      cat('likelihood.new is:',llk.new,'\n')
      record.llk[iter]=llk.new
    }#####20200312 add
  }
  return (list(S.new,gamma.O.new,iter,llk.new,record.llk,proc.time()-tt))
} #


#note: there is some old part that slows downs the speed that has not been as updated as the fid.D.new, that is:
# strongly suggest
# sum(
#   rep(lgamma(sum(alpha.O.new[index])+mle.beta.O[index][1]),I)+
#     rowSums(matrix((alpha.O.new[index]-1),nrow=I,ncol=L,byrow=TRUE)*log_y[,index])-
#     rowSums(matrix((alpha.O.new[index]+1),nrow=I,ncol=L,byrow=TRUE)*log_1minusy[,index])-
#     rep(lgamma(mle.beta.O[index][1]),I)-rep(sum(lgamma(alpha.O.new[index])),I)-
#     rep(sum(alpha.O.new[index])+mle.beta.O[index][1],I)*log_1pluslogit[,b]
# ) +
# TO BE
# I*lgamma(sum(alpha.O.new[index])+mle.beta.O[index][1])+
#   sum(log_y_colsum[index]*(alpha.O.new[index]-1))-
#   sum(log_1minusy_colsum[index]*(alpha.O.new[index]+1))-
#   I*lgamma(mle.beta.O[index][1])-I*sum(lgamma(alpha.O.new[index]))-
#   (sum(alpha.O.new[index])+mle.beta.O[index][1])*log_1pluslogit_colsum[b]+
# But has not verified this.
#instead of using gibbs sampler, directly use multi-dimension probability to calculate
fit.D.noGibbs<-function(cluster.data,initial,L,max.iter,fast=1){
  tt=proc.time()
  
  #constant variables
  y  = cluster.data[[1]]
  Z1 = cluster.data[[2]]
  Z2 = cluster.data[[3]]
  
  I=nrow(y)
  J=ncol(y)
  B=J/L #assume L divides J 
  
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
  #ini.S.B=matrix(ini.S,nrow=B,ncol=L) #block
  ini.alpha.O=mle.beta.O*exp(rowSums(cbind(rep(1,J),mle.M,mle.F)*ini.gamma.O[ini.S,]))
  #ini.alpha.O.B=matrix(ini.alpha.O,nrow=B,ncol=L) #block
  
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
  if (L!=1){
    for (b in 1:B){
      log_1pluslogit[,b]=log(rowSums(y[,L*(b-1)+1:L]/(1-y[,L*(b-1)+1:L]))+rep(1,I))
    }
  }else{
    for (b in 1:B){
      log_1pluslogit[,b]=log(y[,b]/(1-y[,b])+rep(1,I))
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
    Loop=1
    record=matrix(0,nrow=Loop,ncol=J)
    if (L==1){
      all.comb=matrix(1:K,nrow=K,ncol=1)
    }else{
      all.comb=as.matrix(expand.grid(rep(list(1:K),L)))
    }
    for (loop in 1:Loop){
      for (b in 1:B){
        index=L*(b-1)+1:L
        record.joint.prob=rep(0,nrow(all.comb))
        for (ct in 1:nrow(all.comb)){ #count
          block.vec = unlist(all.comb[ct,]) # vector of assignment in block b
          alpha.O.new[index]=mle.beta.O[index]*exp(rowSums(cbind(rep(1,L),mle.M[index],mle.F[index])*gamma.O.new[block.vec,])) # gai le
          
          #alpha.O.new[loc]=mle.beta.O[loc]*exp(sum(c(1,mle.M[loc],mle.F[loc])*gamma.O.new[k,]))
          record.joint.prob[ct]=
            # sum(
            #   rep(lgamma(sum(alpha.O.new[index])+mle.beta.O[index][1]),I)+
            #     rowSums(matrix((alpha.O.new[index]-1),nrow=I,ncol=L,byrow=TRUE)*log_y[,index])-
            #     rowSums(matrix((alpha.O.new[index]+1),nrow=I,ncol=L,byrow=TRUE)*log_1minusy[,index])-
            #     rep(lgamma(mle.beta.O[index][1]),I)-rep(sum(lgamma(alpha.O.new[index])),I)-
            #     rep(sum(alpha.O.new[index])+mle.beta.O[index][1],I)*log_1pluslogit[,b]
            # ) +
            I*lgamma(sum(alpha.O.new[index])+mle.beta.O[index][1])+
            sum(log_y_colsum[index]*(alpha.O.new[index]-1))-
            sum(log_1minusy_colsum[index]*(alpha.O.new[index]+1))-
            I*lgamma(mle.beta.O[index][1])-I*sum(lgamma(alpha.O.new[index]))-
            (sum(alpha.O.new[index])+mle.beta.O[index][1])*log_1pluslogit_colsum[b]+
            sum(log(Pi.new[block.vec]))
        }
        record.joint.prob=record.joint.prob-max(record.joint.prob)  
        record.joint.prob=exp(record.joint.prob)
        S.new[index]=all.comb[sample(x=1:nrow(all.comb),size=1,replace=TRUE,prob=record.joint.prob),]
        alpha.O.new[index]=mle.beta.O[index]*exp(rowSums(cbind(rep(1,L),mle.M[index],mle.F[index])*gamma.O.new[S.new[index],])) #update alpha_O
        
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
          index=L*(b-1)+1:L
          res=res+
            I*lgamma(sum(alpha.O.Mstep[index])+mle.beta.O[index][1])+
            sum(log_y_colsum[index]*(alpha.O.Mstep[index]-1))-
            sum(log_1minusy_colsum[index]*(alpha.O.Mstep[index]+1))-
            I*lgamma(mle.beta.O[index][1])-I*sum(lgamma(alpha.O.Mstep[index]))-
            (sum(alpha.O.Mstep[index])+mle.beta.O[index][1])*log_1pluslogit_colsum[b]
        }
        return (res)
      }
      
      # Mstep.gamma.grad<-function(param){ #param is coefficient matrix of gamma.O
      #   res=0
      #   alpha.O.true = mle.beta.O*exp(rowSums(cbind(rep(1,J),mle.M,mle.F)*param[cluster_assign,]))
      #   alpha.O.2dev = rep(list(matrix(nrow=3*K,ncol=3*K)),J) #second derivative
      #   alpha.O.1dev = rep(list(matrix(nrow=3*K,ncol=1)),J)   #first derivative
      #   alpha.O.dev.sq = rep(list(matrix(nrow=3*K,ncol=3*K)),J)
      #   for (j in 1:J){
      #     temp.matA=matrix(0,nrow=K,ncol=K)
      #     temp.matA[cluster_assign[j],cluster_assign[j]]=1
      #     temp.matB=matrix(c(1,mle.M[j],mle.F[j]),nrow=3,ncol=1)%*%matrix(c(1,mle.M[j],mle.F[j]),nrow=1,ncol=3)*alpha.O.true[j]^2/mle.beta.O[j]
      #     temp.matC=matrix(c(1,mle.M[j],mle.F[j]),nrow=3,ncol=1)%*%matrix(c(1,mle.M[j],mle.F[j]),nrow=1,ncol=3)*alpha.O.true[j]^2
      #     alpha.O.2dev[[j]]=kronecker(temp.matA,temp.matB)
      #     alpha.O.dev.sq[[j]]=kronecker(temp.matA,temp.matC)
      #     temp.vec=matrix(0,nrow=K,ncol=1)
      #     temp.vec[cluster_assign[j]]=1
      #     alpha.O.1dev[[j]]=kronecker(temp.vec,alpha.O.true[j]*matrix(c(1,mle.M[j],mle.F[j]),nrow=3,ncol=1))
      #   }
      #   #precalculate certain values to avoid repeated computation:
      #   AA=0
      #   BB=0
      #   CC=0
      #   DD=0
      #   for (b in 1:B){
      #     index=L*(b-1)+1:L
      #     #bs.alpha.O.2dev=matrix(0,nrow=3*K,ncol=3*K)
      #     bs.alpha.O.1dev=matrix(0,nrow=3*K,1) #bs=block sum
      #     #bs.alpha.O.2dev.B=matrix(0,nrow=3*K,ncol=3*K)
      #     #bs.alpha.O.dev.sq.B=matrix(0,nrow=3*K,ncol=3*K) #bs=block sum
      #     
      #     for (loc in index){
      #       #bs.alpha.O.2dev=bs.alpha.O.2dev+alpha.O.2dev[[loc]]
      #       bs.alpha.O.1dev=bs.alpha.O.1dev+alpha.O.1dev[[loc]]
      #       BB=BB+digamma(alpha.O.true[loc])*alpha.O.1dev[[loc]]
      #       CC=CC+logit_y_colsum[loc]*alpha.O.1dev[[loc]]
      #     }
      #     AA=AA+digamma(sum(alpha.O.true[index])+mle.beta.O[index][1])*bs.alpha.O.1dev
      #     DD=DD+log_1pluslogit_colsum[b]*bs.alpha.O.1dev 
      #   }
      #   
      #   
      #   res=I*(AA-BB)+CC-DD
      #   return (res)
      # }
      
      
      if (fast==1){
        #newtonraphson.new=maxNM(Mstep.gamma,start=ini.gamma.O.Mstep,control=list(printLevel=0))
        #newtonraphson.new=maxNM(fn=Mstep.gamma,start=ini.gamma.O.Mstep,control=list(printLevel=0))
        #gamma.O.new=newtonraphson.new$estimate
        gamma.O.new=ini.gamma.O.Mstep
        
        #llk.new=newtonraphson.new$maximum+sum(log(Pi.new[S.new]))
        llk.new=Mstep.gamma(gamma.O.new)+sum(log(Pi.new[S.new]))
      }else{
        newtonraphson.new=maxNM(Mstep.gamma,start=ini.gamma.O.Mstep,control=list(printLevel=0))
        gamma.O.new=newtonraphson.new$estimate
        llk.new=newtonraphson.new$maximum+sum(log(Pi.new[S.new]))
        
      }
      
      
      cat('iteration: ', iter, 'new - old likelihood=',llk.new-llk.old,'\n')
      cond= as.logical(abs(llk.new-llk.old)>1e-7) & as.logical(iter<max.iter)
      # cat('Gibbs sampler trapped or not:',length.record,'\n')
      llk.old=llk.new
      cat('likelihood.new is:',llk.new,'\n')
      record.llk[iter]=llk.new
    }#####20200312 add
  }
  return (list(S.new,gamma.O.new,iter,llk.new,record.llk,proc.time()-tt))
} #


