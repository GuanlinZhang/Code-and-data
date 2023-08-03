library(gtools)
library(MASS)
library(binaryLogic)
library(mnormt)
library(mvtnorm)
library(FusionLearn)
library(data.table)





set.seed(1)
start.time = proc.time()
n = 800                            #data size
p = 1000                              #number of total predictors
p1 = p
sn = 10                              #size of true model
N = 100                            #number of iterations
M = 10
gammaT = c(rep(1,sn),rep(0,p-sn))   #true active predictors are the first sn ones
sig1 = 1                            #true variance of y1
c1_0 = rep(1,p)                     #tuning factor of variance of true beta1
phox1 = 0.2                         #correlation of x1
sig2 = 1                            #true variance of y1
c2_0 = rep(1,p)                     #tuning factor of variance of true beta2
phox2 = 0.2                         #correlation of x2

sig3 = 1                            #true variance of y1
c3_0 = rep(1,p)                     #tuning factor of variance of true beta1
phox3 = 0.2                         #correlation of x1
sig4 = 1                            #true variance of y1
c4_0 = rep(1,p)                     #tuning factor of variance of true beta2
phox4 = 0.2                         #correlation of x2
corr = 0.6                          #correlation of y1 y2 y3 y4
degree = 2                             #degree of freedom of t-distribution

an = 2*sn
cut_off = 0.75



#prior 
v1 = 2                               #degree of freedom of 1/sig1
v2 = 2                               #degree of freedom of 1/sig2
v3 = 2                               #degree of freedom of 1/sig3
v4 = 2                               #degree of freedom of 1/sig4

Sigma1 = diag(4,p)                   #prior of c1 
Sigma2 = diag(10,p)                  #prior of c2
Sigma3 = diag(4,p)                   #prior of c3 
Sigma4 = diag(10,p)                  #prior of c4


res = 0
res_positive = 0
res_negative = 0


res1 = 0
res1_positive = 0
res1_negative = 0

res0 = 0
res0_positive = 0
res0_negative = 0


res.fusion = 0
res_positive.fusion = 0
res_negative.fusion = 0



#True model
#generate beta1 and x1
rd = sample(1:p,sn)
true.model = rep(0,p)
true.model[rd]=1
beta1_0 = rep(0,p)
for(i in 1:sn){
  beta1_0[rd[i]] = runif(1,0.5,1)
  # beta1_0[i] = rnorm(1,0,c1_0[i]*sig1)
  # if(beta1_0[i]>=0){
  #   beta1_0[i] = beta1_0[i] + 100
  # }else{
  #   beta1_0[i] = beta1_0[i] - 100
  # }
  if(abs(beta1_0[rd[i]])<1/n){i = i-1}
}

sigx1 = array(rep(phox1,p^2),dim=c(p,p))
sigx1 = sigx1+(1-phox1)*diag(p)
mux1 = rep(0,p)
x1 = mvrnorm(n,mux1, sigx1, tol = 1e-6, empirical = FALSE)

#generate beta2 and x2

beta2_0 = rep(0,p)
for(i in 1:sn){
  beta2_0[rd[i]] = runif(1,0.5,1)
  # beta2_0[i] = rnorm(1,0,c2_0[i]*sig2)
  # if(beta2_0[i]>=0){
  #   beta2_0[i] = beta2_0[i] + 100
  # }else{
  #   beta2_0[i] = beta2_0[i] - 100
  # }
  if(abs(beta2_0[rd[i]])<1/n){i = i-1}
}

sigx2 = array(rep(phox2,p^2),dim=c(p,p))
sigx2 = sigx2+(1-phox2)*diag(p)
mux2 = rep(0,p)
x2 = mvrnorm(n,mux2, sigx2, tol = 1e-6, empirical = FALSE)

#generate beta3 and x3

beta3_0 = rep(0,p)
for(i in 1:sn){
  beta3_0[rd[i]] = runif(1,0.5,1)
  # beta3_0[i] = rnorm(1,0,c3_0[i]*sig3)
  # if(beta2_0[i]>=0){
  #   beta2_0[i] = beta2_0[i] + 100
  # }else{
  #   beta2_0[i] = beta2_0[i] - 100
  # }
  if(abs(beta3_0[rd[i]])<1/n){i = i-1}
}

sigx3 = array(rep(phox3,p^2),dim=c(p,p))
sigx3 = sigx3+(1-phox3)*diag(p)
mux3 = rep(0,p)
x3 = mvrnorm(n,mux3, sigx3, tol = 1e-6, empirical = FALSE)

#generate beta4 and x4

beta4_0 = rep(0,p)
for(i in 1:sn){
  beta4_0[rd[i]] = runif(1,0.5,1)
  # beta4_0[i] = rnorm(1,0,c4_0[i]*sig4)
  # if(beta2_0[i]>=0){
  #   beta2_0[i] = beta2_0[i] + 100
  # }else{
  #   beta2_0[i] = beta2_0[i] - 100
  # }
  if(abs(beta4_0[rd[i]])<1/n){i = i-1}
}

sigx4 = array(rep(phox4,p^2),dim=c(p,p))
sigx4 = sigx2+(1-phox4)*diag(p)
mux4 = rep(0,p)
x4 = mvrnorm(n,mux4, sigx4, tol = 1e-6, empirical = FALSE)



indx = list()
indx0 = list()
indx.seperate1 = indx.seperate2 = indx.seperate3 = indx.seperate4 = list()
beta1_est = matrix(rep(0,p*N),nrow = p)
beta2_est = matrix(rep(0,p*N),nrow = p)
beta3_est = matrix(rep(0,p*N),nrow = p)
beta4_est = matrix(rep(0,p*N),nrow = p)
indx.fusion = list()

for(iter in 1:N){
  # iter = 1
  #generate y1, y2, y3 and y4
  
  sig_cor = matrix(c(sig1,corr,corr,corr,corr,sig2,corr,corr,corr,corr,sig3,corr,corr,corr,corr,sig4),nrow = 4)
  y = mvrnorm(n, mu = c(0,0,0,0), Sigma=sig_cor)+mvrnorm(n, mu = c(0,0,0,0), Sigma=sig_cor)^2+ cbind(x1%*%beta1_0,x2%*%beta2_0,x3%*%beta3_0,x4%*%beta4_0)
  y1 = y[,1]
  y2 = y[,2]
  y3 = y[,3]
  y4 = y[,4]
  
  
  y.list = list()
  y.list[[1]] = y1
  y.list[[2]] = y2
  y.list[[3]] = y3
  y.list[[4]] = y4
  x.list = list()
  x.list[[1]] = x1
  x.list[[2]] = x2
  x.list[[3]] = x3
  x.list[[4]] = x4
  
  
  
  #using ranked t to search starting point
  t.value = numeric() 
  for(i in 1:p){
    t.value = c(t.value , summary(lm(y1~x1[,i]-1))$coefficients[,3])
  }
  t.rank = sort(t.value,decreasing = T)
  starting.point = numeric()
  for (i in 1:an) {
    starting.point = c(starting.point, which(t.value==t.rank[i]))
  }
  
  # starting.point
  label = rep(0,p)
  label[starting.point]=1
  # p=1000
  # sn=10
  p = p1
  # label = c(c(1,1,1,1,1,1,1,1,0,1),rep(0,(p-sn)))
  
  gamma1 = which(label == 1)         # start from nested model
  label1 = label2 = label3 = label4 = label
  gamma.seperate1 = which(label1==1)
  gamma.seperate2 = which(label2==1)
  gamma.seperate3 = which(label3==1)
  gamma.seperate4 = which(label4==1)
  l = length(gamma1)
  l1 = length(gamma.seperate1)
  l2 = length(gamma.seperate2)
  l3 = length(gamma.seperate3)
  l4 = length(gamma.seperate4)
  if(l==0){
    poster = ((n+v1)/2)*log(2/(1+t(y1)%*%y1))+((n+v2)/2)*log(2/(1+t(y2)%*%y2))+((n+v3)/2)*log(2/(1+t(y3)%*%y3))+((n+v4)/2)*log(2/(1+t(y4)%*%y4))
    poster1 = ((n+v1)/2)*log(2/(1+t(y1)%*%y1))
    poster2 = ((n+v2)/2)*log(2/(1+t(y2)%*%y2))
    poster3 = ((n+v3)/2)*log(2/(1+t(y3)%*%y3))
    poster4 = ((n+v4)/2)*log(2/(1+t(y4)%*%y4))
  }else{
    x1_g = as.matrix(x1[,gamma1])
    Sigma1_g = diag(diag(Sigma1)[gamma1],nrow = l)
    x2_g = as.matrix(x2[,gamma1])
    Sigma2_g = diag(diag(Sigma2)[gamma1],nrow = l)
    x3_g = as.matrix(x3[,gamma1])
    Sigma3_g = diag(diag(Sigma3)[gamma1],nrow = l)
    x4_g = as.matrix(x4[,gamma1])
    Sigma4_g = diag(diag(Sigma4)[gamma1],nrow = l)
    u1_g = pd.solve(Sigma1_g)+ t(x1_g)%*%x1_g
    w1_g = Sigma1_g%*%u1_g
    u2_g = pd.solve(Sigma2_g)+ t(x2_g)%*%x2_g
    w2_g = Sigma2_g%*%u2_g
    u3_g = pd.solve(Sigma3_g)+ t(x3_g)%*%x3_g
    w3_g = Sigma3_g%*%u3_g
    u4_g = pd.solve(Sigma4_g)+ t(x4_g)%*%x4_g
    w4_g = Sigma4_g%*%u4_g
    
    poster1 = -1/2*sum(log(eigen(w1_g)$value))+((n+v1)/2)*log(2/(1+t(y1)%*%(diag(n)-x1_g%*%solve(u1_g)%*%t(x1_g))%*%y1))
    poster2 = -1/2*sum(log(eigen(w2_g)$value))+((n+v2)/2)*log(2/(1+t(y2)%*%(diag(n)-x2_g%*%solve(u2_g)%*%t(x2_g))%*%y2))
    poster3 = -1/2*sum(log(eigen(w3_g)$value))+((n+v3)/2)*log(2/(1+t(y3)%*%(diag(n)-x3_g%*%solve(u3_g)%*%t(x3_g))%*%y3))
    poster4 = -1/2*sum(log(eigen(w4_g)$value))+((n+v4)/2)*log(2/(1+t(y4)%*%(diag(n)-x4_g%*%solve(u4_g)%*%t(x4_g))%*%y4))
    poster = poster1 + poster2 + poster3 + poster4
    
  }
  indx[[iter]] = rep(0,p)
  indx0[[iter]] = rep(0,p)
  indx[[iter]][gamma1] = 1
  indx0[[iter]][gamma1] = 1
  indx.seperate1[[iter]] = indx.seperate2[[iter]] = indx.seperate3[[iter]] = indx.seperate4[[iter]] = rep(0,p)
  indx.seperate1[[iter]][gamma1] = indx.seperate2[[iter]][gamma1] = indx.seperate3[[iter]][gamma1] = indx.seperate4[[iter]][gamma1] = 1
  for(m in 1:M){
    for(j in 1:p){
      label1_up = label1
      label1_up[j] = 1-label1[j]
      gamma1_up = which(label1_up == 1)
      l1 = length(gamma1_up)
      
      if(l1==0){
        poster1_new = ((n+v1)/2)*log(2/(1+t(y1)%*%y1))
      }else{
        x1_g = as.matrix(x1[,gamma1_up])
        Sigma1_g = diag(diag(Sigma1)[gamma1_up],nrow = l1)
        
        u1_g = solve(Sigma1_g)+ t(x1_g)%*%x1_g
        w1_g = Sigma1_g%*%u1_g
        
        
        poster1_new = -1/2*sum(log(eigen(w1_g)$value))+((n+v1)/2)*log(2/(1+t(y1)%*%(diag(n)-x1_g%*%solve(u1_g)%*%t(x1_g))%*%y1))        
        
      }
      if(label1_up[j]==1){
        success = 1/(exp((poster1-poster1_new))+1)
      }else{
        success = 1/(exp((poster1_new-poster1))+1)
      }
      label1.temp = label1
      label1.temp[j] = rbinom(1,1,success)     
      if(sum(label1.temp)<=an){label1 = label1.temp}
      # label1[j] = rbinom(1,1,success)      # adjust by the square of success
      gamma.seperate1 = which(label1==1)
      if(label1[j]==label1_up[j]){
        poster1 = poster1_new
      }
      
      label2_up = label2
      label2_up[j] = 1-label2[j]
      gamma2_up = which(label2_up == 1)
      l2 = length(gamma2_up)
      
      if(l2==0){
        poster2_new = ((n+v2)/2)*log(2/(1+t(y2)%*%y2))
      }else{
        x2_g = as.matrix(x2[,gamma2_up])
        Sigma2_g = diag(diag(Sigma2)[gamma2_up],nrow = l2)
        
        u2_g = solve(Sigma2_g)+ t(x2_g)%*%x2_g
        w2_g = Sigma2_g%*%u2_g
        
        
        poster2_new = -1/2*sum(log(eigen(w2_g)$value))+((n+v2)/2)*log(2/(1+t(y2)%*%(diag(n)-x2_g%*%solve(u2_g)%*%t(x2_g))%*%y2))        
        
      }
      if(label2_up[j]==1){
        success = 1/(exp((poster2-poster2_new))+1)
      }else{
        success = 1/(exp((poster2_new-poster2))+1)
      }
      label2.temp = label2
      label2.temp[j] = rbinom(1,1,success)     
      if(sum(label2.temp)<=an){label2 = label2.temp}
      # label2[j] = rbinom(1,1,success)      # adjust by the square of success
      gamma.seperate2 = which(label2==1)
      if(label2[j]==label2_up[j]){
        poster2 = poster2_new
      }
      
      
      label3_up = label3
      label3_up[j] = 1-label3[j]
      gamma3_up = which(label3_up == 1)
      l3 = length(gamma3_up)
      
      if(l3==0){
        poster3_new = ((n+v3)/2)*log(2/(1+t(y3)%*%y3))
      }else{
        x3_g = as.matrix(x3[,gamma3_up])
        Sigma3_g = diag(diag(Sigma3)[gamma3_up],nrow = l3)
        
        u3_g = solve(Sigma3_g)+ t(x3_g)%*%x3_g
        w3_g = Sigma3_g%*%u3_g
        
        
        poster3_new = -1/2*sum(log(eigen(w3_g)$value))+((n+v3)/2)*log(2/(1+t(y3)%*%(diag(n)-x3_g%*%solve(u3_g)%*%t(x3_g))%*%y3))        
        
      }
      if(label3_up[j]==1){
        success = 1/(exp((poster3-poster3_new))+1)
      }else{
        success = 1/(exp((poster3_new-poster3))+1)
      }
      label3.temp = label3
      label3.temp[j] = rbinom(1,1,success)     
      if(sum(label3.temp)<=an){label3 = label3.temp}
      # label3[j] = rbinom(1,1,success)      # adjust by the square of success
      gamma.seperate3 = which(label3==1)
      if(label3[j]==label3_up[j]){
        poster3 = poster3_new
      }
      
      label4_up = label4
      label4_up[j] = 1-label4[j]
      gamma4_up = which(label4_up == 1)
      l4 = length(gamma4_up)
      
      if(l4==0){
        poster4_new = ((n+v4)/2)*log(2/(1+t(y4)%*%y4))
      }else{
        x4_g = as.matrix(x4[,gamma4_up])
        Sigma4_g = diag(diag(Sigma4)[gamma4_up],nrow = l4)
        
        u4_g = solve(Sigma4_g)+ t(x4_g)%*%x4_g
        w4_g = Sigma4_g%*%u4_g
        
        
        poster4_new = -1/2*sum(log(eigen(w4_g)$value))+((n+v4)/2)*log(2/(1+t(y4)%*%(diag(n)-x4_g%*%solve(u4_g)%*%t(x4_g))%*%y4))        
        
      }
      if(label4_up[j]==1){
        success = 1/(exp((poster4-poster4_new))+1)
      }else{
        success = 1/(exp((poster4_new-poster4))+1)
      }
      label4.temp = label4
      label4.temp[j] = rbinom(1,1,success)     
      if(sum(label4.temp)<=an){label4 = label4.temp}
      # label4[j] = rbinom(1,1,success)      # adjust by the square of success
      gamma.seperate4 = which(label4==1)
      if(label4[j]==label4_up[j]){
        poster4 = poster4_new
      }
      
      label_up = label
      label_up[j] = 1-label[j]
      gamma_up = which(label_up == 1)
      l = length(gamma_up)
      if(l==0){
        poster_new = ((n+v1)/2)*log(2/(1+t(y1)%*%y1))+((n+v2)/2)*log(2/(1+t(y2)%*%y2))+ ((n+v3)/2)*log(2/(1+t(y3)%*%y3))+((n+v4)/2)*log(2/(1+t(y4)%*%y4))
      }else{
        x1_g = as.matrix(x1[,gamma_up])
        Sigma1_g = diag(diag(Sigma1)[gamma_up],nrow = l)
        x2_g = as.matrix(x2[,gamma_up])
        Sigma2_g = diag(diag(Sigma2)[gamma_up],nrow = l)
        x3_g = as.matrix(x3[,gamma_up])
        Sigma3_g = diag(diag(Sigma3)[gamma_up],nrow = l)
        x4_g = as.matrix(x4[,gamma_up])
        Sigma4_g = diag(diag(Sigma4)[gamma_up],nrow = l)
        u1_g = pd.solve(Sigma1_g)+ t(x1_g)%*%x1_g
        w1_g = Sigma1_g%*%u1_g
        u2_g = pd.solve(Sigma2_g)+ t(x2_g)%*%x2_g
        w2_g = Sigma2_g%*%u2_g
        u3_g = pd.solve(Sigma3_g)+ t(x3_g)%*%x3_g
        w3_g = Sigma3_g%*%u3_g
        u4_g = pd.solve(Sigma4_g)+ t(x4_g)%*%x4_g
        w4_g = Sigma4_g%*%u4_g
        # poster_new = 1/sqrt(det(w1_g))*(2/(1+t(y1)%*%(diag(n)-x1_g%*%pd.solve(u1_g)%*%t(x1_g))%*%y1))^((n+v1)/2)*1/sqrt(det(w2_g))*(2/(1+t(y2)%*%(diag(n)-x2_g%*%pd.solve(u2_g)%*%t(x2_g))%*%y2))^((n+v2)/2)
        # poster_new = -1/2*log(det(w1_g))+((n+v1)/2)*log(2/(1+t(y1)%*%(diag(n)-x1_g%*%pd.solve(u1_g)%*%t(x1_g))%*%y1))-1/2*log(det(w2_g))+((n+v2)/2)*log(2/(1+t(y2)%*%(diag(n)-x2_g%*%pd.solve(u2_g)%*%t(x2_g))%*%y2))
        poster_new = -1/2*sum(log(eigen(w1_g)$value))+((n+v1)/2)*log(2/(1+t(y1)%*%(diag(n)-x1_g%*%pd.solve(u1_g)%*%t(x1_g))%*%y1))-1/2*sum(log(eigen(w2_g)$value))+((n+v2)/2)*log(2/(1+t(y2)%*%(diag(n)-x2_g%*%pd.solve(u2_g)%*%t(x2_g))%*%y2))
        +-1/2*sum(log(eigen(w3_g)$value))+((n+v3)/2)*log(2/(1+t(y3)%*%(diag(n)-x3_g%*%pd.solve(u3_g)%*%t(x3_g))%*%y3))-1/2*sum(log(eigen(w4_g)$value))+((n+v4)/2)*log(2/(1+t(y4)%*%(diag(n)-x4_g%*%pd.solve(u4_g)%*%t(x4_g))%*%y4))
        
        
      }
      if(label_up[j]==1){
        success = 1/(exp((poster-poster_new))+1)
      }else{
        success = 1/(exp((poster_new-poster))+1)
      }
      label.temp = label
      label.temp[j] = rbinom(1,1,success)     
      if(sum(label.temp)<=an){label = label.temp}
      # label[j] = rbinom(1,1,success)      # adjust by the square of success
      gamma1 = which(label==1)
      if(label[j]==label_up[j]){
        poster = poster_new
      }
      
      
    }
    ma1 = rep(0,p)
    ma1[gamma.seperate1] = 1
    indx.seperate1[[iter]] =cbind(indx.seperate1[[iter]],ma1)
    ma2 = rep(0,p)
    ma2[gamma.seperate2] = 1
    indx.seperate2[[iter]] =cbind(indx.seperate2[[iter]],ma2)
    ma3 = rep(0,p)
    ma3[gamma.seperate3] = 1
    indx.seperate3[[iter]] =cbind(indx.seperate3[[iter]],ma3)
    ma4 = rep(0,p)
    ma4[gamma.seperate4] = 1
    indx.seperate4[[iter]] =cbind(indx.seperate4[[iter]],ma4)
    ma = rep(0,p)
    ma[gamma1] = 1
    indx[[iter]] = cbind(indx[[iter]],ma)
    
  }
  
  ################################
  #using fusionlearn package
  ################################
  
  model1 = fusionbase.fit(x.list,y.list,seq(0,5,length.out = 20),n,p,4,depen = 'CORR')
  lambda <- model1[which.min(model1[,2]),1]
  result <- fusionbase(x.list,y.list,lambda,n,p,4)
  indx.fusion[[iter]] = result
  
  # plot(beta1_0,result[[iter]][,1])
  # plot(beta2_0,result[[iter]][,2])
  # plot(beta3_0,result[[iter]][,3])
  # plot(beta4_0,result[[iter]][,4])
  
  # l_fusion =  (rowSums(result$beta==0)==0)
  # res.fusion = res.fusion + (sum((l_fusion-true.model)^2)==0)
  # res_negative.fusion = res_negative.fusion + (sum((l_fusion-true.model)[rd]^2)>=0)
  # res_positive.fusion = res_positive.fusion + (sum((l_fusion-true.model)[-rd]^2)<=0)
  
  
  
  
  # dt0 = data.table(t(indx0[[iter]][,-(1:50)*p]))
  # dt = data.table(t(indx[[iter]][,-(1:50)]))
  # head(dt[, list(repeats=.N, id=.I[[1]]), by=names(dt)][order(repeats, decreasing=T)], 1)
  # l_est1 = as.numeric(dt[, list(repeats=.N, id=.I[[1]]), by=names(dt)][order(repeats, decreasing=T)][1,])[1:p]
  # l_est0 = as.numeric(dt0[, list(repeats=.N, id=.I[[1]]), by=names(dt)][order(repeats, decreasing=T)][1,])[1:p]
  # deci = rowMeans(indx[[iter]][,(M/2+0.5):(M)])
  # l_est= (deci>=cut_off)
  # if(l_est != 0){
  #   x1_g = as.matrix(x1[,gamma_est])
  #   Sigma1_g = diag(diag(Sigma1)[gamma_est],nrow = l_est)
  #   x2_g = as.matrix(x2[,gamma_est])
  #   Sigma2_g = diag(diag(Sigma2)[gamma_est],nrow = l_est)
  #   u1_g = pd.solve(Sigma1_g)+ t(x1_g)%*%x1_g
  #   w1_g = Sigma1_g%*%u1_g
  #   u2_g = pd.solve(Sigma2_g)+ t(x2_g)%*%x2_g
  #   w2_g = Sigma2_g%*%u2_g
  #   beta1_gest = pd.solve(u1_g)%*%t(x1_g)%*%y1
  #   beta2_gest = pd.solve(u2_g)%*%t(x2_g)%*%y2
  #   beta1_est[gamma_est,iter] = beta1_gest
  #   beta2_est[gamma_est,iter] = beta2_gest
  # }
  # res = res + (sum((l_est-true.model)^2)==0)
  # res_positive = res_positive + (sum((l_est-true.model)[rd]^2)>=0)
  # res_negative = res_negative + (sum((l_est-true.model)[-rd]^2)>=0)
  # res1 = res1 + (sum((l_est1-c(rep(1,sn),rep(0,p-sn)))^2)==0)
  # res1_positive = res1_positive + (sum((l_est1-true.model)[rd]^2)>=0)
  # res1_negative = res1_negative + (sum((l_est1-true.model)[-rd]^2)>=0)
  # res0 = res0 + (sum((l_est0-true.model)^2)==0)
  # res0_positive = res0_positive + (sum((l_est0-true.model)[rd]^2)>=0)
  # res0_negative = res0_negative + (sum((l_est0-true.model)[-rd]^2)>=0)
}




total.time = proc.time()-start.time
# print(res/N)
# print(res.fusion/N)
# print(res1/N)
# print(res0/N)


ps = function(n,iter,indx){
  x = indx[[n]][,iter]
  resu = sum(x[rd]==1)
  return(resu)
}

# ps(1,1,indx)

psr = function(n,iter,indx){
  resu = ps(n,iter,indx)/sn
  return(resu)
}

# psr(1,1,indx)

fd = function(n,iter,indx){
  x = indx[[n]][,iter]
  resu = sum(x[-rd]!=0)
  return(resu)
}

# fd(1,1,indx)

fdr = function(n,iter,indx){
  resu = fd(n,iter,indx)/sum(indx[[n]][,iter])
  return(resu)
}

# fdr(1,1,indx)

indx.and = list()
indx.or = list()
for (i in 1:N) {
  indx.and[[i]] = as.matrix((indx.seperate1[[i]] + indx.seperate2[[i]] + indx.seperate3[[i]] + indx.seperate4[[i]])==4)
  indx.or[[i]] = as.matrix((indx.seperate1[[i]] + indx.seperate2[[i]] + indx.seperate3[[i]] + indx.seperate4[[i]])!=0)
}

#=================intersection model===================

psr.and.res = numeric()
psr.and.cum = 0
ps.and.list = matrix(nrow = M,ncol = N)
ps.and.cum.list = matrix(nrow = M,ncol = N)
fdr.and.res = numeric()
fdr.and.cum = 0
fd.and.length = numeric()
fd.and.list = matrix(nrow = M,ncol = N)
fd.and.cum.list = matrix(nrow = M,ncol = N)

for(i in 1:N){
  for(j in 1:M){
    psr.and.unit = psr(i,j,indx.and)
    psr.and.cum = psr.and.cum + psr.and.unit*sn
    psr.and.res = c(psr.and.res,psr.and.unit)
    ps.and.list[j,i] = psr.and.unit
    ps.and.cum.list[j,i] = psr.and.cum/(length(psr.and.res)*sn+sn)
    fdr.and.unit = fdr(i,j,indx.and)
    fd.and.length = c(fd.and.length,sum(indx.and[[i]][,j]))
    fdr.and.cum = fdr.and.cum + fdr.and.unit*fd.and.length[(i-1)*M+j]
    fdr.and.res = c(fdr.and.res,fdr.and.unit)
    fd.and.list[j,i] = fdr.and.unit
    fd.and.cum.list[j,i] = fdr.and.cum/(sum(fd.and.length))
  }
}

# plot(rowMeans(fd.and.list),type = 'l')
# plot(rowMeans(ps.and.list),type = 'l')
# plot(rowMeans(fd.and.cum.list),type = 'l')
# plot(rowMeans(ps.and.cum.list),type = 'l')
# 
# acf(rowMeans(fd.and.list),type = 'correlation',lag.max = M)
# acf(rowMeans(ps.and.list),type = 'correlation',lag.max = M)
# 
# mean(ps.and.list[(M/2):M,],na.rm = T)
# mean(fd.and.list[(M/2):M,],na.rm = T)


#=====================union model===============


psr.or.res = numeric()
psr.or.cum = 0
ps.or.list = matrix(nrow = M,ncol = N)
ps.or.cum.list = matrix(nrow = M,ncol = N)
fdr.or.res = numeric()
fdr.or.cum = 0
fd.or.length = numeric()
fd.or.list = matrix(nrow = M,ncol = N)
fd.or.cum.list = matrix(nrow = M,ncol = N)

for(i in 1:N){
  for(j in 1:M){
    psr.or.unit = psr(i,j,indx.or)
    psr.or.cum = psr.or.cum + psr.or.unit*sn
    psr.or.res = c(psr.or.res,psr.or.unit)
    ps.or.list[j,i] = psr.or.unit
    ps.or.cum.list[j,i] = psr.or.cum/(length(psr.or.res)*sn+sn)
    fdr.or.unit = fdr(i,j,indx.or)
    fd.or.length = c(fd.or.length,sum(indx.or[[i]][,j]))
    fdr.or.cum = fdr.or.cum + fdr.or.unit*fd.or.length[(i-1)*M+j]
    fdr.or.res = c(fdr.or.res,fdr.or.unit)
    fd.or.list[j,i] = fdr.or.unit
    fd.or.cum.list[j,i] = fdr.or.cum/(sum(fd.or.length))
  }
}

# plot(rowMeans(fd.or.list),type = 'l')
# plot(rowMeans(ps.or.list),type = 'l')
# plot(rowMeans(fd.or.cum.list),type = 'l')
# plot(rowMeans(ps.or.cum.list),type = 'l')
# 
# acf(rowMeans(fd.or.list),type = 'correlation',lag.max = M)
# acf(rowMeans(ps.or.list),type = 'correlation',lag.max = M)
# 
# mean(ps.or.list[(M/2):M,],na.rm = T)
# mean(fd.or.list[(M/2):M,],na.rm = T)


################
#result of y1~x1
################

psr.res1 = numeric()
psr.cum1 = 0
ps.list1 = matrix(nrow = M,ncol = N)
ps.cum.list1 = matrix(nrow = M,ncol = N)
fdr.res1 = numeric()
fdr.cum1 = 0
fd.length1 = numeric()
fd.list1 = matrix(nrow = M,ncol = N)
fd.cum.list1 = matrix(nrow = M,ncol = N)

for(i in 1:N){
  for(j in 1:M){
    psr.unit1 = psr(i,j,indx.seperate1)
    psr.cum1 = psr.cum1 + psr.unit1*sn
    psr.res1 = c(psr.res1,psr.unit1)
    ps.list1[j,i] = psr.unit1
    ps.cum.list1[j,i] = psr.cum1/(length(psr.res1)*sn+sn)
    fdr.unit1 = fdr(i,j,indx.seperate1)
    fd.length1 = c(fd.length1,sum(indx.seperate1[[i]][,j]))
    fdr.cum1 = fdr.cum1 + fdr.unit1*fd.length1[(i-1)*M+j]
    fdr.res1 = c(fdr.res1,fdr.unit1)
    fd.list1[j,i] = fdr.unit1
    fd.cum.list1[j,i] = fdr.cum1/(sum(fd.length1))
  }
}

# plot(rowMeans(fd.list1),type = 'l')
# plot(rowMeans(ps.list1),type = 'l')
# plot(rowMeans(fd.cum.list1),type = 'l')
# plot(rowMeans(ps.cum.list1),type = 'l')
# 
# acf(rowMeans(fd.list1),type = 'correlation',lag.max = M)
# acf(rowMeans(ps.list1),type = 'correlation',lag.max = M)
# 
# mean(ps.list1[(M/2):M,],na.rm = T)
# mean(fd.list1[(M/2):M,],na.rm = T)
################
#result of y1~x2
################

psr.res2 = numeric()
psr.cum2 = 0
ps.list2 = matrix(nrow = M,ncol = N)
ps.cum.list2 = matrix(nrow = M,ncol = N)
fdr.res2 = numeric()
fdr.cum2 = 0
fd.length2 = numeric()
fd.list2 = matrix(nrow = M,ncol = N)
fd.cum.list2 = matrix(nrow = M,ncol = N)

for(i in 1:N){
  for(j in 1:M){
    psr.unit2 = psr(i,j,indx.seperate2)
    psr.cum2 = psr.cum2 + psr.unit2*sn
    psr.res2 = c(psr.res2,psr.unit2)
    ps.list2[j,i] = psr.unit2
    ps.cum.list2[j,i] = psr.cum2/(length(psr.res2)*sn+sn)
    fdr.unit2 = fdr(i,j,indx.seperate2)
    fd.length2 = c(fd.length2,sum(indx.seperate2[[i]][,j]))
    fdr.cum2 = fdr.cum2 + fdr.unit2*fd.length2[(i-1)*M+j]
    fdr.res2 = c(fdr.res2,fdr.unit2)
    fd.list2[j,i] = fdr.unit2
    fd.cum.list2[j,i] = fdr.cum2/(sum(fd.length2))
  }
}

# plot(rowMeans(fd.list2),type = 'l')
# plot(rowMeans(ps.list2),type = 'l')
# plot(rowMeans(fd.cum.list2),type = 'l')
# plot(rowMeans(ps.cum.list2),type = 'l')
# 
# acf(rowMeans(fd.list2),type = 'correlation',lag.max = M)
# acf(rowMeans(ps.list2),type = 'correlation',lag.max = M)
# 
# mean(ps.list2[(M/2):M,],na.rm = T)
# mean(fd.list2[(M/2):M,],na.rm = T)
################
#result of y1~x3
################

psr.res3 = numeric()
psr.cum3 = 0
ps.list3 = matrix(nrow = M,ncol = N)
ps.cum.list3 = matrix(nrow = M,ncol = N)
fdr.res3 = numeric()
fdr.cum3 = 0
fd.length3 = numeric()
fd.list3 = matrix(nrow = M,ncol = N)
fd.cum.list3 = matrix(nrow = M,ncol = N)

for(i in 1:N){
  for(j in 1:M){
    psr.unit3 = psr(i,j,indx.seperate3)
    psr.cum3 = psr.cum3 + psr.unit3*sn
    psr.res3 = c(psr.res3,psr.unit3)
    ps.list3[j,i] = psr.unit3
    ps.cum.list3[j,i] = psr.cum3/(length(psr.res3)*sn+sn)
    fdr.unit3 = fdr(i,j,indx.seperate3)
    fd.length3 = c(fd.length3,sum(indx.seperate3[[i]][,j]))
    fdr.cum3 = fdr.cum3 + fdr.unit3*fd.length3[(i-1)*M+j]
    fdr.res3 = c(fdr.res3,fdr.unit3)
    fd.list3[j,i] = fdr.unit3
    fd.cum.list3[j,i] = fdr.cum3/(sum(fd.length3))
  }
}

# plot(rowMeans(fd.list3),type = 'l')
# plot(rowMeans(ps.list3),type = 'l')
# plot(rowMeans(fd.cum.list3),type = 'l')
# plot(rowMeans(ps.cum.list3),type = 'l')
# 
# acf(rowMeans(fd.list3),type = 'correlation',lag.max = M)
# acf(rowMeans(ps.list3),type = 'correlation',lag.max = M)
# 
# mean(ps.list3[(M/2):M,],na.rm = T)
# mean(fd.list3[(M/2):M,],na.rm = T)
################
#result of y1~x4
################

psr.res4 = numeric()
psr.cum4 = 0
ps.list4 = matrix(nrow = M,ncol = N)
ps.cum.list4 = matrix(nrow = M,ncol = N)
fdr.res4 = numeric()
fdr.cum4 = 0
fd.length4 = numeric()
fd.list4 = matrix(nrow = M,ncol = N)
fd.cum.list4 = matrix(nrow = M,ncol = N)

for(i in 1:N){
  for(j in 1:M){
    psr.unit4 = psr(i,j,indx.seperate4)
    psr.cum4 = psr.cum4 + psr.unit4*sn
    psr.res4 = c(psr.res4,psr.unit4)
    ps.list4[j,i] = psr.unit4
    ps.cum.list4[j,i] = psr.cum4/(length(psr.res4)*sn+sn)
    fdr.unit4 = fdr(i,j,indx.seperate4)
    fd.length4 = c(fd.length4,sum(indx.seperate4[[i]][,j]))
    fdr.cum4 = fdr.cum4 + fdr.unit4*fd.length4[(i-1)*M+j]
    fdr.res4 = c(fdr.res4,fdr.unit4)
    fd.list4[j,i] = fdr.unit4
    fd.cum.list4[j,i] = fdr.cum4/(sum(fd.length4))
  }
}

# plot(rowMeans(fd.list4),type = 'l')
# plot(rowMeans(ps.list4),type = 'l')
# plot(rowMeans(fd.cum.list4),type = 'l')
# plot(rowMeans(ps.cum.list4),type = 'l')
# 
# acf(rowMeans(fd.list4),type = 'correlation',lag.max = M)
# acf(rowMeans(ps.list4),type = 'correlation',lag.max = M)
# 
# mean(ps.list4[(M/2):M,],na.rm = T)
# mean(fd.list4[(M/2):M,],na.rm = T)


psr.res = numeric()
psr.cum = 0
psr.cum.res = numeric()
ps.list = matrix(nrow = M,ncol = N)
ps.cum.list = matrix(nrow = M,ncol = N)
fdr.res = numeric()
fdr.cum = 0
fdr.cum.res = numeric()
fd.length = numeric()
fd.list = matrix(nrow = M,ncol = N)
fd.cum.list = matrix(nrow = M,ncol = N)


for(i in 1:N){
  for(j in 5:M){
    psr.unit = psr(i,j,indx)
    psr.cum = psr.cum + psr.unit*sn
    psr.cum.res = c(psr.cum.res,psr.cum/(length(psr.res)*sn+sn))
    psr.res = c(psr.res,psr.unit)
    ps.list[j,i] = psr.unit
    ps.cum.list[j,i] = psr.cum/(length(psr.res)*sn+sn)
    fdr.unit = fdr(i,j,indx)
    fd.length = c(fd.length,sum(indx[[i]][,j]))
    fdr.cum = fdr.cum + fdr.unit*fd.length[(i-1)*M+j]
    fdr.cum.res = c(fdr.cum.res,fdr.cum/(sum(fd.length)))
    fdr.res = c(fdr.res,fdr.unit)
    fd.list[j,i] = fdr.unit
    fd.cum.list[j,i] = fdr.cum/(sum(fd.length))
  }
}

# plot(rowMeans(fd.list),type = 'l')
# plot(rowMeans(ps.list),type = 'l')
# plot(rowMeans(fd.cum.list),type = 'l')
# plot(rowMeans(ps.cum.list),type = 'l')
# 
# acf(rowMeans(fd.list),type = 'correlation',lag.max = M)
# acf(rowMeans(ps.list),type = 'correlation',lag.max = M)

# mean(ps.list[(M/2):M,],na.rm = T)
# mean(fd.list[(M/2):M,],na.rm = T)

psr.fusion = rep(0,N)
fdr.fusion = rep(0,N)
for(i in 1:N){
  x =  indx.fusion[[i]]$beta[,1]
  x = (x!=0)
  psr.fusion[i] = sum(x[rd]==1)/sn
  fdr.fusion[i] = sum(x[-rd]!=0)/sum(x)
}
# plot(psr.fusion,type = 'l')
# plot(fdr.fusion,type = 'l')

# mean(psr.fusion,na.rm = T)
# mean(fdr.fusion,na.rm = T)

summary.table = array(c('','bayesian','group lasso','intersection','union','marginal1','marginal2','marginal3','marginal4','psr',mean(ps.list[(M/2):M,],na.rm = T),mean(psr.fusion,na.rm = T),mean(ps.and.list[(M/2):M,],na.rm = T),mean(ps.or.list[(M/2):M,],na.rm = T),mean(ps.list1[(M/2):M,],na.rm = T),mean(ps.list2[(M/2):M,],na.rm = T),mean(ps.list3[(M/2):M,],na.rm = T),mean(ps.list4[(M/2):M,],na.rm = T),'fdr',mean(fd.list[(M/2):M,],na.rm = T),mean(fdr.fusion,na.rm = T),mean(fd.and.list[(M/2):M,],na.rm = T),mean(fd.or.list[(M/2):M,],na.rm = T),mean(fd.list1[(M/2):M,],na.rm = T),mean(fd.list2[(M/2):M,],na.rm = T),mean(fd.list3[(M/2):M,],na.rm = T),mean(fd.list4[(M/2):M,],na.rm = T)),dim = c(9,3))

