library(readr)


stockdata = read.csv("finance.csv",header = T) # read data from file, 

##dim of data 1 is 700 52

###filling the missing data

library(zoo)
newstock<-array(dim=c(700,51))
for (i in 2:52){
  tem<-na.approx(stockdata[,i],na.rm = FALSE)
  newstock[,(i-1)]<-tem
  print(i)
}

cleanstock<-newstock[1:698,-c(26,45)]  ###dim is 698, 49

names<-colnames(stockdata[,-c(27,46)])[-1]



today<-cleanstock[2:698,]
yest<-cleanstock[1:697,]
return<-log(yest/today)*100

acf(return[,1])
acf(return[,2])
acf(return[,3])
# Autocorrelation at All lags are insignificant at 5\% level for three responses.


oldn<-696
n<-232
s<-3
p<-46 # The last 46 columns are covariates. the first three are the three outcomes.
oldx<-list()
oldx[[1]] = oldx[[2]] = oldx[[3]] = return[,4:49]

###     

####standardization of x 

for (j in 1:p){
  for (i in 1:s){
    oldx[[i]][,j]<-oldx[[i]][,j]/sqrt(sum(oldx[[i]][,j]^2))*sqrt(n)
  }
}
###
subse<-seq(1,696,3)
x = list()
x[[1]]<-oldx[[1]][subse,]
x[[2]]<-oldx[[2]][subse,]
x[[3]]<-oldx[[3]][subse,]

oldy<-list()
oldy[[1]]<-return[1:696,1]
oldy[[2]]<-return[1:696,2]
oldy[[3]]<-return[1:696,3]

y = list()
y[[1]]<-oldy[[1]][subse]
y[[2]]<-oldy[[2]][subse]
y[[3]]<-oldy[[3]][subse]



xdat1<-data.frame(cbind(y[[1]],x[[1]]))
xdat2<-data.frame(cbind(y[[2]],x[[2]]))
xdat3<-data.frame(cbind(y[[3]],x[[3]]))

colnames(xdat1)[1]<-"y1"
colnames(xdat1)[2:(p+1)]<-paste("X", 1:p, sep="")

colnames(xdat2)[1]<-"y2"
colnames(xdat2)[2:(p+1)]<-paste("X", 1:p, sep="")

colnames(xdat3)[1]<-"y3"
colnames(xdat3)[2:(p+1)]<-paste("X", 1:p, sep="")


testx = list()
testx[[1]] = oldx[[1]][subse+2,]
testx[[2]] = oldx[[2]][subse+2,]
testx[[3]] = oldx[[3]][subse+2,]

testy = list()
testy[[1]] = oldy[[1]][subse+2]
testy[[2]] = oldy[[2]][subse+2]
testy[[3]] = oldy[[3]][subse+2]

########################
#bayesian method
########################
library(data.table)


an = 20
M = 10
Sigma1 = diag(4,p)                   #prior of c1 
Sigma2 = diag(4,p)                  #prior of c2
Sigma3 = diag(4,p)                   #prior of c3 

x1 = x[[1]]
x2 = x[[2]]
x3 = x[[3]]

y1 = y[[1]]
y2 = y[[2]]
y3 = y[[3]]


x1.test = testx[[1]]
x2.test = testx[[2]]
x3.test = testx[[3]]

y1.test = testy[[1]]
y2.test = testy[[2]]
y3.test = testy[[3]]

p = length(x1[1,])
n = length(x1[,1])
v1 = v2 = v3 = 2

indx = rep(0,p)
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



gamma1 = which(label == 1)         # start from nested model

l = length(gamma1)
if(l==0){
  poster = ((n+v1)/2)*log(2/(1+t(y1)%*%y1))+((n+v2)/2)*log(2/(1+t(y2)%*%y2))+((n+v3)/2)*log(2/(1+t(y3)%*%y3))
}else{
  x1_g = as.matrix(x1[,gamma1])
  Sigma1_g = diag(diag(Sigma1)[gamma1],nrow = l)
  x2_g = as.matrix(x2[,gamma1])
  Sigma2_g = diag(diag(Sigma2)[gamma1],nrow = l)
  x3_g = as.matrix(x3[,gamma1])
  Sigma3_g = diag(diag(Sigma3)[gamma1],nrow = l)
  
  u1_g = solve(Sigma1_g)+ t(x1_g)%*%x1_g
  w1_g = Sigma1_g%*%u1_g
  u2_g = solve(Sigma2_g)+ t(x2_g)%*%x2_g
  w2_g = Sigma2_g%*%u2_g
  u3_g = solve(Sigma3_g)+ t(x3_g)%*%x3_g
  w3_g = Sigma3_g%*%u3_g
  
  poster = -1/2*sum(log(eigen(w1_g)$value))+((n+v1)/2)*log(2/(1+t(y1)%*%(diag(n)-x1_g%*%solve(u1_g)%*%t(x1_g))%*%y1))-1/2*sum(log(eigen(w2_g)$value))+((n+v2)/2)*log(2/(1+t(y2)%*%(diag(n)-x2_g%*%solve(u2_g)%*%t(x2_g))%*%y2))+ -1/2*sum(log(eigen(w3_g)$value))+((n+v3)/2)*log(2/(1+t(y3)%*%(diag(n)-x3_g%*%solve(u3_g)%*%t(x3_g))%*%y3))
  
}

indx[gamma1] = 1

for(m in 1:M){
  for(j in 1:p){
    label_up = label
    label_up[j] = 1-label[j]
    gamma_up = which(label_up == 1)
    l = length(gamma_up)
    if(l==0){
      poster_new = ((n+v1)/2)*log(2/(1+t(y1)%*%y1))+((n+v2)/2)*log(2/(1+t(y2)%*%y2))+ ((n+v3)/2)*log(2/(1+t(y3)%*%y3))
    }else{
      x1_g = as.matrix(x1[,gamma_up])
      Sigma1_g = diag(diag(Sigma1)[gamma_up],nrow = l)
      x2_g = as.matrix(x2[,gamma_up])
      Sigma2_g = diag(diag(Sigma2)[gamma_up],nrow = l)
      x3_g = as.matrix(x3[,gamma_up])
      Sigma3_g = diag(diag(Sigma3)[gamma_up],nrow = l)
      
      u1_g = solve(Sigma1_g)+ t(x1_g)%*%x1_g
      w1_g = Sigma1_g%*%u1_g
      u2_g = solve(Sigma2_g)+ t(x2_g)%*%x2_g
      w2_g = Sigma2_g%*%u2_g
      u3_g = solve(Sigma3_g)+ t(x3_g)%*%x3_g
      w3_g = Sigma3_g%*%u3_g
      
      
      poster_new = -1/2*sum(log(eigen(w1_g)$value))+((n+v1)/2)*log(2/(1+t(y1)%*%(diag(n)-x1_g%*%solve(u1_g)%*%t(x1_g))%*%y1))-1/2*sum(log(eigen(w2_g)$value))+((n+v2)/2)*log(2/(1+t(y2)%*%(diag(n)-x2_g%*%solve(u2_g)%*%t(x2_g))%*%y2))+-1/2*sum(log(eigen(w3_g)$value))+((n+v3)/2)*log(2/(1+t(y3)%*%(diag(n)-x3_g%*%solve(u3_g)%*%t(x3_g))%*%y3))
      
      
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
    ma = rep(0,p)
    ma[gamma1] = 1
    indx =cbind(indx,ma)
  }
}
dt = data.table(t(indx[,-(1:50)]))

l_est = as.numeric(dt[, list(repeats=.N, id=.I[[1]]), by=names(dt)][order(repeats, decreasing=T)][1,])[1:p]

subname<-names[4:49]
subname[(1:46)*l_est]


model1 = glm(y1~x1[,which(l_est==1)],family=gaussian(link=identity))
summary(model1)
model2 = glm(y2~x2[,which(l_est==1)],family=gaussian(link=identity))
summary(model2)
model3 = glm(y3~x3[,which(l_est==1)],family=gaussian(link=identity))
summary(model3)
predy1<-cbind(rep(1,232),testx[[1]][,which(l_est==1)])%*%model1$coefficients
predy2<-cbind(rep(1,232),testx[[2]][,which(l_est==1)])%*%model2$coefficients
predy3<-cbind(rep(1,232),testx[[3]][,which(l_est==1)])%*%model3$coefficients
b1 = mean((testy[[1]]-predy1)^2)
b2 = mean((testy[[2]]-predy2)^2)
b3 = mean((testy[[3]]-predy3)^2)


# =================full model===============

modelfully1 = glm(y1~x1,family=gaussian(link=identity))
summary(modelfully1)
modelfully2 = glm(y2~x2,family=gaussian(link=identity))
summary(modelfully2)
modelfully3 = glm(y3~x3,family=gaussian(link=identity))
summary(modelfully3)
predfullyy1<-cbind(rep(1,232),testx[[1]])%*%modelfully1$coefficients
predfullyy2<-cbind(rep(1,232),testx[[2]])%*%modelfully2$coefficients
predfullyy3<-cbind(rep(1,232),testx[[3]])%*%modelfully3$coefficients
fm1 = mean((testy[[1]]-predfullyy1)^2)
fm2 = mean((testy[[2]]-predfullyy2)^2)
fm3 = mean((testy[[3]]-predfullyy3)^2)





#================Group lasso======================

model.lasso <- fusionbase.fit(x,y,seq(0,10,length.out = 100),232,46,3,depen="CORR")
lambda <- model.lasso[which.min(model.lasso[,2]),1]
result <- fusionbase(x,y,lambda,232,46,3)

id <-which(result$beta[,1]!=0)
colnames(stockindexVIX)[id+1] 

model.lm1 = lm(y[[1]]~x[[1]][,id])
summary(model.lm1)
model.lm2 = lm(y[[2]]~x[[2]][,id])
summary(model.lm2)
model.lm3 = lm(y[[3]]~x[[3]][,id])
summary(model.lm3)

predgly1<-cbind(rep(1,232),testx[[1]][,id])%*%model.lm1$coefficients
predgly2<-cbind(rep(1,232),testx[[2]][,id])%*%model.lm2$coefficients
predgly3<-cbind(rep(1,232),testx[[3]][,id])%*%model.lm3$coefficients
gl1 = mean((testy[[1]]-predgly1)^2)
gl2 =mean((testy[[2]]-predgly2)^2)
gl3 =mean((testy[[3]]-predgly3)^2)

result.summary = matrix(c(gl1,gl2,gl3,b1,b2,b3,fm1,fm2,fm3),nrow = 3,ncol = 3)
rownames(result.summary) = c('VIX','GSPC','DJI')
colnames(result.summary) = c('group lasso','Bayesian','full model')
print(result.summary)
print(c('group lasso',colnames(stockindexVIX)[id+1] ))
print(c('bayesian',subname[(1:46)*l_est]))
