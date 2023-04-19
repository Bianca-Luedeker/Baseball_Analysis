library(MCMCpack)
library(grDevices)
library(gtools)
library(diagram)


#####################################
#List of Functions

biasdirichlet<-function(a,n){
p<-length(a)
aplus<-sum(a)
one<-rep(1,p)

b <- n * trigamma(aplus)
Binv <- - diag(1/trigamma(a))/n
Hinv <- Binv - (Binv %*% one) %*% (t(one) %*%
Binv)/(1/b + c(t(one) %*% Binv %*% one))

H<-(-n)*diag(trigamma(a))+one%*%t(one)*(n*trigamma(aplus))

deriv.H<-matrix(as.vector(matrix(rep(n*psigamma(aplus,2),p^3),p,p*p)+t(apply(diag(-n*psigamma(a,2)),2,diag))),p*p,p)


bias<-.5*(-Hinv)%*%t(deriv.H)%*%as.vector(-Hinv)

return(bias)
}




dirichletmlebias<-function(X, ind, N)
{
# Aa <- MLE.se.Dirichlet.NR(X, ind, N)
# --------------------------------------------------------
# Aim: Finding the MLEs of the Dirichlet parameter
# vector a = (a_1, ..., a_n) and
# their standard errors
# Method: Using the Newton-Raphson algorithm (2.69)
# Input: X: an m x n observed data matrix
# ind=1: using (2.71) as the initial values
# ind=2: using (2.73) as the initial values
# N: the number of iterations required for
# the Newton-Raphson algorithm
# Output: A --- N by n matrix whose t-th row is a?(t)
# a --- the MLE of a = (a_1, ..., a_n)
# COV --- the estimated asymptotic covariance
# matrix of the MLE \hat{a}
# se --- the estimated standard errors of \hat{a}
# aCI90 - the 90% asymptotic CI of a
# --------------------------------------------------------
m <- dim(X)[1]
n <- dim(X)[2]
one <- rep(1, n)
xmean <- apply(X, 2, mean)
G <- (apply(X, 2, prod))^(1/m)
logG<-apply(log(X),2,sum)/m
if(ind == 1) {
a <- min(X) * one
}
if(ind == 2) {
ga <- - digamma(1)
b <- xmean
de0 <- sum(b * log(G))
aplus <- ((n - 1) * ga)/(sum(b * log(b)) - de0)
a <- b * aplus
}
A <- matrix(0, N, n)
for(tt in 1:N) {
aplus <- sum(a)
#g <- m * (digamma(aplus) - digamma(a) + log(G))
g <- m * (digamma(aplus) - digamma(a) + logG)
b <- m * trigamma(aplus)
Binv <- - diag(1/trigamma(a))/m
Hinv <- Binv - (Binv %*% one) %*% (t(one) %*%
Binv)/(1/b + c(t(one) %*% Binv %*% one))
a <- a - c(Hinv %*% g)
A[tt, ] <- a
}

bias<-biasdirichlet(a,m)
a<-as.vector(a-bias)


#Reminder m-sample size n-variables
#I<-solve(Hinv)
B<--m*diag(trigamma(a))
littleb<-matrix(rep(1,n*n),n,n)*m*digamma(sum(a))
Iinv<- solve(B) - (solve(B) %*% one) %*% (t(one) %*%
Binv)/(1/littleb + c(t(one) %*% solve(B) %*% one))


muhat<-a/sum(a)
shat<-sum(a)

Jmat<-rbind(cbind(matrix(shat*diag(n-1),n-1,n-1),muhat[-n]),c(rep(-shat,n-1),muhat[n]))
COVmus<--solve(Jmat)%*%Iinv%*%solve(t(Jmat))
contrast<-rbind(cbind(matrix(diag(n-1),n-1,n-1),rep(0,n-1)),c(rep(1,n-1),0))

COVmu<-contrast%*%COVmus%*%t(contrast)

return(list(muhat,COVmu))
}







loglike<-function(z,data,s){

N<-dim(data)[1]
k<-dim(data)[2]
logG<-apply(log(data),2,mean)
mu<-z/sum(z)

return(-N*sum( mu*s*logG-lgamma(s*mu)))
}


loglike2<-function(z,data1,s1,data2,s2){

N1<-dim(data1)[1]
N2<-dim(data2)[1]

logG1<-apply(log(data1),2,mean)
logG2<-apply(log(data2),2,mean)
mu<-z/sum(z)

loglike<-N1*sum( mu*s1*logG1-lgamma(s1*mu))+N2*sum( mu*s2*logG2-lgamma(s2*mu))
return(-loglike)
}




mleprecision<-function(data1,mu,N){
N1<-dim(data1)[1]
K<-dim(data1)[2]

logG<-apply(log(data1),2,mean)
s<-(K-1)/2/(-sum(mu*(logG-log(mu))))

for(i in 1:N){
num<-N1*digamma(s)-N1*sum(mu*(digamma(s*mu)))+N1*sum(mu*logG)
den<-N1*trigamma(s)-N1*sum(mu*mu*trigamma(s*mu))
invs<-(1/s)+((1/s)^2)*num/den
s<-1/invs
}

return(s)
}



dirichletnull<-function(d1,d2,N){
mu1<-apply(d1,2,mean)
mu2<-apply(d2,2,mean)
k<-dim(d1)[2]
s_1<-mleprecision(d1,mu1,20)
s_2<-mleprecision(d2,mu2,20)
start<-apply(rbind(d1,d2),2,mean)
x<-optim(par=start,loglike2,method="L-BFGS-B",lower=rep(0.01,k),data1=d1,s1=s_1,data2=d2,s2=s_2)
mu<-x$par/sum(x$par)


for(i in 1:N){

s_1<-mleprecision(d1,mu,20)
s_2<-mleprecision(d2,mu,20)
start<-mu*(s_1+s_2)
x<-optim(par=start,loglike2,method="L-BFGS-B",lower=rep(0.01,k),data1=d1,s1=s_1,data2=d2,s2=s_2)
mu<-x$par/sum(x$par)
}

return(list(mu,s_1,s_2))

}




dirichlet.test<-function(data1,data2){

k<-dim(data1)[2]
samplesize1<-dim(data1)[1]
samplesize2<-dim(data2)[1]

x<-dirichletnull(data1,data2,50)
muhat<-x[[1]]
s1hat<-x[[2]]
s2hat<-x[[3]]

alpha1<-muhat*s1hat
alpha2<-muhat*s2hat
one <- rep(1, k)

B<- -samplesize1*diag(trigamma(alpha1))
b<-matrix(rep(1,k*k),k,k)*samplesize1*digamma(sum(alpha1))
H=B+b
Hinv<- solve(B) - (solve(B) %*% one) %*% (t(one) %*%
solve(B))/(1/b + c(t(one) %*% solve(B) %*% one))



Jmat1<-rbind(cbind(matrix(s1hat*diag(k-1),k-1,k-1),muhat[-k]),c(rep(-s1hat,k-1),muhat[k]))
COVmus1<- -solve(Jmat1)%*%Hinv%*%solve(t(Jmat1))
contrast<-rbind(cbind(matrix(diag(k-1),k-1,k-1),rep(0,k-1)),c(rep(1,k-1),0))
COVmu1<-contrast%*%COVmus1%*%t(contrast)




B<- -samplesize2*diag(trigamma(alpha2))
b<-matrix(rep(1,k*k),k,k)*samplesize2*digamma(sum(alpha2))
H=B+b
Hinv<- solve(B) - (solve(B) %*% one) %*% (t(one) %*%
solve(B))/(1/b + c(t(one) %*% solve(B) %*% one))



Jmat2<-rbind(cbind(matrix(s2hat*diag(k-1),k-1,k-1),muhat[-k]),c(rep(-s2hat,k-1),muhat[k]))
COVmus2<- -solve(Jmat2)%*%Hinv%*%solve(t(Jmat2))
contrast2<-rbind(cbind(matrix(diag(k-1),k-1,k-1),rep(0,k-1)),c(rep(1,k-1),0))
COVmu2<-contrast2%*%COVmus2%*%t(contrast2)

finalCov<-COVmu1+COVmu2
dif<-dirichletmle(data1,ind=1,N=40)[[1]]-dirichletmle(data2,ind=1,N=40)[[1]]

return(list(dif,finalCov))

}





dirichletmle<-function(X, ind, N)
{
# Aa <- MLE.se.Dirichlet.NR(X, ind, N)
# --------------------------------------------------------
# Aim: Finding the MLEs of the Dirichlet parameter
# vector a = (a_1, ..., a_n) and
# their standard errors
# Method: Using the Newton-Raphson algorithm (2.69)
# Input: X: an m x n observed data matrix
# ind=1: using (2.71) as the initial values
# ind=2: using (2.73) as the initial values
# N: the number of iterations required for
# the Newton-Raphson algorithm
# Output: A --- N by n matrix whose t-th row is a?(t)
# a --- the MLE of a = (a_1, ..., a_n)
# COV --- the estimated asymptotic covariance
# matrix of the MLE \hat{a}
# se --- the estimated standard errors of \hat{a}
# aCI90 - the 90% asymptotic CI of a
# --------------------------------------------------------
m <- dim(X)[1]
n <- dim(X)[2]
one <- rep(1, n)
xmean <- apply(X, 2, mean)
G <- (apply(X, 2, prod))^(1/m)
logG<-apply(log(X),2,sum)/m
if(ind == 1) {
a <- min(X) * one
}
if(ind == 2) {
ga <- - digamma(1)
b <- xmean
de0 <- sum(b * log(G))
aplus <- ((n - 1) * ga)/(sum(b * log(b)) - de0)
a <- b * aplus
}
A <- matrix(0, N, n)
for(tt in 1:N) {
aplus <- sum(a)
g <- m * (digamma(aplus) - digamma(a) + logG)
b <- m * trigamma(aplus)
Binv <- - diag(1/trigamma(a))/m
Hinv <- Binv - (Binv %*% one) %*% (t(one) %*%
Binv)/(1/b + c(t(one) %*% Binv %*% one))
a <- a - c(Hinv %*% g)
A[tt, ] <- a
}
COV <- - Hinv
se <- sqrt(diag(COV))
aCI90 <- matrix(0, n, 2)
# 90% Confidence Interval of the parameter vector a
aCI90[, 1] <- a - qnorm(0.95) * se
aCI90[, 2] <- a + qnorm(0.95) * se

#Reminder m-sample size n-variables
#I<-solve(Hinv)
B<--m*diag(trigamma(a))
littleb<-matrix(rep(1,n*n),n,n)*m*digamma(sum(a))
Iinv<- solve(B) - (solve(B) %*% one) %*% (t(one) %*%
Binv)/(1/littleb + c(t(one) %*% solve(B) %*% one))


muhat<-a/sum(a)
shat<-sum(a)

Jmat<-rbind(cbind(matrix(shat*diag(n-1),n-1,n-1),muhat[-n]),c(rep(-shat,n-1),muhat[n]))
COVmus<--solve(Jmat)%*%Iinv%*%solve(t(Jmat))
contrast<-rbind(cbind(matrix(diag(n-1),n-1,n-1),rep(0,n-1)),c(rep(1,n-1),0))

COVmu<-contrast%*%COVmus%*%t(contrast)

return(list(muhat,COVmu))
}




dirichletmlealphas<-function(X, ind, N)
{
# Aa <- MLE.se.Dirichlet.NR(X, ind, N)
# --------------------------------------------------------
# Aim: Finding the MLEs of the Dirichlet parameter
# vector a = (a_1, ..., a_n) and
# their standard errors
# Method: Using the Newton-Raphson algorithm (2.69)
# Input: X: an m x n observed data matrix
# ind=1: using (2.71) as the initial values
# ind=2: using (2.73) as the initial values
# N: the number of iterations required for
# the Newton-Raphson algorithm
# Output: A --- N by n matrix whose t-th row is a?(t)
# a --- the MLE of a = (a_1, ..., a_n)
# COV --- the estimated asymptotic covariance
# matrix of the MLE \hat{a}
# se --- the estimated standard errors of \hat{a}
# aCI90 - the 90% asymptotic CI of a
# --------------------------------------------------------
m <- dim(X)[1]
n <- dim(X)[2]
one <- rep(1, n)
xmean <- apply(X, 2, mean)
G <- (apply(X, 2, prod))^(1/m)
logG<-apply(log(X),2,sum)/m
if(ind == 1) {
a <- min(X) * one
}
if(ind == 2) {
ga <- - digamma(1)
b <- xmean
de0 <- sum(b * log(G))
aplus <- ((n - 1) * ga)/(sum(b * log(b)) - de0)
a <- b * aplus
}
A <- matrix(0, N, n)
for(tt in 1:N) {
aplus <- sum(a)
#g <- m * (digamma(aplus) - digamma(a) + log(G))
g <- m * (digamma(aplus) - digamma(a) + logG)
b <- m * trigamma(aplus)
Binv <- - diag(1/trigamma(a))/m
Hinv <- Binv - (Binv %*% one) %*% (t(one) %*%
Binv)/(1/b + c(t(one) %*% Binv %*% one))
a <- a - c(Hinv %*% g)
A[tt, ] <- a
}
COV <- - Hinv
se <- sqrt(diag(COV))
aCI90 <- matrix(0, n, 2)
# 90% Confidence Interval of the parameter vector a
aCI90[, 1] <- a - qnorm(0.95) * se
aCI90[, 2] <- a + qnorm(0.95) * se

#Reminder m-sample size n-variables
#I<-solve(Hinv)
B<--m*diag(trigamma(a))
littleb<-matrix(rep(1,n*n),n,n)*m*trigamma(sum(a))
Iinv<- solve(B) - (solve(B) %*% one) %*% (t(one) %*%
Binv)/(1/littleb + c(t(one) %*% solve(B) %*% one))


muhat<-a/sum(a)
shat<-sum(a)

#Jmat<-rbind(cbind(matrix(shat*diag(n-1),n-1,n-1),muhat[-n]),c(rep(-shat,n-1),muhat[n]))
#COVmus<--solve(Jmat)%*%Iinv%*%solve(t(Jmat))
#contrast<-rbind(cbind(matrix(diag(n-1),n-1,n-1),rep(0,n-1)),c(rep(1,n-1),0))

#COVmu<-contrast%*%COVmus%*%t(contrast)

return(a)
}




ttest<-function(x,n1,n2){
return(t.test(x[1:n1],x[(n1+1):(n1+n2)])$p.value)}



wilcoxtest<-function(x,n1,n2){
return(wilcox.test(x[1:n1],x[(n1+1):(n1+n2)],paired=FALSE)$p.value)}



###########################################################
#Simulation of power



powersim1<-function(alpha1,size1,alpha2,size2,tot){

z<-c()
z2<-c()
z3<-c()
z4<-c()
z5<-c()

for (i in 1:tot){

sample1<-rdirichlet(size1,alpha1)
sample2<-rdirichlet(size2,alpha2)
x<-dirichletmlebias(sample1,ind=1,N=40)
y<-dirichletmlebias(sample2,ind=1,N=40)

z<-cbind(z,   ifelse(   (2-2*pnorm(abs((x[[1]]-y[[1]])/sqrt(diag(x[[2]]+y[[2]]))),0,1))<.05,1,0)     )
z2<-cbind(z2,ifelse(apply(rbind(sample1,sample2),2,ttest,n1=size1,n2=size2)<.05,1,0))


x1<-dirichletmle(sample1,ind=1,N=40)
y1<-dirichletmle(sample2,ind=1,N=40)

z3<-cbind(z3,   ifelse(  ( 2-2*pnorm(abs((x1[[1]]-y1[[1]])/sqrt(diag(x1[[2]]+y1[[2]]))),0,1))<.05,1,0)    )

x<-dirichlet.test(sample1,sample2)
z4<-cbind(z4,   ifelse(    (2-2*pnorm( abs(x[[1]])/sqrt(diag(x[[2]])),0,1))<.05,1,0))

z5<-cbind(z5,ifelse(apply(rbind(sample1,sample2),2,wilcoxtest,n1=size1,n2=size2)<.05,1,0))


}
power1<-apply(z,1,mean) #bias corrected
power2<-apply(z2,1,mean)#t test
power3<-apply(z3,1,mean)#normal dirichlet
power4<-apply(z4,1,mean)#mle conditioned on nuissance
power5<-apply(z5,1,mean)#rank sum test

return(list(power1,power2,power3,power4,power5))
}



cd28sim1<-powersim1(c(16.097,56.247),26,c(5.088,15.424),63,tot=10000)

cd28sim2<-powersim1(c(23.14,2.243,1.809,4.526),26,c(5.491,1.249,1.293,1.559),63,tot=5000)


cd28sim3<-powersim1(c(18.157,4.13,4.276,5.155),26,c(5.491,1.249,1.293,1.559),63,tot=5000)


cd28sim4<-powersim1(c(23.14,2.243,1.809,4.526),10,c(5.491,1.249,1.293,1.559),20,tot=5000)





x1<-powersim1(c(8.5,11.5),20,c(10,10),20,tot=5000)
x2<-powersim1(c(9,11),20,c(10,10),20,tot=5000)
x3<-powersim1(c(9.5,10.5),20,c(10,10),20,tot=5000)
x4<-powersim1(c(10,10),20,c(10,10),20,tot=5000)
x5<-powersim1(c(10.5,9.5),20,c(10,10),20,tot=5000)
x6<-powersim1(c(11,9),20,c(10,10),20,tot=5000)
x7<-powersim1(c(11.5,8.5),20,c(10,10),20,tot=5000)


x11<-powersim1(c(8.5,11.5),10,c(10,10),30,tot=5000)
x21<-powersim1(c(9,11),10,c(10,10),30,tot=5000)
x31<-powersim1(c(9.5,10.5),10,c(10,10),30,tot=5000)
x41<-powersim1(c(10,10),10,c(10,10),30,tot=5000)
x51<-powersim1(c(10.5,9.5),10,c(10,10),30,tot=5000)
x61<-powersim1(c(11,9),10,c(10,10),30,tot=5000)
x71<-powersim1(c(11.5,8.5),10,c(10,10),30,tot=5000)


x12<-powersim1(c(20,11.5,3.5,5),20,c(20,10,5,5)/2,20,tot=5000)
x22<-powersim1(c(20,11,4,5),20,c(20,10,5,5)/2,20,tot=5000)
x32<-powersim1(c(20,10.5,4.5,5),20,c(20,10,5,5)/2,20,tot=5000)
x42<-powersim1(c(20,10,5,5),20,c(20,10,5,5)/2,20,tot=5000)
x52<-powersim1(c(20,9.5,5.5,5),20,c(20,10,5,5)/2,20,tot=5000)
x62<-powersim1(c(20,9,6,5),20,c(20,10,5,5)/2,20,tot=5000)
x72<-powersim1(c(20,8.5,6.5,5),20,c(20,10,5,5)/2,20,tot=5000)



x13<-powersim1(c(20,11.5,3.5,5),20,c(20,10,5,5)/2,10,tot=5000)
x23<-powersim1(c(20,11,4,5),20,c(20,10,5,5)/2,10,tot=5000)
x33<-powersim1(c(20,10.5,4.5,5),20,c(20,10,5,5)/2,10,tot=5000)
x43<-powersim1(c(20,10,5,5),20,c(20,10,5,5)/2,10,tot=5000)
x53<-powersim1(c(20,9.5,5.5,5),20,c(20,10,5,5)/2,10,tot=5000)
x63<-powersim1(c(20,9,6,5),20,c(20,10,5,5)/2,10,tot=5000)
x73<-powersim1(c(20,8.5,6.5,5),20,c(20,10,5,5)/2,10,tot=5000)



x14<-powersim1(c(20,11.5,3.5,5),20,c(20,10,5,5),20,tot=5000)
x24<-powersim1(c(20,11,4,5),20,c(20,10,5,5),20,tot=5000)
x34<-powersim1(c(20,10.5,4.5,5),20,c(20,10,5,5),20,tot=5000)
x44<-powersim1(c(20,10,5,5),20,c(20,10,5,5),20,tot=5000)
x54<-powersim1(c(20,9.5,5.5,5),20,c(20,10,5,5),20,tot=5000)
x64<-powersim1(c(20,9,6,5),20,c(20,10,5,5),20,tot=5000)
x74<-powersim1(c(20,8.5,6.5,5),20,c(20,10,5,5),20,tot=5000)


x15<-powersim1(c(20,11.5,3.5,5),100,c(20,10,5,5)/2,100,tot=5000)
x25<-powersim1(c(20,11,4,5),100,c(20,10,5,5)/2,100,tot=5000)
x35<-powersim1(c(20,10.5,4.5,5),100,c(20,10,5,5)/2,100,tot=5000)
x45<-powersim1(c(20,10,5,5),100,c(20,10,5,5)/2,100,tot=5000)
x55<-powersim1(c(20,9.5,5.5,5),100,c(20,10,5,5)/2,100,tot=5000)
x65<-powersim1(c(20,9,6,5),100,c(20,10,5,5)/2,100,tot=5000)
x75<-powersim1(c(20,8.5,6.5,5),100,c(20,10,5,5)/2,100,tot=5000)


x16<-powersim1(c(1.25,8.75),20,c(2,8),20,tot=5000)
x26<-powersim1(c(1.5,8.5),20,c(2,8),20,tot=5000)
x36<-powersim1(c(1.75,8.25),20,c(2,8),20,tot=5000)
x46<-powersim1(c(2,8),20,c(2,8),20,tot=5000)
x56<-powersim1(c(2.25,7.75),20,c(2,8),20,tot=5000)
x66<-powersim1(c(2.5,7.5),20,c(2,8),20,tot=5000)
x76<-powersim1(c(2.75,7.25),20,c(2,8),20,tot=5000)


x17<-powersim1(c(1.25,8.75),20,c(2,8)*2,20,tot=5000)
x27<-powersim1(c(1.5,8.5),20,c(2,8)*2,20,tot=5000)
x37<-powersim1(c(1.75,8.25),20,c(2,8)*2,20,tot=5000)
x47<-powersim1(c(2,8),20,c(2,8)*2,20,tot=5000)
x57<-powersim1(c(2.25,7.75),20,c(2,8)*2,20,tot=5000)
x67<-powersim1(c(2.5,7.5),20,c(2,8)*2,20,tot=5000)
x77<-powersim1(c(2.75,7.25),20,c(2,8)*2,20,tot=5000)


x18<-powersim1(c(20,11.5,3.5,5),200,c(20,10,5,5)/2,100,tot=5000)
x28<-powersim1(c(20,11,4,5),200,c(20,10,5,5)/2,100,tot=5000)
x38<-powersim1(c(20,10.5,4.5,5),200,c(20,10,5,5)/2,100,tot=5000)
x48<-powersim1(c(20,10,5,5),200,c(20,10,5,5)/2,100,tot=5000)
x58<-powersim1(c(20,9.5,5.5,5),200,c(20,10,5,5)/2,100,tot=5000)
x68<-powersim1(c(20,9,6,5),200,c(20,10,5,5)/2,100,tot=5000)
x78<-powersim1(c(20,8.5,6.5,5),200,c(20,10,5,5)/2,100,tot=5000)


x19<-powersim1(c(20,11.5,3.5,5),10,c(20,10,5,5)/2,10,tot=5000)
x29<-powersim1(c(20,11,4,5),10,c(20,10,5,5)/2,10,tot=5000)
x39<-powersim1(c(20,10.5,4.5,5),10,c(20,10,5,5)/2,10,tot=5000)
x49<-powersim1(c(20,10,5,5),10,c(20,10,5,5)/2,10,tot=5000)
x59<-powersim1(c(20,9.5,5.5,5),10,c(20,10,5,5)/2,10,tot=5000)
x69<-powersim1(c(20,9,6,5),20,c(20,10,5,5)/2,10,tot=5000)
x79<-powersim1(c(20,8.5,6.5,5),10,c(20,10,5,5)/2,10,tot=5000)



#Gathering Results
xDbias<-rbind(x1[[1]],x2[[1]],x3[[1]],x4[[1]],x5[[1]],x6[[1]],x7[[1]])
xT<-rbind(x1[[2]],x2[[2]],x3[[2]],x4[[2]],x5[[2]],x6[[2]],x7[[2]])
xD<-rbind(x1[[3]],x2[[3]],x3[[3]],x4[[3]],x5[[3]],x6[[3]],x7[[3]])
xcondD<-rbind(x1[[4]],x2[[4]],x3[[4]],x4[[4]],x5[[4]],x6[[4]],x7[[4]])
xW<-rbind(x1[[5]],x2[[5]],x3[[5]],x4[[5]],x5[[5]],x6[[5]],x7[[5]])


dif<-seq(8.5,11.5,.5)/20-.5

plot(dif,xT[,1],type="l")
#points(dif,xD[,1])
lines(dif,xT[,1],col="red")
lines(dif,xD[,1],col="blue")
lines(dif,xcondD[,1],col="green")
lines(dif,xW[,1],col="black")
lines(dif,xDbias[,1],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")

xD1bias<-rbind(x11[[1]],x21[[1]],x31[[1]],x41[[1]],x51[[1]],x61[[1]],x71[[1]])
xT1<-rbind(x11[[2]],x21[[2]],x31[[2]],x41[[2]],x51[[2]],x61[[2]],x71[[2]])
xD1<-rbind(x11[[3]],x21[[3]],x31[[3]],x41[[3]],x51[[3]],x61[[3]],x71[[3]])
xcondD1<-rbind(x11[[4]],x21[[4]],x31[[4]],x41[[4]],x51[[4]],x61[[4]],x71[[4]])
xW1<-rbind(x11[[5]],x21[[5]],x31[[5]],x41[[5]],x51[[5]],x61[[5]],x71[[5]])

dif<-seq(8.5,11.5,.5)/20-.5

plot(dif,xT1[,1],type="l")
#points(dif,xD1[,1])
lines(dif,xT1[,1],col="red")
lines(dif,xD1[,1],col="blue")
lines(dif,xcondD1[,1],col="green")
lines(dif,xW1[,1],col="black")
lines(dif,xD1bias[,1],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")




xD2bias<-rbind(x12[[1]],x22[[1]],x32[[1]],x42[[1]],x52[[1]],x62[[1]],x72[[1]])
xT2<-rbind(x12[[2]],x22[[2]],x32[[2]],x42[[2]],x52[[2]],x62[[2]],x72[[2]])
xD2<-rbind(x12[[3]],x22[[3]],x32[[3]],x42[[3]],x52[[3]],x62[[3]],x72[[3]])
xcondD2<-rbind(x12[[4]],x22[[4]],x32[[4]],x42[[4]],x52[[4]],x62[[4]],x72[[4]])
xW2<-rbind(x12[[5]],x22[[5]],x32[[5]],x42[[5]],x52[[5]],x62[[5]],x72[[5]])



dif1<-rep(0,7)
dif2<-seq(8.5,11.5,.5)/40-.25
dif3<-seq(3.5,6.5,.5)/40-.125
dif4<-rep(0,7)


plot(dif2,xT2[7:1,2],type="l",ylim=c(.05,.40))
#points(dif2,xD2[7:1,2])
lines(dif2,xT2[7:1,2],col="red")
lines(dif2,xD2[7:1,2],col="blue")
lines(dif2,xcondD2[7:1,2],col="green")
lines(dif2,xW2[7:1,2],col="black")
lines(dif2,xD2bias[7:1,2],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")




plot(dif3,xD2[,3],type="l",ylim=c(.05,.7))
#points(dif3,xT2[,3])
lines(dif3,xT2[,3],col="red")
lines(dif3,xD2[,3],col="blue")
lines(dif3,xcondD2[,3],col="green")
lines(dif3,xW2[,3],col="black")
lines(dif3,xD2bias[,3],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")





xD3bias<-rbind(x13[[1]],x23[[1]],x33[[1]],x43[[1]],x53[[1]],x63[[1]],x73[[1]])
xT3<-rbind(x13[[2]],x23[[2]],x33[[2]],x43[[2]],x53[[2]],x63[[2]],x73[[2]])
xD3<-rbind(x13[[3]],x23[[3]],x33[[3]],x43[[3]],x53[[3]],x63[[3]],x73[[3]])
xcondD3<-rbind(x13[[4]],x23[[4]],x33[[4]],x43[[4]],x53[[4]],x63[[4]],x73[[4]])
xW3<-rbind(x13[[5]],x23[[5]],x33[[5]],x43[[5]],x53[[5]],x63[[5]],x73[[5]])

dif1<-rep(0,7)
dif2<-seq(8.5,11.5,.5)/40-.25
dif3<-seq(3.5,6.5,.5)/40-.125
dif4<-rep(0,7)


plot(dif2,xT3[7:1,2],type="l",ylim=c(.05,.275))
#points(dif2,xD3[7:1,2])
lines(dif2,xT3[7:1,2],col="red")
lines(dif2,xD3[7:1,2],col="blue")
lines(dif2,xcondD3[7:1,2],col="green")
lines(dif2,xW3[7:1,2],col="black")
lines(dif2,xD3bias[7:1,2],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")





plot(dif3,xD3[,3],type="l",ylim=c(0.05,.5))
#points(dif3,xT3[,3])
lines(dif3,xT3[,3],col="red")
lines(dif3,xD3[,3],col="blue")
lines(dif3,xcondD3[,3],col="green")
lines(dif3,xW3[,3],col="black")
lines(dif3,xD3bias[,3],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")







xD4bias<-rbind(x14[[1]],x24[[1]],x34[[1]],x44[[1]],x54[[1]],x64[[1]],x74[[1]])
xT4<-rbind(x14[[2]],x24[[2]],x34[[2]],x44[[2]],x54[[2]],x64[[2]],x74[[2]])
xD4<-rbind(x14[[3]],x24[[3]],x34[[3]],x44[[3]],x54[[3]],x64[[3]],x74[[3]])
xcondD4<-rbind(x14[[4]],x24[[4]],x34[[4]],x44[[4]],x54[[4]],x64[[4]],x74[[4]])
xW4<-rbind(x14[[5]],x24[[5]],x34[[5]],x44[[5]],x54[[5]],x64[[5]],x74[[5]])

dif1<-rep(0,7)
dif2<-seq(8.5,11.5,.5)/40-.25
dif3<-seq(3.5,6.5,.5)/40-.125
dif4<-rep(0,7)


plot(dif2,xD4[7:1,2],type="l")
#points(dif2,xT4[7:1,2])
lines(dif2,xT4[7:1,2],col="red")
lines(dif2,xD4[7:1,2],col="blue")
lines(dif2,xcondD4[7:1,2],col="green")
lines(dif2,xW4[7:1,2],col="black")
lines(dif2,xD4bias[7:1,2],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")



plot(dif3,xD4[,3],type="l")
#points(dif3,xT4[,3])
lines(dif3,xT4[,3],col="red")
lines(dif3,xD4[,3],col="blue")
lines(dif3,xcondD4[,3],col="green")
lines(dif3,xW4[,3],col="black")
lines(dif3,xD4bias[,3],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")








xD5bias<-rbind(x15[[1]],x25[[1]],x35[[1]],x45[[1]],x55[[1]],x65[[1]],x75[[1]])
xT5<-rbind(x15[[2]],x25[[2]],x35[[2]],x45[[2]],x55[[2]],x65[[2]],x75[[2]])
xD5<-rbind(x15[[3]],x25[[3]],x35[[3]],x45[[3]],x55[[3]],x65[[3]],x75[[3]])
xcondD5<-rbind(x15[[4]],x25[[4]],x35[[4]],x45[[4]],x55[[4]],x65[[4]],x75[[4]])
xW5<-rbind(x15[[5]],x25[[5]],x35[[5]],x45[[5]],x55[[5]],x65[[5]],x75[[5]])

dif1<-rep(0,7)
dif2<-seq(8.5,11.5,.5)/40-.25
dif3<-seq(3.5,6.5,.5)/40-.125
dif4<-rep(0,7)


plot(dif2,xD5[7:1,2],type="l",ylim=c(0,1))
#points(dif2,xT5[7:1,2])
lines(dif2,xT5[7:1,2],col="red")
lines(dif2,xD5[7:1,2],col="blue")
lines(dif2,xcondD5[7:1,2],col="green")
lines(dif2,xW5[7:1,2],col="black")
lines(dif2,xD5bias[7:1,2],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")



plot(dif3,xD5[,3],type="l",ylim=c(0,1))
#points(dif3,xT5[,3])
lines(dif3,xT5[,3],col="red")
lines(dif3,xD5[,3],col="blue")
lines(dif3,xcondD5[,3],col="green")
lines(dif3,xW5[,3],col="black")
lines(dif3,xD5bias[,3],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")






xD6bias<-rbind(x16[[1]],x26[[1]],x36[[1]],x46[[1]],x56[[1]],x66[[1]],x76[[1]])
xT6<-rbind(x16[[2]],x26[[2]],x36[[2]],x46[[2]],x56[[2]],x66[[2]],x76[[2]])
xD6<-rbind(x16[[3]],x26[[3]],x36[[3]],x46[[3]],x56[[3]],x66[[3]],x76[[3]])
xcondD6<-rbind(x16[[4]],x26[[4]],x36[[4]],x46[[4]],x56[[4]],x66[[4]],x76[[4]])
xW6<-rbind(x16[[5]],x26[[5]],x36[[5]],x46[[5]],x56[[5]],x66[[5]],x76[[5]])

dif<-seq(1.25,2.75,.25)/20-.1

plot(dif,xT6[,1],type="l",ylim=c(0,.75))
#points(dif,xD6[,1])
lines(dif,xT6[,1],col="red")
lines(dif,xD6[,1],col="blue")
lines(dif,xcondD6[,1],col="green")
lines(dif,xW6[,1],col="black")
lines(dif,xD6bias[,1],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")




xD7bias<-rbind(x17[[1]],x27[[1]],x37[[1]],x47[[1]],x57[[1]],x67[[1]],x77[[1]])
xT7<-rbind(x17[[2]],x27[[2]],x37[[2]],x47[[2]],x57[[2]],x67[[2]],x77[[2]])
xD7<-rbind(x17[[3]],x27[[3]],x37[[3]],x47[[3]],x57[[3]],x67[[3]],x77[[3]])
xcondD7<-rbind(x17[[4]],x27[[4]],x37[[4]],x47[[4]],x57[[4]],x67[[4]],x77[[4]])
xW7<-rbind(x17[[5]],x27[[5]],x37[[5]],x47[[5]],x57[[5]],x67[[5]],x77[[5]])

dif<-seq(1.25,2.75,.25)/20-.1

plot(dif,xT7[,1],type="l",ylim=c(0,.75))
#points(dif,xD7[,1])
lines(dif,xT7[,1],col="red")
lines(dif,xD7[,1],col="blue")
lines(dif,xcondD7[,1],col="green")
lines(dif,xW7[,1],col="black")
lines(dif,xD7bias[,1],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")








xD8bias<-rbind(x18[[1]],x28[[1]],x38[[1]],x48[[1]],x58[[1]],x68[[1]],x78[[1]])
xT8<-rbind(x18[[2]],x28[[2]],x38[[2]],x48[[2]],x58[[2]],x68[[2]],x78[[2]])
xD8<-rbind(x18[[3]],x28[[3]],x38[[3]],x48[[3]],x58[[3]],x68[[3]],x78[[3]])
xcondD8<-rbind(x18[[4]],x28[[4]],x38[[4]],x48[[4]],x58[[4]],x68[[4]],x78[[4]])
xW8<-rbind(x18[[5]],x28[[5]],x38[[5]],x48[[5]],x58[[5]],x68[[5]],x78[[5]])

dif1<-rep(0,7)
dif2<-seq(8.5,11.5,.5)/40-.25
dif3<-seq(3.5,6.5,.5)/40-.125
dif4<-rep(0,7)


plot(dif2,xD8[7:1,2],type="l",ylim=c(0,1))
#points(dif2,xT8[7:1,2])
lines(dif2,xT8[7:1,2],col="red")
lines(dif2,xD8[7:1,2],col="blue")
lines(dif2,xcondD8[7:1,2],col="green")
lines(dif2,xW8[7:1,2],col="black")
lines(dif2,xD8bias[7:1,2],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")



plot(dif3,xD8[,3],type="l",ylim=c(0,1))
#points(dif3,xT8[,3])
lines(dif3,xT8[,3],col="red")
lines(dif3,xD8[,3],col="blue")
lines(dif3,xcondD8[,3],col="green")
lines(dif3,xW8[,3],col="black")
lines(dif3,xD8bias[,3],col="purple",lty=3)
lines(-1:1,rep(.05,3),col="black")




xD9bias<-rbind(x19[[1]],x29[[1]],x39[[1]],x49[[1]],x59[[1]],x69[[1]],x79[[1]])
xT9<-rbind(x19[[2]],x29[[2]],x39[[2]],x49[[2]],x59[[2]],x69[[2]],x79[[2]])
xD9<-rbind(x19[[3]],x29[[3]],x39[[3]],x49[[3]],x59[[3]],x69[[3]],x79[[3]])
xcondD9<-rbind(x19[[4]],x29[[4]],x39[[4]],x49[[4]],x59[[4]],x69[[4]],x79[[4]])
xW9<-rbind(x19[[5]],x29[[5]],x39[[5]],x49[[5]],x59[[5]],x69[[5]],x79[[5]])





#######################
#Plots for Dissertation


#windows()    for windows
quartz()    #for mac


m <- matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE)

layout(mat = m,heights = c(0.8,0.2))



plot(seq(0,1,.001),dbeta(seq(0,1,.001),8.5,31.5),ylim=c(0,7),xlab="X1",ylab="Density",main=expression(delta==-.0375),type="l",lty=2,col="red")
lines(seq(0,1,.001),dbeta(seq(0,1,.001),5,15),lty=1,col="black")


plot(seq(0,1,.001),dbeta(seq(0,1,.001),10,30),ylim=c(0,7),xlab="X1",ylab="Density",main=expression(delta==0),type="l",lty=2,col="red")
lines(seq(0,1,.001),dbeta(seq(0,1,.001),5,15),lty=1,col="black")


plot(seq(0,1,.001),dbeta(seq(0,1,.001),11.5,28.5),ylim=c(0,7),xlab="X1",ylab="Density",main=expression(delta==.0375),type="l",lty=2,col="red")
lines(seq(0,1,.001),dbeta(seq(0,1,.001),5,15),lty=1,col="black")



plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("red","black")
legend(x = "top",inset = 0,
        legend = c("Population 1 (Prec = 40)", "Population 2 (Prec=20)"), 
        col=plot_colors,lty=c(2,1), lwd=2, cex=1, box.col="white", horiz = TRUE)










quartz()  #only for mac
#windows()  #only for windows.  Swap as necessary.
m <- matrix(c(1,2,3,4,4,4),nrow = 2,ncol = 3,byrow = TRUE)

layout(mat = m,heights = c(0.8,0.2))


plot(seq(0,1,.001),dbeta(seq(0,1,.001),3.5,36.5),ylim=c(0,10),xlab="X2",ylab="Density",main=expression(delta==-.0375),type="l",lty=2,col="red")
lines(seq(0,1,.001),dbeta(seq(0,1,.001),2.5,17.5),lty=1,col="black")


plot(seq(0,1,.001),dbeta(seq(0,1,.001),5,35),ylim=c(0,10),xlab="X2",ylab="Density",main=expression(delta==0),type="l",lty=2,col="red")
lines(seq(0,1,.001),dbeta(seq(0,1,.001),2.5,17.5),lty=1,col="black")


plot(seq(0,1,.001),dbeta(seq(0,1,.001),6.5,33.5),ylim=c(0,10),xlab="X2",ylab="Density",main=expression(delta==.0375),type="l",lty=2,col="red")
lines(seq(0,1,.001),dbeta(seq(0,1,.001),2.5,17.5),lty=1,col="black")


plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- c("red","black")
legend(x = "top",inset = 0,
        legend = c("Population 1 (Prec = 40)", "Population 2 (Prec=20)"), 
        col=plot_colors,lty=c(2,1), lwd=2, cex=1, box.col="white", horiz = TRUE)







dif1<-rep(0,7)
dif2<-seq(8.5,11.5,.5)/40-.25
dif3<-seq(3.5,6.5,.5)/40-.125
dif4<-rep(0,7)


plot(dif2,xD2[7:1,2],type="l",ylim=c(0,.40),main="Power Sim. for X1 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif2,xT2[7:1,2])
#lines(dif2,xT2[7:1,2],col="red")
lines(dif2,xD2[7:1,2],col="blue",lty=1,lwd=2)
points(dif2,xD2[7:1,2],col="blue",pch=14)

lines(dif2,xcondD2[7:1,2],col="green",lty=2,lwd=2)
points(dif2,xcondD2[7:1,2],col="green",pch=15)

#lines(dif2,xW2[7:1,2],col="black")
lines(dif2,xD2bias[7:1,2],col="purple",lty=3,lwd=2)
points(dif2,xD2bias[7:1,2],col="purple",pch=16)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(alpha==.05),side=2,line=1,at=.05)
mtext("Equal Sample Sizes n=20",side=3,line=.5)
legend("top",legend=c("Dirichlet", "DirNuiss","DirBias"),lty=c(1,2,3),col=c("blue","green","purple"),pch=c(14,15,16),cex=1.5)





#plotting to establish bias corrected versions will be used
windows()
par(mfrow=c(1,3))


plot(dif3,xD9[,3],type="l",ylim=c(0,.45),main="Power Sim. for X2 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif3,xT9[,3])
#lines(dif3,xT9[,3],col="red")
lines(dif3,xD9[,3],col="blue",lty=1,lwd=2)
points(dif3,xD9[,3],col="blue",pch=15)

lines(dif3,xcondD9[,3],col="green",lty=2,lwd=2)
points(dif3,xcondD9[,3],col="green",pch=16)

#lines(dif3,xW9[,3],col="black")
lines(dif3,xD9bias[,3],col="purple",lty=3,lwd=2)
points(dif3,xD9bias[,3],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-6,adj=1,at=-.045)
mtext("Equal Sample Sizes n=10",side=3,cex=.75,line=.5)
legend("top",legend=c("MLE", "MLE Null","MLE BC"),lty=c(1,2,3),lwd=2,col=c("blue","green","purple"),pch=c(15,16,17),cex=1.5)



plot(dif3,xD2[,3],type="l",ylim=c(0,.60),main="Power Sim. for X2 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif3,xT2[,3])
#lines(dif3,xT2[,3],col="red")
lines(dif3,xD2[,3],col="blue",lty=1,lwd=2)
points(dif3,xD2[,3],col="blue",pch=15)

lines(dif3,xcondD2[,3],col="green",lty=2,lwd=2)
points(dif3,xcondD2[,3],col="green",pch=16)

#lines(dif3,xW2[,3],col="black")
lines(dif3,xD2bias[,3],col="purple",lty=3,lwd=2)
points(dif3,xD2bias[,3],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-5,adj=1,at=-.045)
mtext("Equal Sample Sizes n=20",cex=.75,side=3,line=.5)
legend("top",legend=c("MLE", "MLE Null","MLE BC"),lty=c(1,2,3),lwd=2,col=c("blue","green","purple"),pch=c(15,16,17),cex=1.5)



plot(dif3,xD5[,3],type="l",ylim=c(0,1),main="Power Sim. for X2 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif3,xT5[,3])
#lines(dif3,xT5[,3],col="red")
lines(dif3,xD5[,3],col="blue",lty=1,lwd=2)
points(dif3,xD5[,3],col="blue",pch=15)

lines(dif3,xcondD5[,3],col="green",lty=2,lwd=2)
points(dif3,xcondD5[,3],col="green",pch=16)

#lines(dif3,xW5[,3],col="black")
lines(dif3,xD5bias[,3],col="purple",lty=3,lwd=2)
points(dif3,xD5bias[,3],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-4,adj=1,at=-.045)
mtext("Equal Sample Sizes n=100",side=3,cex=.75,line=.5)
legend("top",legend=c("MLE", "MLE Null","MLE BC"),lty=c(1,2,3),lwd=2,col=c("blue","green","purple"),pch=c(15,16,17),cex=1.5)



















#windows()
quartz()
par(mfrow=c(1,3))


plot(dif3,xT9[,3],type="l",ylim=c(0,.45),main="Power Sim. for X2 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif3,xT9[,3])
lines(dif3,xT9[,3],lty=1,lwd=2,col="red")
points(dif3,xT9[,3],col="red",pch=15)
#lines(dif3,xD9[,3],col="blue",lty=1,lwd=2)
#points(dif3,xD9[,3],col="blue",pch=15)

#lines(dif3,xcondD9[,3],col="green",lty=2,lwd=2)
#points(dif3,xcondD9[,3],col="green",pch=16)

lines(dif3,xW9[,3],lty=2,lwd=2,col="black")
points(dif3,xW9[,3],col="black",pch=16)
lines(dif3,xD9bias[,3],col="purple",lty=3,lwd=2)
points(dif3,xD9bias[,3],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-6,adj=1,at=-.045)
mtext("Equal Sample Sizes n=10",side=3,cex=.75,line=.5)
legend("top",legend=c("T", "RankSum","MLE BC"),lty=c(1,2,3),lwd=2,col=c("red","black","purple"),pch=c(15,16,17),cex=1.5)



plot(dif3,xT2[,3],type="l",ylim=c(0,.60),main="Power Sim. for X2 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif3,xT2[,3])
lines(dif3,xT2[,3],lty=1,lwd=2,col="red")
points(dif3,xT2[,3],col="red",pch=15)
#lines(dif3,xD2[,3],col="blue",lty=1,lwd=2)
#points(dif3,xD2[,3],col="blue",pch=15)

#lines(dif3,xcondD2[,3],col="green",lty=2,lwd=2)
#points(dif3,xcondD2[,3],col="green",pch=16)

lines(dif3,xW2[,3],lty=2,lwd=2,col="black")
points(dif3,xW2[,3],col="black",pch=16)
lines(dif3,xD2bias[,3],col="purple",lty=3,lwd=2)
points(dif3,xD2bias[,3],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-6,adj=1,at=-.045)
mtext("Equal Sample Sizes n=20",side=3,cex=.75,line=.5)
legend("top",legend=c("T", "RankSum","MLE BC"),lty=c(1,2,3),lwd=2,col=c("red","black","purple"),pch=c(15,16,17),cex=1.5)



plot(dif3,xT5[,3],type="l",ylim=c(0,1),main="Power Sim. for X2 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif3,xT5[,3])
lines(dif3,xT5[,3],lty=1,lwd=2,col="red")
points(dif3,xT5[,3],col="red",pch=15)
#lines(dif3,xD5[,3],col="blue",lty=1,lwd=2)
#points(dif3,xD5[,3],col="blue",pch=15)

#lines(dif3,xcondD5[,3],col="green",lty=2,lwd=2)
#points(dif3,xcondD5[,3],col="green",pch=16)

lines(dif3,xW5[,3],lty=2,lwd=2,col="black")
points(dif3,xW5[,3],col="black",pch=16)
lines(dif3,xD5bias[,3],col="purple",lty=3,lwd=2)
points(dif3,xD5bias[,3],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-6,adj=1,at=-.045)
mtext("Equal Sample Sizes n=100",side=3,cex=.75,line=.5)
legend("top",legend=c("T", "RankSum","MLE BC"),lty=c(1,2,3),lwd=2,col=c("red","black","purple"),pch=c(15,16,17),cex=1.5)











#windows()
quartz()
par(mfrow=c(1,3))


plot(dif2,xT9[7:1,2],type="l",ylim=c(0,.2),main="Power Sim. for X1 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif2,xT9[7:1,2])
lines(dif2,xT9[7:1,2],lty=1,lwd=2,col="red")
points(dif2,xT9[7:1,2],col="red",pch=15)
#lines(dif2,xD9[7:1,2],col="blue",lty=1,lwd=2)
#points(dif2,xD9[7:1,2],col="blue",pch=15)

#lines(dif2,xcondD9[7:1,2],col="green",lty=2,lwd=2)
#points(dif2,xcondD9[7:1,2],col="green",pch=16)

lines(dif2,xW9[7:1,2],lty=2,lwd=2,col="black")
points(dif2,xW9[7:1,2],col="black",pch=16)
lines(dif2,xD9bias[7:1,2],col="purple",lty=3,lwd=2)
points(dif2,xD9bias[7:1,2],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
#mtext(expression(.05),cex=.8,side=1,line=-6,adj=1,at=-.045)
mtext("Equal Sample Sizes n=10",side=3,cex=.75,line=.5)
legend("top",legend=c("T", "RankSum","MLE BC"),lty=c(1,2,3),lwd=2,col=c("red","black","purple"),pch=c(15,16,17),cex=1.5)



plot(dif2,xT2[7:1,2],type="l",ylim=c(0,.40),main="Power Sim. for X1 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif2,xT2[7:1,2])
lines(dif2,xT2[7:1,2],lty=1,lwd=2,col="red")
points(dif2,xT2[7:1,2],col="red",pch=15)
#lines(dif2,xD2[7:1,2],col="blue",lty=1,lwd=2)
#points(dif2,xD2[7:1,2],col="blue",pch=15)

#lines(dif2,xcondD2[7:1,2],col="green",lty=2,lwd=2)
#points(dif2,xcondD2[7:1,2],col="green",pch=16)

lines(dif2,xW2[7:1,2],lty=2,lwd=2,col="black")
points(dif2,xW2[7:1,2],col="black",pch=16)
lines(dif2,xD2bias[7:1,2],col="purple",lty=3,lwd=2)
points(dif2,xD2bias[7:1,2],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-7.5,adj=1,at=-.045)
mtext("Equal Sample Sizes n=20",side=3,cex=.75,line=.5)
legend("top",legend=c("T", "RankSum","MLE BC"),lty=c(1,2,3),lwd=2,col=c("red","black","purple"),pch=c(15,16,17),cex=1.5)



plot(dif2,xT5[7:1,2],type="l",ylim=c(0,1),main="Power Sim. for X1 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif2,xT5[7:1,2])
lines(dif2,xT5[7:1,2],lty=1,lwd=2,col="red")
points(dif2,xT5[7:1,2],col="red",pch=15)
#lines(dif2,xD5[7:1,2],col="blue",lty=1,lwd=2)
#points(dif2,xD5[7:1,2],col="blue",pch=15)

#lines(dif2,xcondD5[7:1,2],col="green",lty=2,lwd=2)
#points(dif2,xcondD5[7:1,2],col="green",pch=16)

lines(dif2,xW5[7:1,2],lty=2,lwd=2,col="black")
points(dif2,xW5[7:1,2],col="black",pch=16)
lines(dif2,xD5bias[7:1,2],col="purple",lty=3,lwd=2)
points(dif2,xD5bias[7:1,2],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-4.75,adj=1,at=-.045)
mtext("Equal Sample Sizes n=100",side=3,cex=.75,line=.5)
legend("top",legend=c("T", "RankSum","MLE BC"),lty=c(1,2,3),lwd=2,col=c("red","black","purple"),pch=c(15,16,17),cex=1.5)








#windows()
quartz()

par(mfrow=c(1,2))



plot(dif2,xT3[7:1,2],type="l",ylim=c(0,.3),main="Power Sim. for X1 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif2,xT3[7:1,2])
lines(dif2,xT3[7:1,2],lty=1,lwd=2,col="red")
points(dif2,xT3[7:1,2],col="red",pch=15)
#lines(dif2,xD3[7:1,2],col="blue",lty=1,lwd=2)
#points(dif2,xD3[7:1,2],col="blue",pch=15)

#lines(dif2,xcondD3[7:1,2],col="green",lty=2,lwd=2)
#points(dif2,xcondD3[7:1,2],col="green",pch=16)

lines(dif2,xW3[7:1,2],lty=2,lwd=2,col="black")
points(dif2,xW3[7:1,2],col="black",pch=16)
lines(dif2,xD3bias[7:1,2],col="purple",lty=3,lwd=2)
points(dif2,xD3bias[7:1,2],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
#mtext(expression(.05),cex=.8,side=1,line=-6,adj=1,at=-.045)
mtext("Unequal Sample Sizes 20,10",side=3,cex=.75,line=.5)
legend("top",legend=c("T", "RankSum","MLE BC"),lty=c(1,2,3),lwd=2,col=c("red","black","purple"),pch=c(15,16,17),cex=1.5)




plot(dif3,xT3[,3],type="l",ylim=c(0,.45),main="Power Sim. for X2 under Scenario 1",xlab=expression(delta),ylab="Est. Power")
#points(dif3,xT3[,3])
lines(dif3,xT3[,3],lty=1,lwd=2,col="red")
points(dif3,xT3[,3],col="red",pch=15)
#lines(dif3,xD3[,3],col="blue",lty=1,lwd=2)
#points(dif3,xD3[,3],col="blue",pch=15)

#lines(dif3,xcondD3[,3],col="green",lty=2,lwd=2)
#points(dif3,xcondD3[,3],col="green",pch=16)

lines(dif3,xW3[,3],lty=2,lwd=2,col="black")
points(dif3,xW3[,3],col="black",pch=16)
lines(dif3,xD3bias[,3],col="purple",lty=3,lwd=2)
points(dif3,xD3bias[,3],col="purple",pch=17)
lines(-1:1,rep(.05,3),col="black",lty=4)
mtext(expression(.05),cex=.8,side=1,line=-4,adj=1,at=-.045)
mtext("Unequal Sample Sizes 20,10",side=3,cex=.75,line=.5)
legend("top",legend=c("T", "RankSum","MLE BC"),lty=c(1,2,3),lwd=2,col=c("red","black","purple"),pch=c(15,16,17),cex=1.5)


