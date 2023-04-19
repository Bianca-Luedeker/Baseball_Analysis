loglike3<-function(mu,s,data){
N<-dim(data)[1]
logG<-apply(log(data),2,mean)

return(N*lgamma(s) - N*sum(lgamma(s*mu)) + N*sum((mu*s)*logG) )
}

#s = A, mu = pi 
loglike4<-function(mu,s1,s2,data1,data2){

N1<-dim(data1)[1]
N2<-dim(data2)[1]

logG1<-apply(log(data1),2,mean)
logG2<-apply(log(data2),2,mean)
loglike<-N1*lgamma(s1)+N1*sum( mu*s1*logG1)- N1*sum(lgamma(s1*mu))+N2*lgamma(s2)+N2*sum( mu*s2*logG2)-N2*sum(lgamma(s2*mu))

return(loglike)
}



#LRT test for the one node of a tree

tree.test<-function(data1,data2){

df.full<-2*dim(data1)[2]
df.null<-df.full-dim(data1)[2]-1
mle.null<-dirichletnull(data1,data2,50)
mle.null.mu<-mle.null[[1]]
mle.null.s1<-mle.null[[2]]
mle.null.s2<-mle.null[[3]]

mle.full1<-dirichletmlealphas(data1,ind=1,N=40)
mle.full2<-dirichletmlealphas(data2,ind=1,N=40)

lrt.stat<- loglike4(mle.null.mu,mle.null.s1,mle.null.s2,data1,data2)-loglike3(mle.full1/sum(mle.full1),sum(mle.full1),data1)-loglike3(mle.full2/sum(mle.full2),sum(mle.full2),data2)

lrt.stat<--2*lrt.stat
return(lrt.stat)
#return(list(loglike4(mle.null.mu,mle.null.s1,mle.null.s2,data1,data2),loglike3(mle.full1/sum(mle.full1),sum(mle.full1),data1)+loglike3(mle.full2/sum(mle.full2),sum(mle.full2),data2)))
}








 dummy<-c()
 dummy2<-c()
 dummy3<-c()
 for (i in 1:1000){
 dummy[i]<-tree.test(rdirichlet(5,c(20,10,5,5)/2),rdirichlet(5,c(20,10,5,5)))
 dummy2[i]<-tree.test(rdirichlet(20,c(20,10,5,5)/2),rdirichlet(20,c(20,10,5,5)))
 dummy3[i]<-tree.test(rdirichlet(100,c(20,10,5,5)/2),rdirichlet(100,c(20,10,5,5)))

}

 hist(dummy,freq=F,ylim=c(0,.25))
 lines(seq(0,20,100),dchisq(seq(0,20,100),3))
 
 hist(dummy2,freq=F,ylim=c(0,.25))
 lines(seq(0,20,100),dchisq(seq(0,20,100),3))
 
 hist(dummy3,freq=F,ylim=c(0,.25))
 lines(seq(0,20,100),dchisq(seq(0,20,100),3))

 hist(dummy+dummy2+dummy3,freq=F,ylim=c(0,.1))
 lines(seq(0,30,100),dchisq(seq(0,30,100),9))

