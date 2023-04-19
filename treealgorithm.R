library(MCMCpack)
library(grDevices)
library(gtools)
library(diagram)
#####Case of 4 variables

n<-1000
datasettry1<-rdirichlet(n,c(10,30))
datasettry2<-rdirichlet(n,c(4,4,4))
#datasettry3<-rdirichlet(n,c(25,10))


#datasettry<-rdirichlet(100,c(10,20,30,40))
datasettry<-cbind(datasettry1[,1],datasettry1[,2]*datasettry2)
#datasettry<-cbind(datasettry1[,1]*datasettry3,datasettry1[,2]*datasettry2)
datasettry<-datasettry[,c(2,3,1,4)]

n<-dim(datasettry)[[1]]
n_var<-dim(datasettry)[[2]]

j<-floor(n_var/2)

combin_vec<-c()
for(i in 1:j){
combin_vec[i]<-choose(n_var,i)
}


Simple<-dirichletmlealphas(datasettry,ind=1,N=40)
Criterion<-c()
Criterion[1]<-loglike3(Simple/sum(Simple),sum(Simple),datasettry)


#For additional node
for (i in 1:n_var){
alpha1<-dirichletmlealphas(cbind(datasettry[,i],apply(datasettry[,-i],1,sum)),ind=1,N=40)
alpha2<-dirichletmlealphas(datasettry[,-i]/apply(datasettry[,-i],1,sum),ind=1,N=40)
Criterion[i+1]<-loglike3(alpha1/sum(alpha1),sum(alpha1),cbind(datasettry[,i],apply(datasettry[,-i],1,sum)))
#+loglike3(alpha2/sum(alpha2),sum(alpha2),datasettry[,-i]/apply(datasettry[,-i],1,sum))
}

x<-combinations(4,2)
for(i in 1:6){
alpha1<-dirichletmlealphas(cbind(apply(datasettry[,x[i,]],1,sum),apply(datasettry[,-x[i,]],1,sum)),ind=1,N=40)
alpha2<-dirichletmlealphas(datasettry[,x[i,]]/apply(datasettry[,x[i,]],1,sum),ind=1,N=40)
alpha3<-dirichletmlealphas(datasettry[,-x[i,]]/apply(datasettry[,-x[i,]],1,sum),ind=1,N=40)
Criterion[i+5]<-loglike3(alpha1/sum(alpha1),sum(alpha1),cbind(apply(datasettry[,x[i,]],1,sum),apply(datasettry[,-x[i,]],1,sum)))
#+loglike3(alpha2/sum(alpha2),sum(alpha2),datasettry[,x[i,]]/apply(datasettry[,x[i,]],1,sum))
#+loglike3(alpha3/sum(alpha3),sum(alpha3),datasettry[,-x[i,]]/apply(datasettry[,-x[i,]],1,sum))
}

#k<-c(4,5,5,5,5,6,6,6,6,6,6)
k=2
AIC<-2*k-2*Criterion
AICc<-AIC+2*k*(k+1)/(n-k-1)
BIC<--2*Criterion+k*log(n)

-2*Criterion
AIC
AICc
BIC




####Case of 6 variables
#Simulating a nested dirichlet distribution scenario
n<-100
datasettry1<-rdirichlet(n,c(10,10,10))
datasettry2<-rdirichlet(n,c(3,7)*2)
datasettry3<-rdirichlet(n,c(5,5,10))

dataset<-cbind(datasettry1[,1]*datasettry2,datasettry1[,2]*datasettry3,datasettry1[,3])
#dataset<-cbind(datasettry1[,1]*datasettry2,datasettry1[,2]*datasettry3)
dataset<-as.matrix(dataset)
dataset<-data.frame(dataset)
names(dataset)












library(gtools)
#load your own functions
library(diagram)


split.tree<-function(index,dataset,method="MaxLik"){
dataset.dim<-dim(dataset)[2]
index.vec<-index[!index==0]
indicator.vec<-1:length(index.vec)

if(length(index.vec)<3){return(matrix(rep(0,2*length(index)),nrow=2))}

if(length(index.vec)>2){
datasettry<-dataset[,index.vec]/apply(dataset[,index.vec],1,sum)


n<-dim(datasettry)[[1]]
n_var<-dim(datasettry)[[2]]

j<-floor(n_var/2)

combin_vec<-c()
for(i in 1:j){
combin_vec[i]<-choose(n_var,i)
}

total<-1+sum(combin_vec)
Criterion<-c()

Simple<-dirichletmlealphas(datasettry,ind=1,N=40)
Criterion[1]<-loglike3(Simple/sum(Simple),sum(Simple),datasettry)

start<-1

labels<-matrix(rep(0,j),1,j)

for(l in 1:j){
labels<-rbind(labels,cbind(combinations(n_var,l),matrix(rep(0,(j-l)*dim(combinations(n_var,l))[1]),nrow=dim(combinations(n_var,l))[1])))
}

for(l in 1:j){
x<-combinations(n_var,l)


if(l==1){
for(i in 1:(combin_vec[l])){
alpha1<-dirichletmlealphas(cbind(apply(matrix(datasettry[,x[i,]],nrow=n),1,sum),apply(datasettry[,-x[i,]],1,sum)),ind=1,N=40)
alpha2<-dirichletmlealphas(matrix(datasettry[,-x[i,]],nrow=n)/apply(datasettry[,-x[i,]],1,sum),ind=1,N=40)
Criterion[start+i]<-loglike3(alpha1/sum(alpha1),sum(alpha1),cbind(apply(matrix(datasettry[,x[i,]],nrow=n),1,sum),apply(datasettry[,-x[i,]],1,sum)))+
			  loglike3(alpha2/sum(alpha2),sum(alpha2),matrix(datasettry[,-x[i,]],nrow=n)/apply(datasettry[,-x[i,]],1,sum))

}
}

if(l>1){
for(i in 1:(combin_vec[l])){
alpha1<-dirichletmlealphas(cbind(apply(matrix(datasettry[,x[i,]],nrow=n),1,sum),apply(datasettry[,-x[i,]],1,sum)),ind=1,N=40)
alpha2<-dirichletmlealphas(matrix(datasettry[,x[i,]],nrow=n)/apply(datasettry[,x[i,]],1,sum),ind=1,N=40)
alpha3<-dirichletmlealphas(matrix(datasettry[,-x[i,]],nrow=n)/apply(datasettry[,-x[i,]],1,sum),ind=1,N=40)
Criterion[start+i]<-loglike3(alpha1/sum(alpha1),sum(alpha1),cbind(apply(matrix(datasettry[,x[i,]],nrow=n),1,sum),apply(datasettry[,-x[i,]],1,sum)))+
			 loglike3(alpha2/sum(alpha2),sum(alpha2),matrix(datasettry[,x[i,]],nrow=n)/apply(datasettry[,x[i,]],1,sum))+
			 loglike3(alpha3/sum(alpha3),sum(alpha3),matrix(datasettry[,-x[i,]],nrow=n)/apply(datasettry[,-x[i,]],1,sum))
}
}

start<-length(Criterion)
}
k=c(n_var,rep(n_var+1,n_var),rep(n_var+2,dim(labels)[1]-n_var-1))

neg2Crit<-Criterion*(-2)
AIC<- 2*k-2*Criterion
AICc<- AIC+2*k*(k+1)/(n-k-1)
BIC<- -2*Criterion+k*log(n)

if(method=="MaxLik"){   loc<-which(neg2Crit==min(neg2Crit))[1]}
if(method=="AIC"){   loc<-which(AIC==min(AIC))[1]}
if(method=="AICc"){   loc<-which(AICc==min(AICc))[1]}
if(method=="BIC"){   loc<-which(BIC==min(BIC))[1]}

com.vec<-1:n_var
subset1<-labels[loc,]
subset2<-com.vec[-which(com.vec%in%subset1)]
subset1.final<-index.vec[subset1]
subset2.final<-index.vec[subset2]

mymat<-rbind(c(subset1.final,rep(0,dataset.dim-length(subset1.final))), c(subset2.final,rep(0,dataset.dim-length(subset2.final))))
#print(cbind(labels,-2*Criterion,AIC,AICc,BIC))

return(list(splitmatirx=mymat, criteria = neg2Crit))

}
}



directions.next<-function(x,mat.directions){
      cbind(rep(mat.directions[,2],each=2), (2^(x)):(2^(x+1)-1)  )
}


#directions.next(2,xtry)
#directions.next(3,directions.next(2,xtry))




LDM.alg<-function(dataset,meth="MaxLik"){

#dataset<-healthy.full
data<-as.matrix(dataset)


direction.mat<-c()
splits.mat<-c()
splits.temp<-c()
direction.temp<-c()
dummy<-c()
j=2
first.index<-1:dim(data)[2]


splits.mat<-rbind(splits.mat,split.tree(first.index,data,method=meth))
splits.temp<-splits.mat
total<-sum(splits.mat)
if(total!=0){direction.temp<-matrix(c(1,2,1,3),byrow=T,2,2)
             direction.mat<-direction.temp
             }

while(total!=0){
     dummy<-c()
     for (i in 1:dim(splits.temp)[1]){
         dummy<-rbind(dummy,split.tree(splits.temp[i,],data,method=meth))
          }
     total<-sum(dummy)
     splits.temp<-dummy
     splits.mat<-rbind(splits.mat,splits.temp)
     direction.temp<-directions.next(j,direction.temp)
     direction.mat<-rbind(direction.mat,direction.temp)
     j=j+1
}     


#windows()
quartz()
 par(mar=c(1,1,1,1))

 openplotmat()

 elpos<-coordinates(2^(0:(j-2)))
fromto<-direction.mat[-which(apply(splits.mat,1,sum)==0),]


 nr     <-nrow(fromto)

 arrpos <- matrix(ncol=2,nrow=nr)

 for (i in 1:nr) 
     arrpos[i,]<- straightarrow (to=elpos[fromto[i,2],],from=elpos[fromto[i,1],]
         ,lwd=2,arr.pos=0.6,arr.length=0.5)

 textellipse(elpos[1,],0.05,lab="Total",box.col="grey",shadow.col="black",shadow.size=0.005,cex=.9)

nodes.final<-fromto[which(fromto[,2]%in%fromto[,1]==FALSE),2]

for (i in 1:length(nodes.final)){
 label.name<-names(dataset)[splits.mat[nodes.final-1,][i,][which(splits.mat[nodes.final-1,][i,]!=0)]  ]
 textellipse(elpos[nodes.final[i],],0.05,lab=paste(label.name,collapse=' '),box.col="white",shadow.col="black",shadow.size=0.005,cex=.9)
}

return(splits.mat)
}


LDM.alg(dataset,meth="MaxLik")
LDM.alg(dataset,meth="AIC")
LDM.alg(dataset,meth="BIC")

