################################################################################
## Bianca Luedeker
## June 11, 2022
## Analysis of the Baseball Data Using Nested Dirichlet Techniques
## Updated:  9/19/2022
## The update makes better ternary diagrams.
################################################################################

#### Pacakages ####
library(dplyr)        ## For Data management.
library(xtable)       ## Creates tables for latex.
library(gtools)        ##Generates all possible permutations.  
library(sirt)          ##Estimates Dirichlet alpha parameters
library(compositions)  ##Ternary Diagrams.

#### Organize and Clean the Data ####

## Read in Data For batters
batters.data <- read.table("batters_by_year.txt", header=FALSE, sep=",")
View(batters.data)

##Data cleaning 

## data_2: Intermediate data set with ID, birth year, and number of each type of
## outcome with no difference between handiness of pitcher.  
## I had to combine rows where handiness of pitcher was distinguished.  
## This is no longer in the data_2 data set.
data_2 <- batters.data %>%                                    
  group_by(V1, V2) %>%
  dplyr::summarise( "BIRTH" = mean(V3), "INT" = sum(V7), "HBP"= sum(V8),
                   "BB" = sum(V9), "KO" = sum(V10), "FBHR" = sum(V11), 
                   "FB3B"= sum(V12), "FB2B" = sum(V13), "FB1B" = sum(V14),
                   "FB0" = sum(V15), "GBHR" = sum(V16), "GB3B" = sum(V17),
                   "GB2B" = sum(V18), "GB1B" = sum(V19), "GB0" = sum(V20)) %>% 
  as.data.frame()
View(data_2)

## data_3: Intermediate data set with ID, birth year, season year
## and six outcome counts: Home Run, Third Base, Second Base, First Base,
## Out, and Other.
data_3 <- data_2 %>%
  dplyr::mutate("HR" = FBHR + GBHR, "3B" = FB3B + GB3B, "2B" = FB2B + GB2B,
                "1B" = FB1B + GB1B, "OUT" = KO + FB0 + GB0, 
                "OTHER" = INT + HBP + BB, 
                "TOTAL" = INT+HBP+BB+KO+FBHR+FB3B+FB2B+FB1B+FB0+GBHR+GB3B+
                  GB2B+GB1B+GB0) %>%
  select(1:3, 18:24) %>%
  as.data.frame()
View(data_3)

## data_4:  Intermediate data set with ID, Birth Year, Season Year, and
## proportions of six outcomes:  percent home run, percent third base, 
## percent second base, percent first base, percent out, and percent other.
data_4 <- data_3 %>%
  dplyr::mutate( "AGE" = V2 - BIRTH, "HR" = HR/TOTAL, "3B" = `3B`/TOTAL, 
                 "2B" =`2B`/TOTAL, "1B" = `1B`/TOTAL, "OUT" = OUT/TOTAL, 
                 "OTHER" = OTHER/TOTAL) %>%
  dplyr::rename("ID"= V1, "Year" = V2) %>%
  select(1, 11, 4:10) %>%
  as.data.frame()
View(data_4)

## Look at ages to establish groups.  Make three equal sized groups.
## All ages from 20 - 49 are represented.  Group 1: 20 - 24, Group 2: 25- 34,
## Group 3: 35-49
sort(unique(data_4[,2])) 
temp <- sort(data_4[,2])
length(temp)/3
temp[2700]
temp[5400]



data_4$GROUP <- ifelse(data_4$AGE < 25, 1, ifelse(data_4$AGE <35, 2, 3))
View(data_4)
sum(data_4$GROUP==1)  ##1063 batters
sum(data_4$GROUP==2)  ## 6021 batters
sum(data_4$GROUP==3)  ## 1015 batters

batters.group.dat <- data_4 %>% relocate(GROUP, .before = HR)
View(batters.group.dat)

write.csv(batters.group.dat, "batters_grouped.csv")

################################################################################
## June 18, 2022
## Now that we have the data, we need to find the tree.
################################################################################

batter.dat <- read.csv("batters_grouped.csv")[,2:11]
batter.dat <- batter.dat %>% rename( "TRIPLE" = X3B,
                       "DOUBLE" = X2B,
                       "SINGLE" = X1B)
View(batter.dat)

age.1.data <- as.matrix(batter.dat[batter.dat[,3]==1, 4:9])
age.2.data <- as.matrix(batter.dat[batter.dat[,3]==2, 4:9])
age.3.data <- as.matrix(batter.dat[batter.dat[,3]==3, 4:9])

## Check Correlation Structure
cor(age.1.data)
cor(age.2.data)
cor(age.3.data)
aa <- round(cor(batter.dat[, 4:9]), 4)
xtable(aa, type=latex, digits = c(0, rep(4, 6)))

## Correlation matrices for all three groups.

b1 <-  round(cor(age.1.data), 4)
b2 <-  round(cor(age.2.data), 4)
b3 <-  round(cor(age.3.data), 4)

xtable(b1, type=latex, digits = c(0, rep(4, 6)))
xtable(b2, type=latex, digits = c(0, rep(4, 6)))
xtable(b3, type=latex, digits = c(0, rep(4, 6)))

## Check differences between correlations.
round(cor(age.1.data)-cor(age.2.data), 2)
round(cor(age.1.data)-cor(age.3.data), 2)
round(cor(age.2.data) - cor(age.3.data), 2)

## Create Ternary Diagrams  ##
batter.dat <- read.csv("batters_grouped.csv")[,2:11]
batter.dat <- batter.dat %>% rename( "T" = X3B,
                                     "D" = X2B,
                                     "S" = X1B)
head(batter.dat)
YY <- acomp(batter.dat[, 4:9])
XX <- factor(batter.dat[,3], ordered = FALSE)


water.maze.data <- read.csv("water_maze_data.csv")
YY <- acomp(water.maze.data[,2:5])  ##The composition is the dependent variable YY.
XX <- factor(water.maze.data[,1], ordered=FALSE)  ##The only covariable: group membership.

plot(YY, pch=20, col = c("red", "blue", "black")[XX], main="3 Part Subcompositions for 3tg and Wild Groups")
# legend(locator(1), levels(XX), pch=20, col=c("red", "blue"), xpd=NA,  yjust=0)

## This was completed on the HPc below.  See following code.


#### Find the tree diagram for this example.  ####

## June 19, 2022
## This section contains all the functions necessary to find a tree diagram.

##  The Log-Likelhood function ####

log.like.1 <- function(alpha, data){
  n <- nrow(data)
  log.p.bar <- apply(log(data), 2, mean)
  result <- n*lgamma(sum(alpha)) - n*sum(lgamma(alpha)) + n*sum(alpha*log.p.bar)  ##Correct version
  return(result)
}

## One-Step of the binary tree algorithm ##

tree.finder <- function(index, dataset){
  eps <- 10^-5
  dataset <- (dataset+eps)/(1+2*eps) ## Protects against zero entries.
  dataset <- dataset/rowSums(dataset)
  orig.k <- ncol(dataset)  ##number of parts of the composition.
  index.vec <- index[index!=0]  ## removes 0 entries from index vector.
  indicator.vec <- 1:length(index.vec) ## numbers 1 through length of index vector with 0 removed
  
  ## Note:  The index vector gives us the columns to pull from the data set.  
  ## Which of the parts of the composition to use.  Or How many columns to use.
  ## A zero indicates that part of the composition has been dropped.
  
  if(length(index.vec) <= 2){
    return(matrix(rep(0, 2*length(index)), nrow=2))
  }
  
  else{
    new.data <- dataset[, index.vec]/rowSums(dataset[, index.vec])  ##Rescale to get branch proportions.
    n <- nrow(new.data)
    k <- ncol(new.data)
    j <- floor(k/2)  ## Two groups wiil be created.  k=7, (1, 6), (2, 5), (3,4)
    ## We only need the bottom half since (4,3), (5,2), (6, 1) are the same.  
    ## hence the floor of k/2
    
    combin.vec <- c()
    for(i in 1:j){
      combin.vec[i] <- choose(k, i)  ##Number of combinations for 1, 2, ....
    }
    
    total <- 1 + sum(combin.vec)  ## the extra one is for the standard DD.  
    ##Note: this double counts trees at the midway point since (1,2) (3, 4) is 
    ## the same as (3, 4) (1, 2).  Some inefficiency.
    
    criterion <- c()
    standard.alpha <- sirt::dirichlet.mle(new.data)$alpha
    criterion[1] <- log.like.1(standard.alpha, new.data)
    
    start <- 1
    
    ##The matrix labels tells which nodes are in one group.  Nodes not listed
    ## are in the other group.
    labels <-  matrix(rep(0,j), 1, j)  ##matrix ensures good formating
    for(i in 1:j){
      a <- combinations(k, i)
      labels <- rbind(labels, 
                      cbind(a, matrix(rep(0, (j-i)*nrow(a)), nrow = nrow(a))))
      
      ##Case 1: The loner node case: One internal node.
      if(i==1){
        for(m in 1:combin.vec[1]){
          level.1 <- cbind(new.data[,m], 1-new.data[,m])
          alpha.1 <- sirt::dirichlet.mle(level.1)$alpha
          level.2 <- new.data[, -m]/(1-new.data[,m])
          alpha.2 <- sirt::dirichlet.mle(level.2)$alpha
          criterion[start + m] <- log.like.1(alpha.1, level.1) + log.like.1(alpha.2, level.2)
        }
      }
      
      ## Case 2: Two internal nodes
      else{
        for(m in 1:combin.vec[i]){
          level.1 <- cbind(rowSums(new.data[, a[m,]]),  rowSums(new.data[, -a[m,]]))
          alpha.1 <- sirt::dirichlet.mle(level.1)$alpha
          level.2.left <- new.data[, a[m,]]/rowSums(new.data[, a[m,]])
          alpha.2.left <- sirt::dirichlet.mle(level.2.left)$alpha
          level.2.right <- new.data[, -a[m,]]/rowSums(new.data[, -a[m,]])
          alpha.2.right <- sirt::dirichlet.mle(level.2.right)$alpha
          criterion[start + m] <- log.like.1(alpha.1, level.1) +
            log.like.1(alpha.2.left, level.2.left) + log.like.1(alpha.2.right, level.2.right)
        }
        
      }
      start <- length(criterion)
    }
    
    bb <- c(k, rep(k+1, k), rep(k+2, nrow(labels)-k-1))
    neg2.crit <- -2*criterion
    best <- which(neg2.crit==min(neg2.crit))[1]  ##find the min, use the first if there are ties
    
    com.vec <- 1:k
    left.best <- labels[best ,]  ##nodes in left branch for the best
    right.best <- com.vec[-which(com.vec %in% left.best)]  ##Careful in the case of the non-nested design.  
    ##If no nest is the best you get left.best <- (0,0...)
    ## Right best then is integer 0.  Might need another if statement here.
    left.final <- index.vec[left.best]
    right.final <- index.vec[right.best]  ## Same problems here.
    final.matrix <- rbind( 
      c(left.final, rep(0, orig.k - length(left.final))),
      c(right.final, rep(0, orig.k - length(right.final)))
    )
    #return(list( labels = labels, criteria = neg2.crit, what=bb, loc=best, left=left.best, right=right.best, mat=final.matrix))
    return(final.matrix)
  }
}

##Complete Binary tree function ##
## Uses tree.finder to create a binary cascade tree where only one or two
## nodes are at the end of each branch.


complete.bin.tree <- function(data.set){
  eps <- 10^-5
  data.set <- (data.set+eps)/(1+2*eps)
  data.set <- data.set/rowSums(data.set)
  data <- as.matrix(data.set)  ##puts data into correct form
  
  ## Initialize vectors
  direction.mat <- c()
  splits.mat <- c()
  splits.temp <- c()
  dummy <- c()
  
  j = 1
  first.index <- 1:ncol(data)  ##for the first binary split, we are going to use all the variables.
  
  splits.mat <- rbind(splits.mat, tree.finder(first.index, data))  ##Records left and right splits
  splits.temp <- splits.mat
  total <- sum(splits.mat)  ##If total equals 0, the non-nested design is the winner
  
  ## For cases where a nested design wins.
  
  if( total != 0){
    direction.temp <- matrix( c(1,2,1,3), byrow = TRUE, nrow=2, ncol=2)  ## What does this do?
    direction.mat <- direction.temp
  }
  
  ## The iterative step.
  while ( total != 0){
    dummy <- c() 
    ## for each row of current splits.  Note there are 2*number of divisions in splits.
    ## one for the left branch and one for the right branch
    for (i in 1:nrow(splits.temp)){ 
      dummy <- rbind( dummy, tree.finder(splits.temp[i,], data))  ##does one more iteration on each split
    }
    
    total <- sum(dummy)
    splits.temp <- dummy
    splits.mat <- rbind(splits.mat, splits.temp)
    j = j+1
  }
  non.zero <- which(rowSums(splits.mat) != 0)
  splits.mat <- splits.mat[non.zero ,]
  return(splits.mat)
}

#### Finding the correct tree for the baseball data ####

complete.bin.tree(batter.dat[,4:9])



#### June 26, 2022 ####
## Ternary diagrams again.
## Ternary 

batter.dat <- read.csv("batters_grouped.csv")[,2:11]
batter.dat <- batter.dat %>% rename( "TRI" = X3B,
                                     "DBL" = X2B,
                                     "SGL" = X1B,
                                     "OTR" = OTHER)

age.1.data <- as.matrix(batter.dat[batter.dat[,3]==1, 4:9])
age.2.data <- as.matrix(batter.dat[batter.dat[,3]==2, 4:9])
age.3.data <- as.matrix(batter.dat[batter.dat[,3]==3, 4:9])

YY <- acomp(age.1.data + 0.00001)
jpeg("fig_8_1_A.jpeg", width =800, height = 800, quality = 100)
plot(YY, pch=20, cex=1,  main="Age Group 1 Ternary Diagrams")
dev.off()


YY <- acomp(age.2.data + 0.00001)
jpeg("fig_8_1_B.jpeg", width =800, height = 800, quality = 100)
plot(YY, pch=20, cex=1,  main="Age Group 2 Ternary Diagrams")
dev.off()

YY <- acomp(age.3.data + 0.00001)
jpeg("fig_8_1_C.jpeg", width =800, height = 800, quality = 100)
plot(YY, pch=20, cex=1,  main="Age Group 3 Ternary Diagrams")
dev.off()

#### July 16, 2022  #####
## Now that we have a nesting tree, check the correlation structure.
batter.dat <- read.csv("batters_grouped.csv")[,2:11]
batter.dat <- batter.dat %>% rename( "TRI" = X3B,
                                     "DBL" = X2B,
                                     "SGL" = X1B,
                                     "OTR" = OTHER)
## Determine alpha values.

level.one <- cbind(rowSums(batter.dat[,4:6]), rowSums(batter.dat[,7:9]))
head(level.one)
head(batter.dat)
level.one.alpha <- sirt::dirichlet.mle(level.one)$alpha

level.two.a <- cbind(batter.dat[,4]+batter.dat[,5], batter.dat[,6])
level.two.a.alpha <- sirt::dirichlet.mle(level.two.a)$alpha

level.two.b <- cbind(batter.dat[,8], batter.dat[,7]+batter.dat[,9])
level.two.b.alpha <- sirt::dirichlet.mle(level.two.b)$alpha

level.three.a <- cbind(batter.dat[,4], batter.dat[,5])
level.three.a.alpha <- sirt::dirichlet.mle(level.three.a)$alpha

level.three.b <- cbind(batter.dat[,7], batter.dat[,9])
level.three.b.alpha <- sirt::dirichlet.mle(level.three.b)$alpha

## Create a data set with 10,000 observations.

set.seed(1983)
level.one.dat <- gtools::rdirichlet(10000, level.one.alpha)
level.two.a.dat <- gtools::rdirichlet(10000, level.two.a.alpha)
level.two.b.dat <- gtools::rdirichlet(10000, level.two.b.alpha)
level.three.a.dat <- gtools::rdirichlet(10000, level.three.a.alpha)
level.three.b.dat <- gtools::rdirichlet(10000, level.three.b.alpha)

sim.baseball.dat <-
  cbind("HR" = level.one.dat[,1]*level.two.a.dat[,1]*level.three.a.dat[,1],
        "TRI" = level.one.dat[,1]*level.two.a.dat[,1]*level.three.a.dat[,2],
        "DBL"= level.one.dat[,1]*level.two.a.dat[,2],
        "SGL" = level.one.dat[,2]*level.two.b.dat[,2]*level.three.b.dat[,1],
        "OUT" =  level.one.dat[,2]*level.two.b.dat[,1],
        "OTR" = level.one.dat[,2]*level.two.b.dat[,2]*level.three.b.dat[,2])

## Check the row sums and column means.

rowSums(sim.baseball.dat)[1:200]
colMeans(sim.baseball.dat)
colMeans(batter.dat[,4:9])

## Check correlation structure

aa <- round(cor(batter.dat[, 4:9]), 3)
bb <- round(cor(sim.baseball.dat), 3)

xtable(aa, type=latex, digits = c(0, rep(3, 6)))
xtable(bb, type=latex, digits=c(0, rep(3,6)))

#### Test For Differences ####

## Re-paste all relevant functions here.

##Null Log-likelihood multiple groups ##
null.mean.ll.multi <- function(z, A, n, log.x.bar){
  j <- length(n)
  res <- 0
  pi.est <- z/sum(z)
  for(i in  1:j){
    res <- res+n[i]*sum(A[i]*pi.est*log.x.bar[i,])-n[i]*sum(lgamma(A[i]*pi.est))
  }
  return(-1*res)
}

## Precision Maximization Function ##
precision.max <- function(pi.est, A.est, log.x.bar, convcrit = 1e-05, maxit=100){
  A <- A.est
  conv <- 1
  iter <- 1
  while((conv > convcrit) & (iter < maxit)){
    
    inv.new.A <- 1/A + 
      1/A^2*(digamma(A) - sum(pi.est*(digamma(A*pi.est)-log.x.bar)))/
      (trigamma(A) - sum(pi.est^2*trigamma(A*pi.est)))
    
    new.A <- 1/inv.new.A
    conv <- abs(new.A - A)
    iter <- iter + 1
    A <- new.A
  }
  return(A)
}

## Gradient of Log-Likelihood ##

multi.gradient <- function(z, A, n, log.x.bar){
  j <- length(n)
  k <- ncol(log.x.bar)
  res <- rep(0, k)
  pi.est <- z/sum(z)
  for(i in 1:j){
    temp <- n[i]*A[i]/sum(z)*(log.x.bar[i,] - digamma(A[i]*pi.est)-
                                sum(pi.est*(log.x.bar[i,]-digamma(A[i]*pi.est))))
    res <- res +temp
  }
  return(res)
}

##Maximization function ##

null.ll.max.multi<- function(data.set, convcrit = 1e-05, maxit=50){
  k <- ncol(data.set)-1
  group.id <- sort(unique(data.set[,1]))  ## Added July 18, 2022.  Discovered BUG.
  j <- length(group.id)
  pi.est <- z.est <- colMeans(data.set [, 2:(k+1)])
  n <- vector ()
  log.x.bar <- vector()
  A.est <- vector()
  ## This loops is for determining the values of n, log.x.bar, initial 
  ## precisions, and initial mean vectors.
  #for(i in 1:j){  Replaced July 18, 2022
  for(i in 1:j){
    ID <- group.id[i]
    my.data <- data.set[data.set[,1]==ID, 2:(k+1)]
    n <- c(n, nrow(my.data))
    log.x.bar <- rbind(log.x.bar, colMeans(log(my.data)))
    pi.initial <-  colMeans(my.data)
    p2 <- mean(my.data[,1]^2)
    A.est <- c(A.est, (pi.initial[1] - p2)/(p2 - (pi.initial[1])^2))
  }
  
  conv <- 1
  iter <- 1
  A.new <- rep(0, j)
  while((conv > convcrit) & (iter < maxit)){
    for(i in 1:j){
      A.new[i] <- precision.max(pi.est, A.est[i], log.x.bar[i,])
    }
    z.new <- optim(z.est, null.mean.ll.multi,  gr=multi.gradient, A=A.new, n=n, 
                   log.x.bar=log.x.bar, method = "L-BFGS-B", 
                   lower= rep(0.001, k))$par
    ##Note: z.new was changed on 5/2/2022 when the Nelder-Mead method failed for 
    ## a mean vector with components near zero.
    
    conv <- max(c(abs(A.new - A.est), abs(z.new-z.est)))
    iter <- iter + 1
    A.est <- A.new
    z.est <- z.new
    pi.est <- z.new/sum(z.new)
  }
  ## Compute the max under the null.
  
  return(list(precisions=A.est, mean=pi.est, iter=iter, n=n, 
              log.x.bar=log.x.bar))
}

## LRT Statistic Function ##
## This function was changed on May 31, 2022 to return the LRT statistic
## instead of a p-value for multiple layers.  It was given the name 
## my.lrt.fun.3
my.lrt.fun.3 <- function(data.set){
  temp <- null.ll.max.multi(data.set)
  log.x.bar <- temp$log.x.bar
  n <- temp$n
  group.id <- sort(unique(data.set[,1]))  ## Added July 18, 2022
  j <- length(group.id)
  A <- temp$precisions
  pi.est <- temp$mean
  k <- length(pi.est)
  null.max <- 0
  alt.max <- 0
  for(i in 1:j){
    null.max <- null.max + (n[i]*lgamma(A[i])-n[i]*sum(lgamma(A[i]*pi.est))+
                              n[i]*sum((A[i]*pi.est-1)*log.x.bar[i,]))
    ID <- group.id[i]  ## Added July 18, 2022
    rows.to.pull <- data.set[,1]==ID  ##Changed July 18, 2022
    group.dat <- data.set[rows.to.pull, 2:(k+1)]
    temp2 <- sirt::dirichlet.mle(group.dat)
    alt.A <- temp2$alpha0
    alt.pi <- temp2$xsi
    alt.max <- alt.max + (n[i]*lgamma(alt.A) - n[i]*sum(lgamma(alt.A*alt.pi)) + 
                            n[i]*sum((alt.A*alt.pi-1)*log.x.bar[i,]))
  }
  test.stat <- 2*(alt.max-null.max)
  return(test.stat)
}

#### Perform Test #### 
## July 17, 2022:  Perform test using 5 layer nesting tree. 
## Step 1:  Get rid of non-essential columns.

baseb.dat <- batter.dat[,4:9]

## Step 2:  Fix zero entries.

eps <- 10^-5
baseb.dat <- (baseb.dat + eps)/(1+2*eps)
baseb.dat <- baseb.dat/rowSums(baseb.dat)

##Step 3:  Reattach to groups.
baseb.dat <- cbind(Group = batter.dat[,3], baseb.dat)

## Step 4: create Layers
GM <- baseb.dat[,1]
S1 <- rowSums(baseb.dat[, 2:4])
S2 <- rowSums(baseb.dat[,5:7])
S3 <- rowSums(baseb.dat[,2:3])
S4 <- rowSums(baseb.dat[,c(5, 7)])

L1.dat <- cbind(GM, S1, S2)
a1 <- my.lrt.fun.3(L1.dat)

L2.dat <- cbind(GM, S3/S1, baseb.dat[,4]/S1)
a2 <- my.lrt.fun.3(L2.dat)

L3.dat <- cbind(GM, baseb.dat[,6]/S2, S4/S2)
a3 <- my.lrt.fun.3(L3.dat)

L4.dat <- cbind(GM, baseb.dat[,2]/S3, baseb.dat[,3]/S3)
a4 <- my.lrt.fun.3(L4.dat)

L5.dat <- cbind(GM, baseb.dat[,5]/S4, baseb.dat[,7]/S4)
a5 <- my.lrt.fun.3(L5.dat)

overall.stat <- sum(c(a1, a2, a3, a4, a5))

## The df is 3(10) - 3(5) = 15 
## p-value

pchisq(overall.stat, 10, lower.tail=FALSE)

#### All pairwise tests ####

## Groups 1 and 2
L1.dat.1 <- L1.dat[L1.dat[,1]!=3,]
L2.dat.1 <- L2.dat[L2.dat[,1]!=3 ,]
L3.dat.1 <- L3.dat[L3.dat[,1]!=3 ,]
L4.dat.1 <- L4.dat[L4.dat[,1]!=3 ,]
L5.dat.1 <- L5.dat[L5.dat[,1]!=3 ,]

d1 <- my.lrt.fun.3(L1.dat.1)
e1 <- my.lrt.fun.3(L2.dat.1)
f1 <- my.lrt.fun.3(L3.dat.1)
g1 <- my.lrt.fun.3(L4.dat.1)
h1 <- my.lrt.fun.3(L5.dat.1)
tot1 <- d1+e1+f1+g1+h1

## Groups 1 and 3
L1.dat.2 <- L1.dat[L1.dat[,1]!=2,]
L2.dat.2 <- L2.dat[L2.dat[,1]!=2 ,]
L3.dat.2 <- L3.dat[L3.dat[,1]!=2 ,]
L4.dat.2 <- L4.dat[L4.dat[,1]!=2 ,]
L5.dat.2 <- L5.dat[L5.dat[,1]!=2 ,]

d2 <- my.lrt.fun.3(L1.dat.2)
e2 <- my.lrt.fun.3(L2.dat.2)
f2 <- my.lrt.fun.3(L3.dat.2)
g2 <- my.lrt.fun.3(L4.dat.2)
h2 <- my.lrt.fun.3(L5.dat.2)
tot2 <- d2+e2+f2+g2+h2


## Groups 2 and 3
L1.dat.3 <- L1.dat[L1.dat[,1]!=1,]
L2.dat.3 <- L2.dat[L2.dat[,1]!=1 ,]
L3.dat.3 <- L3.dat[L3.dat[,1]!=1 ,]
L4.dat.3 <- L4.dat[L4.dat[,1]!=1 ,]
L5.dat.3 <- L5.dat[L5.dat[,1]!=1 ,]

d3 <- my.lrt.fun.3(L1.dat.3)
e3 <- my.lrt.fun.3(L2.dat.3)
f3 <- my.lrt.fun.3(L3.dat.3)
g3 <- my.lrt.fun.3(L4.dat.3)
h3 <- my.lrt.fun.3(L5.dat.3)
tot3 <- d3+e3+f3+g3+h3
tot3

pchisq(tot1, 5, lower.tail=FALSE)
pchisq(tot2, 5, lower.tail=FALSE)
pchisq(tot3, 5, lower.tail=FALSE)


#### July 30th, 2022:  Single Layer analysis ####

batter.dat <- read.csv("batters_grouped.csv")[,2:11]
batter.dat <- batter.dat %>% rename( "TRIPLE" = X3B,
                                     "DOUBLE" = X2B,
                                     "SINGLE" = X1B)
View(batter.dat)

baseb.dat <- batter.dat[,4:9]
eps <- 10^-5
baseb.dat <- (baseb.dat + eps)/(1+2*eps)
baseb.dat <- baseb.dat/rowSums(baseb.dat)
baseb.dat <- cbind(Group = batter.dat[,3], baseb.dat)
View(baseb.dat)

overall.stat <- my.lrt.fun.3(baseb.dat)  ## Overall Test
pchisq(overall.stat, 10, lower.tail=FALSE)

## Two by Two tests

## Groups 1 and 2
baseb.dat.1 <- baseb.dat[baseb.dat[,1]!=3,]
my.lrt.fun.3(baseb.dat.1)

## Groups 1 and 3
baseb.dat.2 <- baseb.dat[baseb.dat[,1]!=2,]
aa <- my.lrt.fun.3(baseb.dat.2)
pchisq(aa, 5, lower.tail=FALSE)

## Groups 2 and 3
baseb.dat.3 <- baseb.dat[baseb.dat[,1]!=1,]
my.lrt.fun.3(baseb.dat.3)


#### July 31, 2022:  Collapsing Tree test. ####
batter.dat <- read.csv("batters_grouped.csv")[,2:11]
batter.dat <- batter.dat %>% rename( "TRI" = X3B,
                                     "DBL" = X2B,
                                     "SGL" = X1B,
                                     "OTR" = OTHER)
baseb.dat <- batter.dat[,4:9]
eps <- 10^-5
baseb.dat <- (baseb.dat + eps)/(1+2*eps)
baseb.dat <- baseb.dat/rowSums(baseb.dat)

S1 <- rowSums(baseb.dat[, 1:3])
S2 <- rowSums(baseb.dat[,4:6])
S3 <- rowSums(baseb.dat[,1:2])
S4 <- rowSums(baseb.dat[,c(4, 6)])

L1.dat <- cbind(S1, S2)
L2.dat <- cbind(S3/S1, baseb.dat[,3]/S1)
L3.dat <- cbind(baseb.dat[,5]/S2, S4/S2)
L4.dat <- cbind(baseb.dat[,1]/S3, baseb.dat[,2]/S3)
L5.dat <- cbind(baseb.dat[,4]/S4, baseb.dat[,6]/S4)

##Function that gives a confidence interval for a4-(a1+a2) given two data sets.
alpha.ci.fun <- function(data.set.1, data.set.2, conf=0.95){
  nn <- nrow(data.set.1)
  ones <- matrix(rep(1, 4), nrow=2, ncol=2)
  alpha.1 <- sirt::dirichlet.mle(data.set.1)$alpha
  alpha.2 <- sirt::dirichlet.mle(data.set.2)$alpha
  A1 <- sum(alpha.1)
  A2 <- sum(alpha.2)
  Q1.inv <- diag(1/(-nn*trigamma(alpha.1)))
  Q2.inv <- diag(1/(-nn*trigamma(alpha.2)))
  z1 <- nn*trigamma(A1)
  z2 <- nn*trigamma(A2)
  denom.1 <- 1/z1 + sum(Q1.inv)
  denom.2 <- 1/z2 + sum(Q2.inv)
  info.invs.1 <- 1/denom.1*Q1.inv%*%ones%*%Q1.inv - Q1.inv
  info.invs.2 <- 1/denom.2*Q2.inv%*%ones%*%Q2.inv - Q2.inv
  point.est <- alpha.1[1] - A2
  var.est <- info.invs.1[1, 1] + sum(info.invs.2)
  SE <- -qnorm((1-conf)/2)*var.est
  CI <- c(point.est-SE, point.est+SE)
  return(c(CI, alpha.1, alpha.2))
  }

##Subtree 1
alpha.ci.fun(L1.dat, L2.dat)  ##CI: (-0.887, -0.886)

##Subtree 2
L1.dat.switch <- cbind(S2, S1)
alpha.ci.fun(L1.dat.switch, L3.dat)

##Subtree 3
alpha.ci.fun(L2.dat, L4.dat)

##subtree 4
L3.dat.switch <- cbind(S4/S2, baseb.dat[,5]/S2)
alpha.ci.fun(L3.dat.switch, L5.dat)


## August 8, 2022
## Examine mean vectors between three groups.

baseb.1 <- baseb.dat[baseb.dat[,1]==1,]
baseb.2 <- baseb.dat[baseb.dat[,1]==2,]
baseb.3 <- baseb.dat[baseb.dat[,1]==3,]

Group_1 <- colMeans(baseb.1)[2:7]
Group_2 <- colMeans(baseb.2)[2:7]
Group_3 <- colMeans(baseb.3)[2:7]

temp <- round(cbind(Group_1, Group_2, Group_3), 4)
xtable(temp, digits=4)
