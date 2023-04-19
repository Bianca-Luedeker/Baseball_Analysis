##Bianca Luedeker
## July 27, 2021 - August 3, 2021
## My copy of Turner's tree finding code

####Packages####
library(grDevices)
library(gtools)
library(diagram)
library(sirt)  ##Dirichlet MLE function

#####  Playing Around ####

water.maze.data <- read.csv("water_maze_data.csv")[, 2:5]
 View(water.maze.data)


n <- 1000
data.set.1 <- gtools::rdirichlet(n, c(10, 30))
data.set.2 <- gtools::rdirichlet(n, c(4, 4, 4))
data.set.A <- cbind(data.set.1[,1], data.set.1[,2]*data.set.2)  ##simulated data


##July 29: My nested data set
data.level.1 <- gtools::rdirichlet(n, c(100, 300))
data.2.left <- gtools::rdirichlet(n, c(3, 7))
data.2.right <- gtools::rdirichlet(n, c(6,4))
data.set.B <- cbind(data.level.1[,1]*data.2.left[,1],
                    data.level.1[,1]*data.2.left[,2],
                    data.level.1[,2]*data.2.right[,1],
                    data.level.1[,2]*data.2.right[,2]
                    )


n <- nrow(data.set.A)
k <- n_var <- ncol(data.set.A)
j <- floor(k/2)

## create a null vector to store all possible combinations.
combin.vec <- c()  
for(i in 1:j){
  combin.vec[i] <- choose(k, i)
}


## July 27 note: dirichletmlealphas can be replaced by sirt::dirichlet.mle$alpha

#dirichletmlealphas(data.set, ind=1, N=40)  ## Just comment this bad boy out.
sirt::dirichlet.mle(data.set.A)


#### Log Likelihood FUnctions ####

##July 27: What is loglike3?
## July28: A loglikelihood function based on the mean and precision instead of 
## individual alpha.  This is on Turner page 44.
## There is a mistake in Turner's loglikelihoods as he is missing a -1


## Either of these loglikelihoods should work.

##Input is the alphas.
log.like.1 <- function(alpha, data){
  n <- nrow(data)
  log.p.bar <- apply(log(data), 2, mean)
  #result <- n*lgamma(sum(alpha)) - n*sum(lgamma(alpha)) + n*sum((alpha-1)*log.p.bar)  ##INCorrect version
  result <- n*lgamma(sum(alpha)) - n*sum(lgamma(alpha)) + n*sum(alpha*log.p.bar)  ##Correct version
  return(result)
}

### August 3rd Update.  So in Turner's code, he has alpha and not alpha - 1.
##  I claim mine is correct, but his gives results that match the simulation 
## settings.  
## June 19, 2022:  As you can see. The loglikelihood function has been fixed to 
## Turner's method.

## Input is A the precision (sum of the alphas), mu is the mean vector (alphas/A)
log.like.2 <- function(mu, A, data){
  n <- nrow(data)
  log.p.bar <- apply(log(data), 2, mean)
  result <- n*lgamma(A) - n*sum(lgamma(A*mu)) + n*sum((A*mu-1)*log.p.bar)
  return(result)
}

#### Binary Tree Algorithm ####

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

#### August 3rd 2021 ####

## Comparison of my tree.finder program and Turner's Split Tree Function.

## Simulated data sets

n <- 1000

## Data Set A:
data.set.1 <- gtools::rdirichlet(n, c(10, 30))
data.set.2 <- gtools::rdirichlet(n, c(4, 4, 4))
data.set.A <- cbind(data.set.1[,1], data.set.1[,2]*data.set.2)  ##simulated data
head(data.set.A)

##Nesting Structure of data set A should be {(X1), (X2, X3, X4)}


split.tree(1:4, data.set.A)  ## Gives correct results.
tree.finder(1:4, data.set.A)

## Data Set B:

data.level.1 <- gtools::rdirichlet(n, c(100, 300))
data.2.left <- gtools::rdirichlet(n, c(3, 7))
data.2.right <- gtools::rdirichlet(n, c(6,4))
data.set.B <- cbind(data.level.1[,1]*data.2.left[,1],
                    data.level.1[,1]*data.2.left[,2],
                    data.level.1[,2]*data.2.right[,1],
                    data.level.1[,2]*data.2.right[,2])
head(data.set.B)

## For data set B, we should have {(X1, X2), (X3, X4)}

split.tree(1:4, data.set.B)  ## Gives correct results.
tree.finder(1:4, data.set.B)

split.tree(1:4, data.set.B)$criteria
tree.finder(1:4, data.set.B)$criteria

## Data set C:

data.level.1 <- gtools::rdirichlet(n, c(2, 3))
data.2.left <- gtools::rdirichlet(n, c(3, 7))
data.2.right <- gtools::rdirichlet(n, c(2, 4,6))
data.set.C <- cbind(data.level.1[,1]*data.2.left[,1],
                    data.level.1[,1]*data.2.left[,2],
                    data.level.1[,2]*data.2.right[,1],
                    data.level.1[,2]*data.2.right[,2],
                    data.level.1[,2]*data.2.right[,3])

data.set.C <- data.set.C[, c(3, 1, 4, 2, 5)]  ## scramble the columns

## For data set C, we should have {(X2, X4), (X1, X3, X5)}
head(data.set.C) 

split.tree(1:5, data.set.C)  ## Gives correct results.
tree.finder(1:5, data.set.C)


## Data set D:
data.level.1 <- gtools::rdirichlet(n, c(4, 5))
data.2.left <- gtools::rdirichlet(n, c(1, 2, 3))
data.2.right <- gtools::rdirichlet(n, c(2, 2,6))
data.set.D <- cbind(data.level.1[,1]*data.2.left[,1],
                    data.level.1[,1]*data.2.left[,2],
                    data.level.1[,1]*data.2.left[,3],
                    data.level.1[,2]*data.2.right[,1],
                    data.level.1[,2]*data.2.right[,2],
                    data.level.1[,2]*data.2.right[,3])

data.set.D <- data.set.D[, c(1, 4, 2, 5, 3, 6)]  ## scramble the columns 
## For data set C, we should have {(X1, X3, X5), (X2, X4, X6)}

head(data.set.D) 

split.tree(1:6, data.set.D)  ## Gives correct results.
tree.finder(1:6, data.set.D)

#### Find the best tree with the water maze data ####
water.maze.data <- read.csv("water_maze_data.csv", header = TRUE)[, 2:5]
head(water.maze.data)

tree.finder(1:4, water.maze.data)  ## Groups {(TQ, AQ2), (AQ1, OQ)}  {(1, 4), (2,3)}

## Level 1 MLES
left.totals <- rowSums(water.maze.data[,c(1,4)])
right.totals <- rowSums(water.maze.data[,c(2,3)])
level.1 <- cbind(left.totals, right.totals)
level.2.left  <- water.maze.data[,c(1,4)]/left.totals
level.2.right <- water.maze.data[,c(2,3)]/right.totals

sirt::dirichlet.mle(level.1)$alpha
sirt::dirichlet.mle(level.2.left)$alpha
sirt::dirichlet.mle(level.2.right)$alpha

## Does simulated data look like the water maze data?
n <- 1000
data.level.1 <- gtools::rdirichlet(n, c(11.2, 8.1))
data.2.left <- gtools::rdirichlet(n, c(9.2, 5.6))
data.2.right <- gtools::rdirichlet(n, c(11.6, 10.3))
data.set.WM <- cbind(data.level.1[,1]*data.2.left[,1],
                    data.level.1[,1]*data.2.left[,2],
                    data.level.1[,2]*data.2.right[,1],
                    data.level.1[,2]*data.2.right[,2])

data.set.WM <- data.set.WM[, c(1, 3, 4, 2)]  ##Reorganize the columns to match orignal data set.

cor(data.set.WM)
cor(water.maze.data)

#### August 4th 2021 ####
## Recap:  Both Turner's split.tree and my tree.finder work.  These functions
## only do one binary split.  The input is the data set, the output is 
## the best fitting single binary split.  To get the best fitting tree, 
## this needs to be applied to the data set over and over.

##  Next Function:  directions.next What does this function do?

directions.next <- function( x, direction.mat){
  cbind( rep(direction.mat[,2], each = 2), 2^x:(2^(x+1)-1) )  ##each column has 2^x elements
  ## directions.mat must have nrows = 2^(x-1)
  ## Not certain about what this does yet.
}



#### Complete Binary tree function ####
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


##Update: I think I can get the generalized binary tree for any example now.
## My test example is this
## First Split { (1, 2, 3, 4), (5, 6, 7)}; alpha = (30,70)
## Left Split 1: {(1), (2, 3, 4)}; alpha = (2, 4)
## Left Split 2: {(2), (3, 4)}; alpha = (1, 2),  final alpha = (3,5)
## Right Split 1: {(5), (6,7)}; alpha = (8, 6); final alpha = (1, 2)

## Simulated data set E
set.seed(1983)
n <- 1000
level.1 <- gtools::rdirichlet(n, c(30,70))
left.1 <- gtools::rdirichlet(n, c(2, 4))
left.2 <- gtools::rdirichlet(n, c(1,2))
left.3 <- gtools::rdirichlet(n, c(3,5))
right.1 <- gtools::rdirichlet(n, c(8,6))
right.2 <- gtools::rdirichlet(n, c(1,2))

data.set.E <- cbind(
  V1 = level.1[,1]*left.1[,1],
  V2 = level.1[,1]*left.1[,2]*left.2[,1],
  V3 = level.1[,1]*left.1[,2]*left.2[,2]*left.3[,1],
  V4 = level.1[,1]*left.1[,2]*left.2[,2]*left.3[,2],
  V5 = level.1[,2]*right.1[,1],
  V6 = level.1[,2]*right.1[,2]*right.2[,1],
  V7 = level.1[,2]*right.1[,2]*right.2[,2]
)

head(data.set.E)
head(rowSums(data.set.E))


complete.bin.tree(data.set.E)

##Note:  Using the current set of parameters (from Aug 4) the algorithm gets the
## nesting correct.  However, I can alter the parameters to cause the algorithm 
## to fail.

complete.bin.tree(data.set.D)




## August 5th 
## Pruning attempt using confidence intervals.
##  The smallest example is a tree with three nodes that can be collasped.
## The generalized tree will be {(V1), (V2, V3)}


## Simulation


## Check how variance works with toy example.
## Create m = 1000 data sets each of size n = 100

m <- 1000 ##number of simulations
n <- 100 ## sample size
