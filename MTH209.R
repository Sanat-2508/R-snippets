####libraries to load 

library(pracma)
library(MASS)
library(matlib)
library(mvtnorm)
library(ggplot2)
library(boot)
library(factoextra)
library(ISLR)
library(ggfortify)
library(gridExtra)
library(animation)
library(jpeg)
library(dplyr)
library(gMOIP)
library(lpSolve)
library(igraph)
library(wCorr)
library(plot3D)
library(fMultivar)
library(graphics)
library(MCMCpack)


#######libraries to load ##### end here



###### Lectures o revise

#L6_1
#L12_4 (Bayesian stuff)


#### elementary transformations


rowadd(A, 1, 2, -2)
#A[2, ] -> A[2, ] - 2*A[1, ]


rowmult(A, 1, -2)
#A[1, ] -> - 2*A[1, ]



##### clusters
d <- dist(df, method = "euclidean")
hc <- hclust(d)
fviz_cluster(list(data = df, cluster = sub_grp))

####PCA
pc <- prcomp(auto, scale=TRUE)

####SMV classifictaion plot

z <- svm(default ~ ., Default[,c(1,3:4)])
plot(z, Default[,c(1,3:4)])



####clustering using kmeans
clusters <- kmeans(df, 4, nstart = 10)
fviz_cluster(clusters, df)




####Linear Programmming snippets

#### Ploting


#A <- matrix(c(1, 1, 1, 7, -1, 1, -1, 0, 0, -1), nrow = 5, ncol = 2, byrow = TRUE)
#b <- c(7, 35, 3, 0, 0)
plotPolytope(
  A,
  b,
  type = rep("c", ncol(A)),
  crit = "max",
  faces = rep("c", ncol(A)),
  plotFaces = TRUE,
  plotFeasible = TRUE,
  plotOptimum = FALSE,
  labels = "coord"
)


A <- matrix(c(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 ,
              0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ,
              1 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 ,
              0 , 0 , 0 , 0 , 0 , 1), nrow = 11, ncol = 18, byrow = TRUE)

f.con <- rbind(A, diag(18))
f.rhs <- c(1298, 1948, 465, 605, 451, 338, 260, 183, 282, 127, 535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

f.obj <- c(39.99, 126.27, 102.70, 81.68, 38.81, 71.99, 31.21, 22.28, 321.04, 145.36, 33.82, 154.05, 64.19, 87.90, 107.98, 65.45, 39.08, 167.38)
f.dir <- c("=", "=", "=", "=", "=", "=", "=", "=", "=", "=", "=", 
           ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=", 
           ">=", ">=", ">=", ">=", ">=", ">=", ">=", ">=")

Sol <- lp ("min", f.obj, f.con, f.dir, f.rhs)
Sol$objval



###### Network Ananlysis


# A = Adjacency matrix
# D = degree matrix      // a_ii = d_i = degree of vertex i
#                        // a_ij =0
# L = D - A              // Laplacian Matrix
# P = Transiton matrix   // p_ij = a_ij/d_i
# Normalised Laplacian = I - P


# sum of eigen values of normalised laplacian <= n, equality hold iif there is no isolated node 

# dimension of the nullspace of the Laplacian matrix equals the number of connected components of the corresponding graph.




####### probabilty and statistics

# multinomial ditribution
dmultinom(c(0,1,11),prob=c(0.20, 0.30, 0.50))



### QQ plot
Y1 <- rnorm(1000)
qqnorm(Y1, col='blue', main='Y1 ~ N(0,1)', xlim = c(-3,3), ylim = c(-3,3))



#### optimisation

norm.data <- rnorm(100,1,2)

normal.lik<-function(theta,y){
  mu <- theta[1]
  sigma2 <- theta[2]
  n <- length(y)
  logl<- -n/2*log(2*pi) -n/2*log(sigma2) - (2*sigma2)^(-1)*sum((y-mu)^2)
  return(-logl)
}

optim(c(0,1),normal.lik,y=norm.data,method="BFGS")



#### confidence interval

t.test(dat,conf.level=0.95)$conf.int

t.test(cogbehav,control,var.equal=TRUE,conf.level=0.95)



#### Bootstrap

boot.results<-boot(Books$P,function(x,b){median(x[b])}, 10000) #function for data x
boot.ci(boot.results) # resampling case b

##### parametric bootstrap

# ................ required functions ....................................

y.median <- function(y, i){ # y: data vector; i: indices vector
  return(median(y[i]))}
gamma.rg <- function(y,p){ # function to generate random gamma variates
  rgamma(length(y), shape=p[1], rate=p[2])}

#.........................................................................

y <- c(5.88,5.55,5.40,1.83,2.31,1.32,1.52,6.79,4.99,3.87,1.21,10.44,3.71,1.68,2.53,5.40,0.17,9.00,1.41,3.37,2.99,1.68,1.73,6.43,4.16)
s <- 2; r <- 0.5
p <- c(s,r)
gamma_bootmed = boot(y, y.median, R=10000, sim = "parametric", ran.gen=gamma.rg, mle=p)
gamma_bootmed


###### hypothesis testing

prop.test(524,1008,p=0.50,alt="two.sided",conf.level=0.95, correct=FALSE)

t.test(Polid$ideology[Polid$race=="hispanic"], mu=4.0, alt="two.sided")


#######################################################
### using bayes (no idea what it is doing)#############
###########################################################
fit<-MCMCregress(change[Anorexia$therapy=="cb"]~1,mcmc=5000000,
                 b0=0,B0=10^{-15},c0=10^{-15},d0=10^{-15})
# mean has normal prior dist. with mean b0=0, variance 1/B0 (std.dev. > 31 million)
# variance has inverse gamma prior distribution (c0/2=shape, d0/2=scale)
summary(fit)


fit.bayes<-MCMCregress(after-before ~ factor(therapy), mcmc=10000000, b0=0, B0=10^{-15}, c0=10^{-15}, d0=10^{-15}, data=Anor[-(30:46),])
# mean has normal prior dist. with mean b0=0, variance 1/B0(std.dev.>31million)
# variance has highly disperse inverse gamma prior distribution (tiny c0 and d0)
summary(fit.bayes)







