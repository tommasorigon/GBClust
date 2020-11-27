library(GBClust)
library(ggplot2)
library(mcclust)
library(mvtnorm)
library(reshape2)
library(gridExtra)
library(sn)

rm(list=ls())

# ----------------------------
# Oracle Similarity matrix
# ----------------------------

ProbsOracle <- function(x, xi, Omega, alpha, nu){
  out <- apply(xi, 1, function(means) dmst(x, xi = means, Omega = Omega, alpha = alpha, nu = nu))
  out / rowSums(out)
}

SOracle <- function(x, xi, Omega, alpha, nu){
  probs <- ProbsOracle(x, xi, Omega, alpha, nu)
  out <- matrix(0,nrow(x),nrow(x))
  for(h in 1:ncol(probs)){
    out <- out + probs[,h]%*%t(probs[,h])
  }
  diag(out) <- rep(1,nrow(x))
  out
}

MisclOracle <- function(x, xi, Omega, alpha, nu){
  probs      <- ProbsOracle(x, xi, Omega, alpha, nu)
  pr_miscl   <- 1 - apply(probs,1, max)
  pr_miscl
}

Miscl <- function(S, cluster, medoids){
  
  pr_miscl <- numeric(nrow(S))
  for(i in 1:n){
    Gi <- cluster[i]
    pr_miscl[i] <- 1 - S[i,medoids[cluster[i]]]
  }
  pr_miscl
}

# -----------------
# Initialization
# -----------------

d <- 2
H <- 4
sigma2 <- numeric(3)
n <- 200

# -----------------------------------------
# SCENARIO 1 - GRAPHS
# -----------------------------------------

G0 <- rep(1:H, each = n/H)

xi    <- as.matrix(expand.grid(seq(-2,2,length=2), seq(-2,2,length=2)))
Omega <- matrix(c(1,0,0,1),2,2)
alpha <- c(0,0)
nu    <- 2

# DATASET CREATION --------------------
set.seed(123456)
dataset <- NULL
for(h in 1:H){
  dataset <- rbind(dataset, rmst(n = n/H, xi = xi[h,], Omega = Omega, alpha = alpha, nu = nu))
}

ggplot(data=data.frame(dataset,Cluster=as.factor(G0)), aes(x=X1, y=X2,col=Cluster)) + geom_point() + theme_light() + xlab("") + ylab("")

# ESTIMATION --------------------------
D_kmeans  <- as.matrix(dist(dataset,method = "euclidean"))^2
D_diss    <- as.matrix(dist(dataset,method = "manhattan"))

p1_select_kmeans <- kdiss_select(D_kmeans, k_max = 15,method = "silhouette", nstart = 50) + ggtitle("Squared Euclidean")
fit_kmeans <- kdiss(D_kmeans, H, nstart = 50)
p1_select_kdiss <- kdiss_select(D_diss, k_max = 15, method = "silhouette", nstart = 50) + ggtitle("Manhattan")
fit_kdiss   <- kdiss(D_diss, H, nstart = 50)

grid.arrange(p1_select_kmeans,p1_select_kdiss,ncol=2)
psim2_kdiss <- arrangeGrob(p1_select_kmeans,p1_select_kdiss,ncol=2)


# # Cluster visualization
ggplot(data=data.frame(dataset,Cluster=as.factor(fit_kmeans$cluster)), aes(x=X1,y=X2,col=Cluster)) + geom_point() + theme_light() +xlab("x") + ylab("y")

ggplot(data=data.frame(dataset,Cluster=as.factor(fit_kdiss$cluster)), aes(x=X1, y=X2,col=Cluster)) + geom_point() + theme_light() +xlab("x") + ylab("y")


# MCMC
fit_gibbs_kmeans <- kmeans_gibbs(dataset, k = H, a_lambda = 0, b_lambda = 0, R = 5000, burn_in = 1000, trace=TRUE)
fit_gibbs_kdiss <- Minkowski_gibbs(dataset, k = H, p = 1,  a_lambda = 0, b_lambda = 0, R = 5000, burn_in = 1000, trace=TRUE)

# QUANTITY OF INTERESTS -----------------------

SGibbs_kmeans        <- mcclust::comp.psm(fit_gibbs_kmeans$G)
SGibbs_kdiss         <- mcclust::comp.psm(fit_gibbs_kdiss$G)

# CO-CLUSTERING ----------------------------------

ProbsO    <- ProbsOracle(dataset, xi, Omega, alpha, nu)
SO        <- SOracle(dataset, xi, Omega, alpha, nu)

MAE_S_kmeans  <- mean(abs(c(SO[lower.tri(SO)] - SGibbs_kmeans[lower.tri(SGibbs_kmeans)])))
MAE_S_kdiss   <- mean(abs(c(SO[lower.tri(SO)] - SGibbs_kdiss[lower.tri(SGibbs_kdiss)])))

data_plot <- data.frame(melt(SO), Type="Oracle")
data_plot <- rbind(data_plot, data.frame(melt(SGibbs_kmeans), Type="Squared Euclidean"))
data_plot <- rbind(data_plot, data.frame(melt(SGibbs_kdiss), Type="Manhattan"))

colnames(data_plot) <- c("Var1","Var2","Probability","Type")

p1_co <- ggplot(data=data_plot, aes(x=Var1,y=Var2,fill=Probability)) + geom_tile()+ scale_fill_distiller(palette="Blues",direction=1) + xlab("") + ylab("") + theme_light()+ facet_grid(.~Type) + theme(legend.position = "none")


# -----------------------------------------------
# Save GRAPHS
# -----------------------------------------------

setwd("~/Google Drive/University/Lavori/Bayesian Gibbs Clustering/img")

ggsave("Sim2_Co.pdf",p1_co, width=9, height = 3.5)
ggsave("Sim2_select.pdf",psim2_kdiss, width=9, height = 4)
