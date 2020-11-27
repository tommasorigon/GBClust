library(GBClust)
library(ggplot2)
library(mcclust)
library(mcclust.ext)
library(mvtnorm)
library(reshape2)
library(gridExtra)

rm(list = ls())

# Oracle Similarity matrix ---------------------

ProbsOracle <- function(x, mu, sigma2) {
  out <- apply(mu, 1, function(means) dmvnorm(x, means, sigma2 * diag(ncol(x))))
  out / rowSums(out)
}

SOracle <- function(x, mu, sigma2) {
  probs <- ProbsOracle(x, mu, sigma2)
  out <- matrix(0, nrow(x), nrow(x))
  for (h in 1:ncol(probs)) {
    out <- out + probs[, h] %*% t(probs[, h])
  }
  diag(out) <- rep(1, nrow(x))
  out
}

Miscl <- function(S, cluster, medoids) {
  pr_miscl <- numeric(nrow(S))
  for (i in 1:n) {
    Gi <- cluster[i]
    pr_miscl[i] <- 1 - S[i, medoids[cluster[i]]]
  }
  pr_miscl
}

# Initialization ----------------------

d <- 2
H <- 4
sigma2 <- numeric(3)
n <- 200

# Misclassification
MisclG <- numeric(3)
MisclO <- numeric(3)

# Co-clustering
MAE_S <- numeric(3)
RMSE_S <- numeric(3)

# Variation of information
VIG <- numeric(3)
VIO <- numeric(3)
VIHB <- numeric(3)

# Scenario 1 -----------------------------

sigma2[1] <- 0.75

G0 <- rep(1:H, each = n / H)
mu <- as.matrix(expand.grid(seq(-2, 2, length = 2), seq(-2, 2, length = 2)))

# Dataset creation --------------------

set.seed(123)
dataset <- NULL
for (h in 1:H) {
  dataset <- rbind(dataset, rmvnorm(n = n / H, mean = mu[h, ], sigma = sigma2[1] * diag(d)))
}

p1_g0 <- ggplot(data = data.frame(dataset, Cluster = as.factor(G0)), aes(x = Var1, y = Var2, col = Cluster)) +
  geom_point() +
  theme_light() +
  xlab("") +
  ylab("")

# Estimation --------------------------

D <- as.matrix(dist(dataset, method = "euclidean"))^2

p1_select <- kdiss_select(D, k_max = 10, method = "silhouette", nstart = 5)
fit <- kdiss(D, H, nstart = 10)

fit_gibbs <- kmeans_gibbs(dataset, k = H, a_lambda = 0, b_lambda = 0, R = 5000, burn_in = 1000, trace = TRUE)

# Quantities of interest-----------------------

SGibbs <- mcclust::comp.psm(fit_gibbs$G)
MisclGibbs_pr <- Miscl(SGibbs, fit$cluster, comp_medoids(D, fit$cluster))

ProbsO <- ProbsOracle(dataset, mu, sigma2[1])
SO <- SOracle(dataset, mu, sigma2[1])
MisclO_pr <- Miscl(SO, apply(ProbsO, 1, which.max), comp_medoids(D, apply(ProbsO, 1, which.max))) 

# Misclassification ------------------------------

MisclG[1] <- mean(MisclGibbs_pr)
MisclO[1] <- mean(MisclO_pr)

data_plot <- rbind(
  data.frame(dataset, Misclassification = MisclO_pr, Type = "Oracle"),
  data.frame(dataset, Misclassification = MisclGibbs_pr, Type = "Squared Euclidean")
)

p1_miscl <- ggplot(data = data_plot, aes(x = Var1, y = Var2, col = Misclassification)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) +
  facet_grid(. ~ Type)

# Co-clustering ----------------------------------

MAE_S[1] <- mean(abs(c(SO[lower.tri(SO)] - SGibbs[lower.tri(SGibbs)])))
RMSE_S[1] <- sqrt(mean((c(SO[lower.tri(SO)] - SGibbs[lower.tri(SGibbs)]))^2))

data_plot <- data.frame(melt(SO), Type = "Oracle")
data_plot <- rbind(data_plot, data.frame(melt(SGibbs), Type = "Squared Euclidean"))

colnames(data_plot) <- c("Var1", "Var2", "Probability", "Type")

p1_co <- ggplot(data = data_plot, aes(x = Var1, y = Var2, fill = Probability)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  xlab("") +
  ylab("") +
  theme_light() +
  facet_grid(. ~ Type)

# Horizontal Bounds -------------------------------

Credible <- credibleball(t(fit$cluster), fit_gibbs$G, alpha = 0.05, c.dist = "VI")

data_plot <- rbind(
  data.frame(dataset, Cluster = as.factor(Credible$c.star), Type = "Maximum a posteriori"),
  data.frame(dataset, Cluster = as.factor(Credible$c.horiz), Type = "Horizontal bound")
)

p1_hb <- ggplot(data = data_plot, aes(x = Var1, y = Var2, col = Cluster)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  facet_grid(. ~ Type)

# Variation of Information ------------------------

VIG[1] <- vi.dist(G0, fit$cluster, base = 2)
VIO[1] <- vi.dist(G0, apply(ProbsO, 1, which.max), base = 2)
VIHB[1] <- vi.dist(fit$cluster, Credible$c.horiz)

# Scenario 2 --------------------------------------

sigma2[2] <- 1.5

G0 <- rep(1:H, each = n / H)
mu <- as.matrix(expand.grid(seq(-2, 2, length = 2), seq(-2, 2, length = 2)))

# Dataset creation --------------------

set.seed(123)
dataset <- NULL
for (h in 1:H) {
  dataset <- rbind(dataset, rmvnorm(n = n / H, mean = mu[h, ], sigma = sigma2[2] * diag(d)))
}

p2_g0 <- ggplot(data = data.frame(dataset, Cluster = as.factor(G0)), aes(x = Var1, y = Var2)) +
  geom_point() +
  theme_light() +
  xlab("") +
  ylab("")

# Estimation --------------------------

D <- as.matrix(dist(dataset, method = "euclidean"))^2

p2_select <- kdiss_select(D, k_max = 10, method = "silhouette", nstart = 5)
fit <- kdiss(D, H, nstart = 10)

fit_gibbs <- kmeans_gibbs(dataset, k = H, a_lambda = 0, b_lambda = 0, R = 5000, burn_in = 1000, trace = TRUE)

# Quantities of interests -----------------------

SGibbs <- mcclust::comp.psm(fit_gibbs$G)
MisclGibbs_pr <- Miscl(SGibbs, fit$cluster, comp_medoids(D, fit$cluster))

ProbsO <- ProbsOracle(dataset, mu, sigma2[2])
SO <- SOracle(dataset, mu, sigma2[2])
MisclO_pr <- Miscl(SO, apply(ProbsO, 1, which.max), comp_medoids(D, apply(ProbsO, 1, which.max))) 

# Misclassification ------------------------------

MisclG[2] <- mean(MisclGibbs_pr)
MisclO[2] <- mean(MisclO_pr)

data_plot <- rbind(
  data.frame(dataset, Misclassification = MisclO_pr, Type = "Oracle"),
  data.frame(dataset, Misclassification = MisclGibbs_pr, Type = "Squared Euclidean")
)

p2_miscl <- ggplot(data = data_plot, aes(x = Var1, y = Var2, col = Misclassification)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) +
  facet_grid(. ~ Type)

# Co-clustering ----------------------------------

MAE_S[2] <- mean(abs(c(SO[lower.tri(SO)] - SGibbs[lower.tri(SGibbs)])))
RMSE_S[2] <- sqrt(mean((c(SO[lower.tri(SO)] - SGibbs[lower.tri(SGibbs)]))^2))

data_plot <- data.frame(melt(SO), Type = "Oracle")
data_plot <- rbind(data_plot, data.frame(melt(SGibbs), Type = "Squared Euclidean"))

colnames(data_plot) <- c("Var1", "Var2", "Probability", "Type")

p2_co <- ggplot(data = data_plot, aes(x = Var1, y = Var2, fill = Probability)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  xlab("") +
  ylab("") +
  theme_light() +
  facet_grid(. ~ Type)

# Horizontal Bounds -------------------------------

Credible <- credibleball(t(fit$cluster), fit_gibbs$G, alpha = 0.05, c.dist = "VI")

Credible$c.star <- as.factor(Credible$c.star)
Credible$c.horiz <- as.factor(Credible$c.horiz)

levels(Credible$c.horiz)[c(3, 1, 4, 2)] <- c(1, 2, 3, 4)

data_plot <- rbind(
  data.frame(dataset, Cluster = Credible$c.star, Type = "Maximum a posteriori"),
  data.frame(dataset, Cluster = Credible$c.horiz, Type = "Horizontal bound")
)

p2_hb <- ggplot(data = data_plot, aes(x = Var1, y = Var2, col = Cluster, shape = Cluster)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  facet_grid(. ~ Type)

# Variation of Information ------------------------

VIG[2] <- vi.dist(G0, fit$cluster, base = 2)
VIO[2] <- vi.dist(G0, apply(ProbsO, 1, which.max), base = 2)
VIHB[2] <- vi.dist(fit$cluster, Credible$c.horiz)

# Scenario 3 ------------------------------

sigma2[3] <- 3

G0 <- rep(1:H, each = n / H)
mu <- as.matrix(expand.grid(seq(-2, 2, length = 2), seq(-2, 2, length = 2)))

# Dataset creation --------------------

set.seed(123)
dataset <- NULL
for (h in 1:H) {
  dataset <- rbind(dataset, rmvnorm(n = n / H, mean = mu[h, ], sigma = sigma2[3] * diag(d)))
}

p3_g0 <- ggplot(data = data.frame(dataset, Cluster = as.factor(G0)), aes(x = Var1, y = Var2, col = Cluster)) +
  geom_point() +
  theme_light() +
  xlab("") +
  ylab("")

# Estimation --------------------------

D <- as.matrix(dist(dataset, method = "euclidean"))^2

p3_select <- kdiss_select(D, k_max = 10, method = "silhouette", nstart = 5)
fit <- kdiss(D, H, nstart = 10)

fit_gibbs <- kmeans_gibbs(dataset, k = H, a_lambda = 0, b_lambda = 0, R = 5000, burn_in = 1000, trace = TRUE)

# Quantity of interests -----------------------

SGibbs <- mcclust::comp.psm(fit_gibbs$G)
MisclGibbs_pr <- Miscl(SGibbs, fit$cluster, comp_medoids(D, fit$cluster))

ProbsO <- ProbsOracle(dataset, mu, sigma2[3])
SO <- SOracle(dataset, mu, sigma2[3])
MisclO_pr <- Miscl(SO, apply(ProbsO, 1, which.max), comp_medoids(D, apply(ProbsO, 1, which.max))) # MisclOracle(dataset,mu,sigma2[1])

# Misclassification ------------------------------

MisclG[3] <- mean(MisclGibbs_pr)
MisclO[3] <- mean(MisclO_pr)

data_plot <- rbind(
  data.frame(dataset, Misclassification = MisclO_pr, Type = "Oracle"),
  data.frame(dataset, Misclassification = MisclGibbs_pr, Type = "Squared Euclidean")
)

p3_miscl <- ggplot(data = data_plot, aes(x = Var1, y = Var2, col = Misclassification)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) +
  facet_grid(. ~ Type)

# Co-clustering ----------------------------------

MAE_S[3] <- mean(abs(c(SO[lower.tri(SO)] - SGibbs[lower.tri(SGibbs)])))
RMSE_S[3] <- sqrt(mean((c(SO[lower.tri(SO)] - SGibbs[lower.tri(SGibbs)]))^2))

data_plot <- data.frame(melt(SO), Type = "Oracle")
data_plot <- rbind(data_plot, data.frame(melt(SGibbs), Type = "Squared Euclidean"))

colnames(data_plot) <- c("Var1", "Var2", "Probability", "Type")

p3_co <- ggplot(data = data_plot, aes(x = Var1, y = Var2, fill = Probability)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  xlab("") +
  ylab("") +
  theme_light() +
  facet_grid(. ~ Type)

# Horizontal Bounds -------------------------------

Credible <- credibleball(t(fit$cluster), fit_gibbs$G, alpha = 0.05, c.dist = "VI")

data_plot <- rbind(
  data.frame(dataset, Cluster = as.factor(Credible$c.star), Type = "Maximum a posteriori"),
  data.frame(dataset, Cluster = as.factor(Credible$c.horiz), Type = "Horizontal bound")
)

p3_hb <- ggplot(data = data_plot, aes(x = Var1, y = Var2, col = Cluster)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  facet_grid(. ~ Type)

# Variation of Information ------------------------

VIG[3] <- vi.dist(G0, fit$cluster, base = 2)
VIO[3] <- vi.dist(G0, apply(ProbsO, 1, which.max), base = 2)
VIHB[3] <- vi.dist(fit$cluster, Credible$c.horiz)

# Silhouette
grid.arrange(p1_select, p2_select, p3_select, ncol = 1)

# Co-clustering
grid.arrange(p1_co, p2_co, p3_co, ncol = 1)

# Co-clustering
grid.arrange(p1_miscl, p2_miscl, p3_miscl, ncol = 1)

# Credible balls
grid.arrange(p1_hb, p2_hb, p3_hb, ncol = 1)

# Scenario 1
grid.arrange(p1_co, p1_miscl, p1_hb, ncol = 1)

# Save the graphs

ggsave("Sim1_hb.pdf", p2_hb, width = 9, height = 4)
ggsave("Sim1_Co.pdf", p2_co, width = 9, height = 4)
ggsave("Sim1_miscl.pdf", p2_miscl, width = 9, height = 4)

library(knitr)

kable(data.frame(sigma2, MAE_S, VIHB, VIG, VIO), digits = 4)

library(xtable)
xtable(data.frame(sigma2, MAE_S, VIHB, VIG, VIO), digits = 4)
