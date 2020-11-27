library(GBClust)
library(ggplot2)
library(mcclust)
library(mvtnorm)

rm(list = ls())

# ----------------------------
# Oracle Similarity matrix
# ----------------------------

Miscl <- function(S, cluster, medoids) {
  pr_miscl <- numeric(nrow(S))
  for (i in 1:n) {
    Gi <- cluster[i]
    pr_miscl[i] <- 1 - S[i, medoids[cluster[i]]]
  }
  pr_miscl
}

# -----------------
# Initialization
# -----------------

d <- 2
H <- 2
n <- 200

sigma2 <- 1

G0 <- rep(1:H, each = n / H)
mu <- matrix(c(-1.2, 1.2, -1.2, 1.2), 2, 2)

# DATASET CREATION --------------------

set.seed(1234)
dataset <- NULL
for (h in 1:H) {
  dataset <- rbind(dataset, rmvnorm(n = n / H, mean = mu[h, ], sigma = sigma2 * diag(d)))
}

# ESTIMATION --------------------------

D <- as.matrix(dist(dataset, method = "euclidean"))^2
fit <- kmeans2(dataset, k = H, nstart = 10)
fit_gibbs <- kmeans_gibbs(dataset, k = H, a_lambda = 0, b_lambda = 0, R = 5000, burn_in = 1000, trace = TRUE)

SGibbs <- mcclust::comp.psm(fit_gibbs$G)
MisclGibbs_pr <- Miscl(SGibbs, fit$cluster, comp_medoids(D, fit$cluster))

# Plot 1

p1 <- ggplot(data = data.frame(dataset, Cluster = as.factor(fit$cluster)), aes(x = X1, y = X2, col = Cluster, shape = Cluster)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  scale_color_brewer(palette = "Set1") +
  ggtitle("K-means clustering")

# Plot 2
data_plot <- rbind(data.frame(dataset, Probability = MisclGibbs_pr))

p2 <- ggplot(data = data_plot, aes(x = X1, y = X2, col = Probability)) +
  geom_point(size = 2) +
  theme_light() +
  xlab("") +
  ylab("") +
  scale_color_distiller(palette = "Spectral", direction = -1, limits = c(0, .75)) +
  ggtitle("Misclassification probabilities")


# library(gridExtra)
#
# p3 <- arrangeGrob(p1, p2, ncol = 2)
#
# ggsave("intro.pdf", p3, width = 10, height = 4.2)
