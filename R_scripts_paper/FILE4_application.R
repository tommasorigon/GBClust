# GBClust needs to be installed --------
library(GBClust)

# Other packages
library(reshape2)
library(ggplot2)
library(poLCA)

# Importing the dataset ---------------------

rm(list = ls())

data("carcinoma")
dataset <- carcinoma
dataset <- as.matrix(dataset - 1) > 0

# Selection of K ---------------------------

set.seed(123)
p_select <- kbinary_select(dataset, k_max = 8, nstart = 250)
p_select

# Point estimate ---------------------------

k <- 3
set.seed(123)
fit <- kbinary(x = dataset, k = k, trace = FALSE, nstart = 250)
G_map <- fit$cluster

fit$loss

# Gibbs sampling ---------------------------

set.seed(123)
fit_gibbs <- kbinary_gibbs(dataset, k, lambda = 1, R = 15000, burn_in = 1000, trace = TRUE)

# plot(fit_gibbs$loss,type="l")
# min(fit_gibbs$loss)

# Alternative point estimates -----------------

library(mcclust.ext)
G_VI <- minVI(mcclust::comp.psm(fit_gibbs$G))$cl

table(G_map, G_VI)

loss <- 0
for (j in 1:k) {
  loss <- loss + GBClust:::loss_binary(dataset[G_map == j, ])
}
loss

loss <- 0
m_centroids <- matrix(0, k, ncol(dataset))
for (j in 1:k) {
  nj <- sum(G_VI == j)
  m_centroids[j, ] <- nj / (nj + 1) * colMeans(dataset[G_VI == j, ]) + 0.5 / (nj + 1)
  loss <- loss + GBClust:::loss_binary(dataset[G_VI == j, ])
}
loss

# Centroids --------------------------------

colnames(m_centroids) <- colnames(dataset)
library(knitr)
library(xtable)

kable(m_centroids, digits = 2)
xtable(m_centroids, digits = 2)

# Co-clustering matrix -----------------------

SGibbs <- mcclust::comp.psm(fit_gibbs$G)
SGibbs <- SGibbs[order(G_map), order(G_map)]

# MAP plot
data_plot <- as.matrix(SGibbs)
colnames(data_plot) <- rownames(data_plot) <- 1:NROW(dataset)
data_plot <- melt(data_plot)
data_plot <- rbind(
  data.frame(data_plot, Type = "Maximum a posteriori"),
  data.frame(data_plot, Type = "Minimum VI")
)
colnames(data_plot) <- c("Var1", "Var2", "Probability", "Type")

dummy_data_plot <- rbind(
  data.frame(x = c(38.5, 38 + 29 + 0.5), y = c(38.5, 38 + 29 + 0.5), Type = "Maximum a posteriori"),
  data.frame(x = c(44.5, 44 + 23 + 0.5), y = c(44.5, 44 + 23 + 0.5), Type = "Minimum VI")
)

p1_co <- ggplot(data_plot, aes(Var1, Var2, fill = Probability)) +
  geom_tile() +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  xlab("Slides") +
  ylab("Slides") +
  theme_light() +
  facet_grid(. ~ Type) +
  geom_vline(data = dummy_data_plot, aes(xintercept = x), linetype = "dashed") +
  geom_hline(data = dummy_data_plot, aes(yintercept = y), linetype = "dashed")

p1_co

# New patient probabilities ------------------------

x_new1 <- c(0, 1, 0, 0, 0, 0, 0)
x_new2 <- c(0, 1, 0, 0, 1, 0, 0)
x_new3 <- c(1, 1, 1, 0, 1, 0, 1)

probs1 <- numeric(k)
probs2 <- numeric(k)
probs3 <- numeric(k)

for (j in 1:k) {
  data_new <- rbind(x_new1, dataset[G_VI == j, ])
  probs1[j] <- exp(-(GBClust:::loss_binary(data_new) - GBClust:::loss_binary(dataset[G_VI == j, ])))
}
probs1 <- probs1 / sum(probs1)

for (j in 1:k) {
  data_new <- rbind(x_new2, dataset[G_VI == j, ])
  probs2[j] <- exp(-(GBClust:::loss_binary(data_new) - GBClust:::loss_binary(dataset[G_VI == j, ])))
}
probs2 <- probs2 / sum(probs2)

for (j in 1:k) {
  data_new <- rbind(x_new3, dataset[G_VI == j, ])
  probs3[j] <- exp(-(GBClust:::loss_binary(data_new) - GBClust:::loss_binary(dataset[G_VI == j, ])))
}
probs3 <- probs3 / sum(probs3)

tab <- rbind(probs1, probs2, probs3)
colnames(tab) <- c("Cluster I", "Cluster II", "Cluster III")
rownames(tab) <- c("Slide I", "Slide II", "Slide III")
kable(m_centroids, digits = 2)
kable(tab, digits = 2)
xtable(tab, digits = 2)

# Save graphs -----------------------------------------

ggsave("App_select.pdf", p_select, width = 8, height = 4)
ggsave("App_Co.pdf", p1_co, width = 9, height = 4)
