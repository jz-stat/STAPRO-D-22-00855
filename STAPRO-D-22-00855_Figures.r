## packages are needed for the plots
library(EntropicStatistics)


# function of gse. input are sample frequency counts and order(s) of the gse
generalized.entropy <- function(sample_freq, order){
  results <- vector()
  n <- sum(sample_freq)
  phats <- sample_freq/n
  phats <- phats[phats>0]
  for (i in 1:length(order)) {
    pmhats <- phats^order[i]/sum(phats^order[i])
    results[i] <- sum(-pmhats*log(pmhats))
  }
  return(results)
}


### Examples 1 and 2

# the three distributions have the same H. p1 and p2 GSE are similar when order is more than 1

p1 <- (1:3)^2/(sum((1:3)^2))
p2 <- c(p1, 0.00001)/sum(c(p1, 0.00001))
orders <- seq(0.01, 2, by = 0.01)
generalized.entropy(p1, orders)
generalized.entropy(p2, orders) 


c <- 8.2
K <- 8
p3 <- (1:K)^c/(sum((1:K)^c))


## Figure 2(a)
K <- 8
p4 <- dbinom(0:10, 10, 0.15)
p5 <- dpois(1:K, 2)/sum(dpois(1:K, 2))
orders <- seq(0.2, 2.5, by = 0.001)
diff <- round(generalized.entropy(p4, orders), 8) - round(generalized.entropy(p5, orders), 8)
op <- par(mar=c(5, 6, 4, 2) + 0.1, cex.lab = 1.5)
plot(diff~orders, type = "l", ylab = expression("H"^"(m)"~"(p)"~"- H"^"(m)"~"(q)"), xlab = "m (Orders of GSE)")
abline(h=0, lty = 2)
par(op)

## Figure 2(b)
K <- 8
bp <- 0.25
c <- 2.05
p4 <- dbinom(0:(K-1), K-1, bp)
p5 <- dpois(1:K, c)/sum(dpois(1:K, c))
orders <- seq(0.01, 15, by = 0.01)
diff <- round(generalized.entropy(p4, orders), 8) - round(generalized.entropy(p5, orders), 8)
op <- par(mar=c(5, 6, 4, 2) + 0.1, cex.lab = 1.5)
plot(diff~orders, type = "l", ylab = expression("H"^"(m)"~"(p)"~"- H"^"(m)"~"(q)"), xlab = "m (Orders of GSE)")
abline(h=0, lty = 2)
par(op)


#################
## heatmap for Figure 1

library(latticeExtra) 

generalized.entropy <- function(sample_freq, order){
  results <- vector()
  n <- sum(sample_freq)
  phats <- sample_freq/n
  phats <- phats[phats>0]
  for (i in 1:length(order)) {
    pmhats <- phats^order[i]/sum(phats^order[i])
    results[i] <- sum(-pmhats*log(pmhats))
  }
  return(results)
}

orders <- seq(0.50, 3, by = 0.01)

binom_n <- 10
sample_size <- 10000
sample_1 <- table(rbinom(size=binom_n, n=sample_size, 0.1))
sample_2 <- table(rbinom(size=binom_n, n=sample_size, 0.2))
sample_3 <- table(rbinom(size=binom_n, n=sample_size, 0.3))
sample_4 <- table(rbinom(size=binom_n, n=sample_size, 0.4))
sample_5 <- table(rbinom(size=binom_n, n=sample_size, 0.5))
sample_6 <- table(rbinom(size=binom_n, n=sample_size, 0.6))
sample_7 <- table(rbinom(size=binom_n, n=sample_size, 0.7))
sample_8 <- table(rbinom(size=binom_n, n=sample_size, 0.8))
sample_9 <- table(rbinom(size=binom_n, n=sample_size, 0.9))

sample_poisson_1 <- rpois(10000, 1)
sample_poisson_2 <- rpois(10000, 2)
sample_poisson_3 <- rpois(10000, 3)
sample_poisson_4 <- rpois(10000, 4)
sample_poisson_5 <- rpois(10000, 5)
sample_poisson_6 <- rpois(10000, 6)
sample_poisson_7 <- rpois(10000, 7)
sample_poisson_8 <- rpois(10000, 8)
sample_poisson_9 <- rpois(10000, 9)


gse_1 <- generalized.entropy(sample_1, orders)
gse_2 <- generalized.entropy(sample_2, orders)
gse_3 <- generalized.entropy(sample_3, orders)
gse_4 <- generalized.entropy(sample_4, orders)
gse_5 <- generalized.entropy(sample_5, orders)
gse_6 <- generalized.entropy(sample_6, orders)
gse_7 <- generalized.entropy(sample_7, orders)
gse_8 <- generalized.entropy(sample_8, orders)
gse_9 <- generalized.entropy(sample_9, orders)
gse_10 <- generalized.entropy(sample_poisson_1, orders)
gse_11 <- generalized.entropy(sample_poisson_2, orders)
gse_12 <- generalized.entropy(sample_poisson_3, orders)
gse_13 <- generalized.entropy(sample_poisson_4, orders)
gse_14 <- generalized.entropy(sample_poisson_5, orders)
gse_15 <- generalized.entropy(sample_poisson_6, orders)
gse_16 <- generalized.entropy(sample_poisson_7, orders)
gse_17 <- generalized.entropy(sample_poisson_8, orders)
gse_18 <- generalized.entropy(sample_poisson_9, orders)



gse_samples <- function(sample, orders){
  generalized.entropy <- function(sample_freq, order){
    results <- vector()
    n <- sum(sample_freq)
    phats <- sample_freq/n
    phats <- phats[phats>0]
    for (i in 1:length(order)) {
      pmhats <- phats^order[i]/sum(phats^order[i])
      results[i] <- sum(-pmhats*log(pmhats))
    }
    return(results)
  }
  results <- vector("list", length(sample))
  results <- lapply(sample, function(x) generalized.entropy(x, orders))
}


df1 <- t(data.frame(binom_0.1 = gse_1, binom_0.2 = gse_2, binom_0.3 = gse_3, binom_0.4 = gse_4, binom_0.5 = gse_5, binom_0.6 = gse_6, binom_0.7 = gse_7, binom_0.8 = gse_8, binom_0.9 = gse_9, Poisson_1 = gse_10, Poisson_2 = gse_11, Poisson_3 = gse_12, Poisson_4 = gse_13, Poisson_5 = gse_14, Poisson_6 = gse_15, Poisson_7 = gse_16, Poisson_8 = gse_17, Poisson_9 = gse_18))
colnames(df1) <- orders

library(ggplot2)
library(tidyverse)

# Figure 1 (a)

as.data.frame(df1) |>
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  ggplot(aes(factor(name, unique(name)), rowname, fill = value)) + 
  ggtitle(expression("GSEs with orders between 0.5 and 3.0: Binominal and Poisson")) +
  scale_x_discrete(labels = ~., breaks = ~ c(seq(6, 28, 2)/10)) + theme(plot.title = element_text(hjust = 0.5, size = 25), text = element_text(size=25)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("blue4", "white", "red3")) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = 2))

# Figure 1 (b)

as.data.frame(df1[1:9,]) |>
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  ggplot(aes(factor(name, unique(name)), rowname, fill = value)) +
  ggtitle(expression("GSEs with orders between 0.5 and 3.0: Binomial distribution p = 0.1, 0.2, ... ,0.9")) +
  scale_x_discrete(labels = ~., breaks = ~ c(seq(6, 28, 2)/10)) + theme(plot.title = element_text(hjust = 0.5, size = 25), text = element_text(size=25)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("blue4", "white", "red3")) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = 2))

# Figure 1 (c)

as.data.frame(df1[10:18,]) |>
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  ggplot(aes(factor(name, unique(name)), rowname, fill = value)) +
  ggtitle(expression("GSEs with orders between 0.5 and 3.0: Poisson distributions "~lambda~"= 1, 2, .., 9")) +
  scale_x_discrete(labels = ~., breaks = ~ c(seq(6, 28, 2)/10)) + theme(plot.title = element_text(hjust = 0.5, size = 25), text = element_text(size=25)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("blue4", "white", "red3")) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = 2))

# Figure 1 (d)

as.data.frame(df1[10:13,]) |>
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  ggplot(aes(factor(name, unique(name)), rowname, fill = value)) +
  ggtitle(expression("GSEs with orders between 0.5 and 3.0: Poisson distributions "~lambda~"= 1, 2, 3, 4")) +
  scale_x_discrete(labels = ~., breaks = ~ c(seq(6, 28, 2)/10)) + theme(plot.title = element_text(hjust = 0.5, size = 25), text = element_text(size=25)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("blue4", "white", "red3")) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = 2))

# Figure 1 (e)

as.data.frame(df1[13:18,]) |>
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  ggplot(aes(factor(name, unique(name)), rowname, fill = value)) +
  ggtitle(expression("GSEs with orders between 0.5 and 3.0: Poisson distributions "~lambda~"= 4, 5, 6, 7, 8, 9")) + 
  scale_x_discrete(labels = ~., breaks = ~ c(seq(6, 28, 2)/10)) + theme(plot.title = element_text(hjust = 0.5, size = 25), text = element_text(size=25)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("blue4", "white", "red3")) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = 2))

as.data.frame(df1[11:14,]) |>
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  ggplot(aes(factor(name, unique(name)), rowname, fill = value)) +
  ggtitle(expression("GSEs with orders between 0.5 and 3.0: Poisson distributions "~lambda~"= 2, 3, 4, 5")) + 
  scale_x_discrete(labels = ~., breaks = ~ c(seq(6, 28, 2)/10)) + theme(plot.title = element_text(hjust = 0.5, size = 25), text = element_text(size=25)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("blue4", "white", "red3")) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = 2))

# Figure 1 (f)

as.data.frame(df1[12:15,]) |>
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  ggplot(aes(factor(name, unique(name)), rowname, fill = value)) +
  ggtitle(expression("GSEs with orders between 0.5 and 3.0: Poisson distributions "~lambda~"= 3, 4, 5, 6")) + 
  scale_x_discrete(labels = ~., breaks = ~ c(seq(6, 28, 2)/10)) + theme(plot.title = element_text(hjust = 0.5, size = 25), text = element_text(size=25)) + 
  geom_tile() +
  scale_fill_gradientn(colours = c("blue4", "white", "red3")) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none", 
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        text = element_text(face = 2))


