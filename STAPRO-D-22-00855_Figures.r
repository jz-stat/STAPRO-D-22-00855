### packages are needed for the plots
library(EntropicStatistics)


### function of gse. input are sample frequency counts and order(s) of the gse
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


### Plots for Examples 1 and 2

## the three distributions have the same H. p1 and p2 GSE are similar when order is more than 1

p1 <- (1:3)^2/(sum((1:3)^2))
p2 <- c(p1, 0.00001)/sum(c(p1, 0.00001))
orders <- seq(0.01, 2, by = 0.01)
generalized.entropy(p1, orders)
generalized.entropy(p2, orders) 


c <- 8.2
K <- 8
p3 <- (1:K)^c/(sum((1:K)^c))


## Figure (a)
K <- 8
p4 <- dbinom(0:10, 10, 0.15)
p5 <- dpois(1:K, 2)/sum(dpois(1:K, 2))
orders <- seq(0.2, 2.5, by = 0.001)
diff <- round(generalized.entropy(p4, orders), 8) - round(generalized.entropy(p5, orders), 8)
op <- par(mar=c(5, 6, 4, 2) + 0.1, cex.lab = 1.5)
plot(diff~orders, type = "l", ylab = expression("H"^"(m)"~"(p)"~"- H"^"(m)"~"(q)"), xlab = "m (Orders of GSE)")
abline(h=0, lty = 2)
par(op)

## Figure (b)
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
### Heatmaps for simulation data constructed with R package EntropicStatistics

## Random Number Generator for Discrete Pareto

rpareto <- function(n, a){
  return(floor(runif(n)^{-1/a}))
}

## Obtain Binomial, Poisson, Pareto Random Data

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

sample_poisson_1 <- rpois(sample_size, 1)
sample_poisson_2 <- rpois(sample_size, 2)
sample_poisson_3 <- rpois(sample_size, 3)
sample_poisson_4 <- rpois(sample_size, 4)
sample_poisson_5 <- rpois(sample_size, 5)
sample_poisson_6 <- rpois(sample_size, 6)
sample_poisson_7 <- rpois(sample_size, 7)
sample_poisson_8 <- rpois(sample_size, 8)
sample_poisson_9 <- rpois(sample_size, 9)

sample_pareto_0.5 <- table(rpareto(sample_size, 0.5))
sample_pareto_0.6 <- table(rpareto(sample_size, 0.6))
sample_pareto_0.7 <- table(rpareto(sample_size, 0.7))
sample_pareto_0.8 <- table(rpareto(sample_size, 0.8))
sample_pareto_0.9 <- table(rpareto(sample_size, 0.9))
sample_pareto_1 <- table(rpareto(sample_size, 1))
sample_pareto_1.1 <- table(rpareto(sample_size, 1.1))
sample_pareto_1.2 <- table(rpareto(sample_size, 1.2))
sample_pareto_1.3 <- table(rpareto(sample_size, 1.3))
sample_pareto_1.4 <- table(rpareto(sample_size, 1.4))
sample_pareto_1.5 <- table(rpareto(sample_size, 1.5))


df <-
  list(
    binom_0.1 = sample_1,
    binom_0.2 = sample_2,
    binom_0.3 = sample_3,
    binom_0.4 = sample_4,
    binom_0.5 = sample_5,
    binom_0.6 = sample_6,
    binom_0.7 = sample_7,
    binom_0.8 = sample_8,
    binom_0.9 = sample_9,
    Poisson_1 = sample_poisson_1,
    Poisson_2 = sample_poisson_2,
    Poisson_3 = sample_poisson_3,
    Poisson_4 = sample_poisson_4,
    Poisson_5 = sample_poisson_5,
    Poisson_6 = sample_poisson_6,
    Poisson_7 = sample_poisson_7,
    Poisson_8 = sample_poisson_8,
    Poisson_9 = sample_poisson_9,
    Pareto_0.5 = sample_pareto_0.5,
    Pareto_0.6 = sample_pareto_0.6,
    Pareto_0.7 = sample_pareto_0.7,
    Pareto_0.8 = sample_pareto_0.8,
    Pareto_0.9 = sample_pareto_0.9,
    Pareto_1 = sample_pareto_1,
    Pareto_1.1 = sample_pareto_1.1,
    Pareto_1.2 = sample_pareto_1.2,
    Pareto_1.3 = sample_pareto_1.3,
    Pareto_1.4 = sample_pareto_1.4,
    Pareto_1.5 = sample_pareto_1.5
  )

## Make the plots

EntropicStatistics::HeatMap(df, title = "All Samples", title_text_size = 40, label_text_size = 40)
EntropicStatistics::HeatMap(df[c(1:9, 19:29)], title = "Binomial and Pareto Samples", title_text_size = 50, label_text_size = 50)
EntropicStatistics::HeatMap(df[10:18], title = "Poisson Samples", title_text_size = 50, label_text_size = 50)
EntropicStatistics::HeatMap(df[1:9], title = "Binomial Samples", title_text_size = 50, label_text_size = 50)
EntropicStatistics::HeatMap(df[19:29], title = "Pareto Samples", title_text_size = 50, label_text_size = 50)
EntropicStatistics::HeatMap(df[10:13], title = "Partial Poisson Samples", title_text_size = 50, label_text_size = 50)
EntropicStatistics::HeatMap(df[13:18], title = "Partial Poisson Samples", title_text_size = 50, label_text_size = 50)

### Heatmaps for Breast Cancer Wisconsin (Diagnostic) Data Set constructed with R package EntropicStatistics

data <- read.csv("~/Downloads/data.csv") #  Download from https://www.kaggle.com/datasets/uciml/breast-cancer-wisconsin-data
data_lists <- as.list(as.data.frame(apply(data[3:32], 2, function(x){table(cut(x, 10))})))
EntropicStatistics::HeatMap(data_lists, title = "Breast Cancer Wisconsin (Diagnostic) Data Set", title_text_size = 40, label_text_size = 40)
EntropicStatistics::HeatMap(data_lists, selection = c(4, 7, 10, 1, 6, 5, 3, 8, 9, 2), title = "Breast Cancer Wisconsin (Diagnostic) Data Set", title_text_size = 40, label_text_size = 40)
EntropicStatistics::HeatMap(data_lists, selection = c(14, 17, 20, 13, 11, 15, 19, 18, 12, 16), title = "Breast Cancer Wisconsin (Diagnostic) Data Set", title_text_size = 40, label_text_size = 40)
EntropicStatistics::HeatMap(data_lists, selection = c(30, 24, 29, 26, 23, 21, 27, 25, 22, 28), title = "Breast Cancer Wisconsin (Diagnostic) Data Set", title_text_size = 40, label_text_size = 40)
