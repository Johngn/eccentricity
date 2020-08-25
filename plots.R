library(tidyverse)

ida <- read_csv("./data/ida_collision.csv")
cn <- read_csv("./data/cn_collision.csv")

data <- read.delim('./data/Outida_7_000100000000.dat') 

names(ida) <- c('a','e','i','m')
names(cn) <- c('a','e','i','m')

ggplot(ida) +
  geom_histogram(aes(i)) +
  geom_density(aes(x = i))

ida %>%
  ggplot(aes(n, a)) +
  geom_col()
  
ggplot(data = ida) +
  geom_bar(mapping = aes(x = a))