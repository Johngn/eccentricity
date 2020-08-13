library(tidyverse)

ida <- read_csv("./data/dataida_mars_4.csv")
ida_small = data[1:400,]
cn <- read_csv("data/datatestcn.csv")

names(ida) <- c('planet_id','time','m','a','e','i')
names(cn) <- c("planet","time","a","eccentricity","inclination")

ggplot() +
  geom_point(data = ida_small, 
             aes(x = planet_id, y = e),
             color = 'orange')

ida %>%
  ggplot(aes(n, a)) +
  geom_col()
  
ggplot(data = ida) +
  geom_bar(mapping = aes(x = a))