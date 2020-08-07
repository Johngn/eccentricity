library(tidyverse)

ida <- read_csv("summerproject/data/datatest.csv")
cn <- read_csv("summerproject/data/datatestcn.csv")

names(ida) <- c("planet","time","a","eccentricity","inclination")
names(cn) <- c("planet","time","a","eccentricity","inclination")

ggplot() +
  geom_point(data = ida, 
             aes(x = a, y = eccentricity),
             color = 'orange') +
  geom_point(data = cn,
             aes(x = a, y = eccentricity),
             color = 'dodgerblue1') +
  geom_smooth(data = ida, 
              aes(x = a, y = eccentricity),
              color = 'orange') +
  geom_smooth(data = cn,
              aes(x = a, y = eccentricity),
              color = 'dodgerblue1')

ida %>%
  ggplot(aes(n, a)) +
  geom_col()
  
ggplot(data = ida) +
  geom_bar(mapping = aes(x = a))