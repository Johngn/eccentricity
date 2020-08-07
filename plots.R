library(tidyverse)

ida <- read_csv("summerproject/data/datatest.csv")
cn <- read_csv("summerproject/data/datatestcn.csv")

names(ida) <- c("planet","time","a","eccentricity","inclination")
names(cn) <- c("planet","time","a","eccentricity","inclination")

ida %>%
  ggplot(aes(x = a, y = eccentricity)) +
  geom_point() +
  geom_point(cn)
  labs(title = "")

ggplot() +
  geom_point(data = ida, 
             aes(x = a, y = eccentricity),
             color = 'orange') +
  geom_point(data = cn,
             aes(x = a, y = eccentricity),
             color = 'blue') +
  geom_smooth(data = ida, 
              aes(x = a, y = eccentricity),
              color = 'orange') +
  geom_smooth(data = cn,
              aes(x = a, y = eccentricity),
              color = 'blue')
