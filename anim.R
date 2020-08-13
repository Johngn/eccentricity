library(tidyverse)
library(gganimate)

data = read_csv('./data/dataida_mars_4.csv')

colnames(data) <- c('planet_id','time','m','a','e','i')

data_small = data[1:400,]

ggplot(data_small, aes(a, e, size = m)) +
  geom_point() +
  transition_time(time)

gganimate(g)

p <- ggplot(data_small, aes(x = a, y = e)) +
  geom_point()

plot(p)
anim <- p +
  transition_states(Species,
                    transition_length = 2,
                    state_length = 1)

anim

library(gapminder)

ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, colour = country)) +
  geom_point(alpha = 0.7, show.legend = FALSE) +
  scale_colour_manual(values = country_colors) +
  scale_size(range = c(2, 12)) +
  scale_x_log10() +
  facet_wrap(~continent) +
  # Here comes the gganimate specific bits
  labs(title = 'Year: {frame_time}', x = 'GDP per capita', y = 'life expectancy') +
  transition_time(year) +
  ease_aes('linear')
