library(ggplot2)
library(tidyverse)
library(ggalluvial)

### ALLUVIAL

df = read.csv("catanal.csv", sep=',')

colorfill1 <- c("#002f6c", "#75B2DD", "#FFC700", "#990000","#A5ACAF")
colorfill2 <- c("#A5ACAF", "#990000", "#FFC700", "#75B2DD","#002f6c")
ggplot(data = df,aes(axis1 = PREVIOUS_CATEGORIZATION, axis2 = CATEGORIZATION, y = Freq)) + geom_alluvium(aes(fill = PREVIOUS_CATEGORIZATION),show.legend = FALSE) +
  geom_stratum(aes(fill = PREVIOUS_CATEGORIZATION), show.legend = FALSE) + geom_text(stat = "stratum", size = 6, aes(label = after_stat(stratum))) +
  scale_fill_manual(values = colorfill1, na.value = colorfill2) +theme_void() +theme( axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


ggsave('alluvial_cat.pdf', width = 10, height = 7, device='png', dpi=600)