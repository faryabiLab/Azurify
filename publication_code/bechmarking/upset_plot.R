library(UpSetR)


#the input is a csv file marking if a mutation is present or absent in a resource
#example 
#Azurify, 1, ClinVar, 1, UPennTest, 0
e2 = read.csv("emergent_pathogenic.csv", sep=',')

pdf('emergent_path.pdf')
upset(e2, sets = c("Azurify", "Reported", "ClinVar"), 
      mainbar.y.label = "Azurify Discordant Pathogenic Variants",order.by = "freq", main.bar.color = 'dark red', mb.ratio = c(0.55, 0.45),text.scale = c(1.5, 1.6, 1.5, 1.5), point.size = 3.5, line.size = 2)
dev.off()