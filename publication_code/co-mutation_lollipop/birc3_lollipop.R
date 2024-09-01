library(dplyr)
library(ggplot2)
library(GenVisR)
library(reshape2)

birc3 = ('file')
azb3 <- read.table(birc3, sep ="\t",  header = TRUE, stringsAsFactors = FALSE, fill = TRUE, quote ="", na.strings ="NA")



azb3$AA_NUM <- as.numeric(gsub("\\D", "", azb3$PCHANGE))
azb3$AA_NUM <- sub("^(\\d{3}).*$", "\\1", azb3$AA_NUM)
azb3$AA_NUM <- as.character(azb3$AA_NUM)
azb3$AA_CNT <- as.numeric(ave(azb3$AA_NUM, azb3$AA_NUM, FUN = length))
azb3$AA_NUM <- as.numeric(azb3$AA_NUM)
azb3$AA_CNT <- as.numeric(azb3$AA_CNT)

train_birc3
cpdtrain <- read.table(train_birc3, sep ="\t",  header = TRUE, stringsAsFactors = FALSE, fill = TRUE, quote ="", na.strings ="NA")
cpdtrain$AA_NUM <- as.numeric(gsub("\\D", "", cpdtrain$PCHANGE))
cpdtrain$AA_NUM <- as.character(cpdtrain$AA_NUM)
cpdtrain$AA_NUM <- sub("^(\\d{3}).*$", "\\1", cpdtrain$AA_NUM)
cpdtrain$AA_CNT <- as.numeric(ave(cpdtrain$AA_NUM, cpdtrain$AA_NUM, FUN = length))



#create a merged dataframe of cpdtrain and azurify results
azb3$method <- "Azurify"
cpdtrain$method <- "Learning Set"
azb3$PLOT_EFFECT <- ifelse(grepl("frameshift", azb3$EFFECT), "Frameshift", 
                         ifelse(grepl("missense", azb3$EFFECT), "Missense", 
                                ifelse(grepl("stop", azb3$EFFECT), "Truncating", 
                                       ifelse(grepl("inframe", azb3$EFFECT), "In frame Del", "Other"))))


cpdtrain$PLOT_EFFECT <- ifelse(grepl("frameshift", cpdtrain$EFFECT), "Frameshift", 
                           ifelse(grepl("missense", cpdtrain$EFFECT), "Missense", 
                                  ifelse(grepl("stop", cpdtrain$EFFECT), "Truncating", 
                                         ifelse(grepl("inframe", cpdtrain$EFFECT), "In frame Del", "Other"))))




az_sub<- subset(azb3, select = c("AA_NUM", "AA_CNT", "method", "PLOT_EFFECT", "comp_pred"))
names(az_sub)[names(az_sub) == 'comp_pred'] <- "CATEGORIZATION"
clin_sub<- subset(cpdtrain, select = c("AA_NUM", "AA_CNT", "method", "PLOT_EFFECT", "CATEGORIZATION"))
combined_df <- rbind(az_sub, clin_sub)
combined_df <- na.omit(combined_df)
combined_df$AA_NUM <- as.numeric(combined_df$AA_NUM)
combined_df$CATEGORIZATION <-sub("Disease Associated", "Pathogenic", combined_df$CATEGORIZATION)
combined_df$CATEGORIZATION <-sub("VOUS", "VUS", combined_df$CATEGORIZATION)
combined_df$CATEGORIZATION <- factor(combined_df$CATEGORIZATION,               
                         levels = c("Pathogenic", "VUS", "Likely Benign"))


#select focused locations
combined_df <- combined_df %>% filter(AA_NUM > 520 & AA_NUM < 575)


p_changes <-ggplot() +
  geom_segment( data = combined_df, aes(x=AA_NUM, xend=AA_NUM, y=0, yend=AA_CNT,linetype=method) ) +
  geom_point(data =combined_df, inherit.aes = FALSE, aes(x=AA_NUM, y=AA_CNT, color=CATEGORIZATION), size=3) + 
  xlab("Amino Acid") +  ylab("# BIRC3 Mutations")  + guides(fill=guide_legend(title="New Legend Title")) + 
  scale_color_manual(values=c("red","#011F5b", "#A5ACAF")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position = "top", legend.title=element_blank())

p_changes
#comment or uncomment for zoom/no zoom
#pc1 <- p_changes + annotate("rect", xmin=527, xmax = 572, ymin = 0, ymax = -0.8)
#pc1

pc1 <- p_changes + annotate("rect", xmin=0, xmax = 650, ymin = 0, ymax = -0.8)
pc2 <- pc1 + annotate("rect", xmin=32, xmax = 97, ymin = 0.1, ymax = -0.9, col="#27AE60", fill="#27AE60")
pc3 <- pc2 + annotate("rect", xmin=172, xmax = 235, ymin = 0.1, ymax = -0.9,col="#27AE60", fill="#27AE60")
pc4 <- pc3 + annotate("rect", xmin=258, xmax = 322, ymin = 0.1, ymax = -0.9, col="#27AE60", fill="#27AE60")
pc5 <- pc4 + annotate("rect", xmin=446, xmax = 526, ymin = 0.1, ymax = -0.9, col="#FF5733", fill="#FF5733")
pc6 <- pc5 + annotate("rect", xmin=554, xmax = 596, ymin = 0.1, ymax = -0.9, col="#2874A6", fill="#2874A6")
pc7 <- pc6 + annotate("text", x = 65, y = -0.4, label = c("BIR"), col="white")
pc7 <- pc7 + annotate("text", x = 203, y = -0.4, label = c("BIR"), col="white")
pc7 <- pc7 + annotate("text", x = 290, y = -0.4, label = c("BIR"), col="white")
pc7 <- pc7 + annotate("text", x = 490, y = -0.4, label = c("CARD"), col="white")
pc7 <- pc7 + annotate("text", x = 575, y = -0.4, label = c("ZF"), col="white")
pc7

#lets create a zoomed in view of the ZF domain
