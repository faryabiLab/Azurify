library(dplyr)
library(ggplot2)
library(GenVisR)
library(reshape2)

srsf2 = ('filepath')
azsr <- read.table(srsf2, sep ="\t",  header = TRUE, stringsAsFactors = FALSE, fill = TRUE, quote ="", na.strings ="NA")


azsr$AA_NUM <- sapply(strsplit(azsr$PCHANGE, "_"), `[`, 1)
azsr$AA_NUM <- as.numeric(gsub("\\D", "", azsr$AA_NUM))
azsr$AA_NUM <- sub("^(\\d{3}).*$", "\\1", azsr$AA_NUM)
azsr$AA_NUM <- as.character(azsr$AA_NUM)
azsr$AA_CNT <- as.numeric(ave(azsr$AA_NUM, azsr$AA_NUM, FUN = length))
azsr$AA_NUM <- as.numeric(azsr$AA_NUM)
azsr$AA_CNT <- as.numeric(azsr$AA_CNT)

train_srsf2 = ('filepath')
cpdtrain <- read.table(train_srsf2, sep ="\t",  header = TRUE, stringsAsFactors = FALSE, fill = TRUE, quote ="", na.strings ="NA")
cpdtrain$AA_NUM <- sapply(strsplit(cpdtrain$PCHANGE, "_"), `[`, 1)
cpdtrain$AA_NUM <- as.numeric(gsub("\\D", "", cpdtrain$AA_NUM))
cpdtrain$AA_NUM <- as.character(cpdtrain$AA_NUM)
cpdtrain$AA_NUM <- sub("^(\\d{3}).*$", "\\1", cpdtrain$AA_NUM)
cpdtrain$AA_CNT <- as.numeric(ave(cpdtrain$AA_NUM, cpdtrain$AA_NUM, FUN = length))



#create a merged dataframe of cpdtrain and azurify results
azsr$method <- "Azurify"
cpdtrain$method <- "Learning Set"
azsr$PLOT_EFFECT <- ifelse(grepl("frameshift", azsr$EFFECT), "Frameshift", 
                         ifelse(grepl("missense", azsr$EFFECT), "Missense", 
                                ifelse(grepl("stop", azsr$EFFECT), "Truncating", 
                                       ifelse(grepl("inframe", azsr$EFFECT), "In frame Del", "Other"))))


cpdtrain$PLOT_EFFECT <- ifelse(grepl("frameshift", cpdtrain$EFFECT), "Frameshift", 
                           ifelse(grepl("missense", cpdtrain$EFFECT), "Missense", 
                                  ifelse(grepl("stop", cpdtrain$EFFECT), "Truncating", 
                                         ifelse(grepl("inframe", cpdtrain$EFFECT), "In frame Del", "Other"))))




az_sub<- subset(azsr, select = c("AA_NUM", "AA_CNT", "method", "PLOT_EFFECT", "comp_pred"))
names(az_sub)[names(az_sub) == 'comp_pred'] <- "CATEGORIZATION"
clin_sub<- subset(cpdtrain, select = c("AA_NUM", "AA_CNT", "method", "PLOT_EFFECT", "CATEGORIZATION"))
combined_df <- rbind(az_sub, clin_sub)
combined_df <- na.omit(combined_df)
combined_df$AA_NUM <- as.numeric(combined_df$AA_NUM)
combined_df$CATEGORIZATION <-sub("Disease Associated", "Pathogenic", combined_df$CATEGORIZATION)
combined_df$CATEGORIZATION <-sub("VOUS", "VUS", combined_df$CATEGORIZATION)
combined_df$CATEGORIZATION <- factor(combined_df$CATEGORIZATION,               
                         levels = c("Pathogenic", "VUS", "Likely Benign"))



p_changes <-ggplot() +
  geom_segment( data = combined_df, aes(x=AA_NUM, xend=AA_NUM, y=0, yend=AA_CNT,linetype=method)) + coord_cartesian(ylim=c(0, 25)) + #ylim(0, 30) +
  geom_jitter(data =combined_df, width=0.2, height=0.2, inherit.aes = FALSE, aes(x=AA_NUM, y=AA_CNT, color=CATEGORIZATION), size=3) +
  xlab("Amino Acid") +  ylab("SRSF2 Mutations")  + guides(fill=guide_legend(title="New Legend Title")) + 
  scale_color_manual(values=c("red","#011F5b", "#A5ACAF")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position = "top", legend.title=element_blank())

pc1 <- p_changes + annotate("rect", xmin=0, xmax = 225, ymin = 0, ymax = -0.8)
pc2 <- pc1 + annotate("rect", xmin=14, xmax = 92, ymin = 0.1, ymax = -0.9, col="#27AE60", fill="#27AE60")
pc3 <- pc2 + annotate("text", x = 55, y = -0.4, label = c("RPM"), col="white")
pc3



#select focused locations
combined_df_zoom <- combined_df %>% filter(AA_NUM > 70 & AA_NUM < 120)

p_changes <-ggplot() +
  geom_segment( data = combined_df_zoom, aes(x=AA_NUM, xend=AA_NUM, y=0, yend=AA_CNT,linetype=method)) + coord_cartesian(ylim=c(0, 25)) + #ylim(0, 30) +
  geom_jitter(data =combined_df_zoom, width=0.2, height=0.2, inherit.aes = FALSE, aes(x=AA_NUM, y=AA_CNT, color=CATEGORIZATION), size=3) +
  xlab("Amino Acid") +  ylab("SRSF2 Mutations")  + guides(fill=guide_legend(title="New Legend Title")) + 
  scale_color_manual(values=c("red","#011F5b", "#A5ACAF")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank(), legend.position = "top", legend.title=element_blank())
p_changes
