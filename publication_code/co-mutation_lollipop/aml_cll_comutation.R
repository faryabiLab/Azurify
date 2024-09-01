library(dplyr)
library(ggplot2)
library(GenVisR)
library(reshape2)

aml_input = ('filepath')
amls <- read.table(aml_input, sep ="\t",  header = TRUE, stringsAsFactors = FALSE, fill = TRUE, quote ="", na.strings ="NA")


amls$PLOT_EFFECT <- ifelse(grepl("inframe_deletion", amls$EFFECT), "In frame Del",
                          ifelse(grepl("inframe_insertion", amls$EFFECT), "In frame Ins",
                                 ifelse(grepl("synonymous_variant", amls$EFFECT), "Silent",
                                        ifelse(grepl("3'", amls$EFFECT), "3' UTR",
                                               ifelse(grepl("5'", amls$EFFECT), "5' UTR",
                                                      ifelse(grepl("splice", amls$EFFECT), "Splice Site",
                                                             ifelse(grepl("frameshift_variant", amls$EFFECT), "Frameshift",
                                                                    ifelse(grepl("missense", amls$EFFECT), "Missense", 
                                                                           ifelse(grepl("\\*", amls$P_CHANGE), "Nonsense", "Other")))))))))


most_deletrious <- c("Nonsense", "Missense", "Frameshift", "In frame Ins", "In frame Del", "Splice Site", "3'UTR", "5' UTR", "Silent", "Other")


subsamp <- subset(amls, select= c("SAMPLE", "GENE", "PLOT_EFFECT", "PCHANGE", "DIAGNOSIS", "FAF", "CATEGORIZATION"))
subsamp <- subsamp[!duplicated(subsamp), ]
names(subsamp) <-c("sample", "gene", "variant_class", "amino_acid_change_WU", "DX", "VAF","CATEGORIZATION")

waterfall(subsamp, plotMutBurden = FALSE, fileType = "Custom", variant_class_order= most_deletrious, mainGrid = FALSE, mainDropMut = TRUE, section_heights = c(0,10), maxGenes = 10)



cll_input = ('filepath')
cll <- read.table(cll_input, sep ="\t",  header = TRUE, stringsAsFactors = FALSE, fill = TRUE, quote ="", na.strings ="NA")

cll$PLOT_EFFECT <- ifelse(grepl("inframe_deletion", cll$EFFECT), "In frame Del",
                           ifelse(grepl("inframe_insertion", cll$EFFECT), "In frame Ins",
                                  ifelse(grepl("synonymous_variant", cll$EFFECT), "Silent",
                                         ifelse(grepl("3'", cll$EFFECT), "3' UTR",
                                                ifelse(grepl("5'", cll$EFFECT), "5' UTR",
                                                       ifelse(grepl("splice", cll$EFFECT), "Splice Site",
                                                              ifelse(grepl("frameshift_variant", cll$EFFECT), "Frameshift",
                                                                     ifelse(grepl("missense", cll$EFFECT), "Missense", 
                                                                            ifelse(grepl("\\*", cll$P_CHANGE), "Nonsense", "Other")))))))))


most_deletrious <- c("Nonsense", "Missense", "Frameshift", "In frame Ins", "In frame Del", "Splice Site", "3'UTR", "5' UTR", "Silent", "Other")


subsamp <- subset(cll, select= c("SAMPLE", "GENE", "PLOT_EFFECT", "PCHANGE", "DIAGNOSIS", "FAF", "CATEGORIZATION"))
subsamp <- subsamp[!duplicated(subsamp), ]
names(subsamp) <-c("sample", "gene", "variant_class", "amino_acid_change_WU", "DX", "VAF","CATEGORIZATION")
main_layer <- theme_grey()
waterfall(subsamp, plotMutBurden = FALSE, fileType = "Custom", variant_class_order= most_deletrious, mainGrid = FALSE, mainDropMut = TRUE, section_heights = c(0,10), maxGenes = 10)
