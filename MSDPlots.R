#install.packages("openxlsx")
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("arsenal")
#install.packages("ggpubr")
#install.packages("gridExtra")
#install.packages("broom")
library(openxlsx)
library(ggplot2)
library(dplyr)
#library(reshape2)
#library(arsenal)
#library(gridExtra)
#library(reshape2)
library(tidyr)
library(cowplot)
#library(stringr)
#library(broom)
library(ggpubr)


#-------------------Functions----------------
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

tukey_hsd <- function(x, ...){
  comparison <- comparison2 <- adj.p.value <- p.adj <-
    term <- group1 <- group2 <- NULL
  res <- TukeyHSD(x) %>%
    broom::tidy() %>%
    mutate(comparison2 = comparison) %>%
    separate(comparison2, into= c("group2", "group1"), sep = "-") %>%
    #rename(p.adj = adj.p.value) %>%
    mutate(adj.p.value = signif(adj.p.value, 2)) %>%
    select(term, group1, group2, everything())
  res
}
roundUp <- function(x,to=10){to*(x%/%to + as.logical(x%%to))}

Mean <- function(x) mean(x, na.rm = T)
SEM <- function(x) sd(x, na.rm = T)/sqrt(length(x[which(!is.na(x))]))

#--------Plot BoxPlots by Tissue and Treatment---------


setwd("/Users/swethagarimalla/Box/Bioinformatics (Swetha Garimalla)/Experiments/PRV-2019-004/MSD")
file <- "190812 PRV-2019-004 MSD V-PLEX 2-fold brain summary.xlsx"
manifest <- "PRE025 Tissue Manifest.xlsx"
study <- "PRV-2018-004"

wkbk <- loadWorkbook(file)
manifest <- loadWorkbook(manifest)

manifest <- as.data.frame(read.xlsx(manifest,
                                    sheet = "PRE025 Inventory",
                                    skipEmptyRows = TRUE,
                                    colNames = T));
manifest<-manifest[which(manifest$Tissue=="Rest Ctx"),]
manifest<-manifest[which(manifest$Segment=="Right A"),]

sheets <- getSheetNames(file) 
blues <- colorRampPalette(c("skyblue2", "skyblue4"))
input <- data.frame()
#sheets <- sheets[grepl("Replicates", sheets)]
for (sheet in sheets){
  
  data <- as.data.frame(read.xlsx(wkbk,
                                  sheet = sheet,
                                  skipEmptyRows = TRUE,
                                  colNames = T));

  data<-merge(data, manifest, by.x = "Sample.Name", by.y = "Sample.ID", all.x = TRUE)

  
  data$Tissue <- sheet
  data$Group<-as.character(data$Group) 
  data$`GCase/mg` <- as.numeric(data$`Normalized.conc.(pg/mg)`)
  
  gd_data <- data %>%
    group_by(Group) %>%
    summarise(mean = Mean(`GCase/mg`), sem = SEM(`GCase/mg`))
  
  data <- merge(data, gd_data, by = "Group", all.x = TRUE)
  data$Group <- as.character(data$Test.Article.ICV)
  data$Group[which(data$Sample.Name=="WT")] <- "WT"
  data$Sample <- as.character(data$Sample.ID)
  data$GCase <- as.numeric(data$`GCase/mg`)
  data<-unique(data)
  
  Combos <- aov(GCase ~ Group, data = data) %>%
    tukey_hsd(... = "Group")
  
  Combos <- Combos[which(Combos$adj.p.value < 0.2),]
  Combos <- rbind(Combos[which(Combos$group1 == "Excipient"),],
                  Combos[which(Combos$group2 == "Excipient"),])
  
  
  from <-max((as.numeric(data$GCase)), na.rm = T)
  num <- dim(Combos)[1] 
  y.position <-seq(from = from + from/10, by = from/10, length.out = num)
  
  data$Group<-factor(data$Group, levels=c("WT",
                                          "Excipient",
                                          "PR006A Virus Dose 1",
                                          "PR006A Virus Dose 2",
                                          "PR006A Virus Dose 3"))
  
  ggplot(data, aes(x = Group, y = GCase)) + 
    geom_point(size = 10, stroke=2, aes(color = Group))  + 
    #geom_point(size = 10, aes(shape = Sample), stroke=2)  + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text=element_text(size=25, face = "bold", hjust = 0.5),
          axis.title=element_text(size=30, face = "bold", hjust = 0.5),
          plot.title=element_text(size=32, face = "bold", hjust = 0.5),
          legend.position="none") +
    stat_pvalue_manual(data = Combos, label = "adj.p.value",
                       y.position = y.position,
                       tip.length = 0.005) +
    geom_point(mapping=aes(x=Group, y=mean, color = Group), 
               shape=95, size = 35) +
    scale_color_manual(values = c("grey", "red", blues(3))) +
    #ylim(0, 2600000) +
    labs(title = paste0("PRV-2018-004: ", sheet), 
         x="Group", 
         y = "Normalized conc (pg/mg)") 
  #  scale_x_discrete(labels = c( bquote("WT"), 
  #                               bquote("Excipient"),
  ##                               bquote("PR006A Virus Dose 1"),
  ##                               bquote("PR006A Virus Dose 2"),
  #                               bquote("PR006A Virus Dose 3")))
  ggsave(paste0(sheet, "_TwoFold.jpeg"))
  input<-rbind(input, data)
}


wkbk <- createWorkbook()
addWorksheet(wb = wkbk, "Table 1 -  Final Data Table", gridLines = TRUE) #add a worksheet to workbook. 
writeDataTable(wkbk, "Table 1 -  Final Data Table", input)

addWorksheet(wb = wkbk, "Table 2 -  Summary Data Table", gridLines = TRUE) #add a worksheet to workbook. 
writeDataTable(wkbk, "Table 2 -  Summary Data Table", unique(input[,c("Group", "Tissue", "mean", "sem")]))

saveWorkbook(wkbk, paste0(getwd(),"/RscriptEdited_",Sys.time(),".xlsx"), 
             overwrite = TRUE)






