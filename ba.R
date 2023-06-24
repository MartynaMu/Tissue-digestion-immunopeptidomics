# Author: Martyna Muszczek
# 20/06/2023
# Immunopeptidomics comparison of clinical tissue cut into 
# chunks vs digested to single cell suspension
# input: ba from netmhcpan from tumor chunks immunopeptidome 8-14 AA long

#Libraries----------------------------------------------------------------------
libs <- c("tidyverse", "ggplot2", "ggpubr", "DEP", "readxl")
lapply(libs, require, character.only = TRUE)

#Plot theme---------------------------------------------------------------------
theme_set(theme_DEP1())

#Load tumor chunks ids----------------------------------------------------------
#BA predictions from NetMHCpan for tumor chunks 8-14-mers
ba_TC <- read_xlsx("data/BA_tumor_chunks_8_14.xlsx")

#Category table-----------------------------------------------------------------
#Creating a function that will:
# - create a new df - "df2"
# - compute a new column based on the BA_rank value in original df - "df"
# - loop the action for all columns in a range specified in "cols" argument
ba_category <- function(df, cols, df2) {
  for (i2 in cols) {
    df2[,ncol(df2)+1] <- case_when(
      0.5 >= df[,i2] & df[,i2] <= 2 ~ "Strong binder",
      2 >= df[,i2] & df[,i2] >= 0.5 ~ "Weak binder",
      .default = "Non-binder"
    )
  }
  return(df2)
}

ba_TC_cont <- ba_TC[,1]
ba_TC_cont <- ba_category(ba_TC, 2:13, ba_TC_cont)
colnames(ba_TC_cont)[2:13] <- colnames(ba_TC)[2:13]

ba_TC_cat_L <- pivot_longer(ba_TC_cont, 2:13,
                           values_to = "Binding category", 
                           names_to = "HLA supertype")

#Load tumor digested ids----------------------------------------------------------
ba_TD <- read_xlsx("data/BA_tumor_digested_8_14.xlsx")

#Category table-----------------------------------------------------------------
ba_TD_cont <- ba_TD[,1]
ba_TD_cont <- ba_category(ba_TD, 2:13, ba_TD_cont)
colnames(ba_TD_cont)[2:13] <- colnames(ba_TD)[2:13]

ba_TD_cat_L <- pivot_longer(ba_TD_cont, 2:13,
                           values_to = "Binding category", 
                           names_to = "HLA supertype")

#Visualization of BA------------------------------------------------------------
#Tumor chunks 8-14-mers
tc_ba <- ggplot(ba_TC_cat_L, aes(`HLA supertype`, fill = `Binding category`))+
  geom_bar(width = .5, position = "fill")+
  scale_fill_grey(start = 0.8, end = 0.2)+
  labs(title = "Binding affinity prediction of tumor chunks 8-14-mers",
       y = "Share")+
  ggpubr::rremove("x.text")+
  ggpubr::rremove("legend")

#Tumor digested 8-14-mers
td_ba <- ggplot(ba_TD_cat_L, aes(`HLA supertype`, fill = `Binding category`))+
  geom_bar(width = .5, position = "fill")+
  scale_fill_grey(start = 0.8, end = 0.2)+
  labs(title = "Binding affinity prediction of digested tumor 8-14-mers",
       y = "Share")+
  ggpubr::rremove("x.text")

ggpubr::ggarrange(tc_ba, td_ba, ncol = 2, widths = c(1, 1.55))
