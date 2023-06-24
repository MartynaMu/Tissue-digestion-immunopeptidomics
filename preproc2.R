
# Author: Martyna Muszczek
# 19/06/2023
# Immunopeptidomics comparison of clinical tissue cut into 
# chunks vs digested to single cell suspension
# input: diann set to norm, spec lib from fragpipe

#Libraries----------------------------------------------------------------------
libs <- c("tidyverse", "ggplot2", "proBatch", "diann", "ggpubr", "PCAtools", "DEP", 
          "ComplexHeatmap", "eulerr", "scrime")
lapply(libs, require, character.only = TRUE)

#Plot theme---------------------------------------------------------------------
theme_set(theme_DEP1())

#Preprocessing------------------------------------------------------------------
##Matrix load===================================================================
dfRaw <- read.delim("data/report.tsv", dec = ".")

##Filter out iRTs and contams===================================================
dfRaw <- filter(dfRaw, !grepl("Biognosys", dfRaw$Protein.Ids))
dfRaw <- filter(dfRaw, !grepl("contam", dfRaw$Protein.Ids))

##maxLFQ========================================================================
peps_maxlfq <- diann_maxlfq(dfRaw[dfRaw$Q.Value <= 0.01 & dfRaw$PG.Q.Value <= 0.01,], 
                                group.header="Modified.Sequence", id.header = "Precursor.Id", 
                                quantity.header = "Precursor.Normalised")
##Tidying=======================================================================
#converting output matrix to df
maxLFQ_ip <- as.data.frame(peps_maxlfq)

#row names to column
maxLFQ_ip <- tibble::rownames_to_column(maxLFQ_ip, "Modified.Sequence")

#reordering columns 1>12
maxLFQ_ip <- maxLFQ_ip[,c(1, 10, 2, 3, 4, 6, 5, 7, 8, 9, 11, 12)]

#adding important columns
#take modified.sequence which will be a shared column and a column to be added, 
#left_join by shared column
#select columns to add, filter for unique modified.seq and join
add <- dfRaw %>% select(Protein.Ids, Genes, Stripped.Sequence, Modified.Sequence)

id <- add %>% group_by(Protein.Ids, Genes, Stripped.Sequence) %>% 
  summarise(Modified.Sequence = unique(Modified.Sequence)) %>% ungroup()

df_new <- id %>% left_join(y = maxLFQ_ip, by = "Modified.Sequence")

#Identifying outlier runs
boxplot(log2(df_new[5:15]))
hist(log2(df_new[5:15]))

#Subsetting runs with more similar TIC
df <- df_new[,-c(5, 8, 9)]
df <- as.data.frame(df)

write_delim(df, "data/df.txt", delim = "\t")

##Tech rep merge================================================================
#Function to impute missing observations from other technical rep, then average them
merge2reps <- function(df, col, x) {
  df[, col] <- if_else( #pasting value from next tech rep in case of NA
    is.na(df[, col]) == TRUE, 
    df[, col+1],
    df[, col])
  df[, col+1] <- if_else( #pasting a value from previous tech rep in case of NA
    is.na(df[, col+1]) == TRUE, 
    df[, col],
    df[, col+1])
  new <- rowMeans(df[, c(col, (col+1))]) #creating a new variable with mean result of every 2 cols
  df[,(ncol(df)+1)] <- new #adding the created variable to the df
  colnames(df)[ncol(df)] <- paste0("Sample", x) 
  return(df)
}

for (i in seq.int(5, 11, 2)) { #seq from 1st num col to the last one by every 2
  e = 1
  df <- merge2reps(df = df, col = i, x = e)
  e <- e + 1
}

df_m <- df[,-(5:12)] #storing averaged tech reps
#Renaming samples
colnames(df_m)[5:8] <- c("Tumor digested", 
                         "Lung chunks", 
                         "Tumor chunks", 
                         "Lung digested")
df_m <- df_m[,c(1, 2, 3, 4, 7, 5, 6, 8)]

#Formatting---------------------------------------------------------------------
dfWLog <- mutate(df_m, log2(df_m[5:8]), .keep = "all")
dfLLog <- pivot_longer(dfWLog, 5:8, values_to = "Intensity",
                       names_to = "Sample", values_drop_na = TRUE)

#Summaries----------------------------------------------------------------------
sum_stats <- dfLLog %>% group_by(Sample) %>% summarise(
            IDS = n(),
            Mean = mean(Intensity),
            Median = median(Intensity),
            SD = sd(Intensity),
            Var = var(Intensity)) %>% ungroup()

#Length distribution------------------------------------------------------------
dfLLog$AA.length <- str_length(dfLLog$Stripped.Sequence)

ggplot(dfLLog, aes(AA.length)) +
  geom_bar()+
  facet_wrap(~Sample)+
  scale_x_continuous(breaks = seq.int(6,30,1))+
  scale_y_continuous(breaks = seq.int(0,700,100))

#Retrieving peptide IDS for all samples
tumor_chunks_ids <- dfLLog$Stripped.Sequence[dfLLog$Sample == "Tumor chunks"]
write_lines(tumor_chunks_ids, "data/tumor_chunks_ids.txt")
tumor_digested_ids <- dfLLog$Stripped.Sequence[dfLLog$Sample == "Tumor digested"]
write_lines(tumor_digested_ids, "data/tumor_digested_ids.txt")
lung_chunks_ids <- dfLLog$Stripped.Sequence[dfLLog$Sample == "Lung chunks"]
write_lines(lung_chunks_ids, "data/lung_chunks_ids.txt")
lung_digested_ids <- dfLLog$Stripped.Sequence[dfLLog$Sample == "Lung digested"]
write_lines(lung_digested_ids, "data/lung_digested_ids.txt")

#Venn diagrams------------------------------------------------------------------
euler_lung <- euler(list(
  "Lung chunks" = unique(lung_chunks_ids), 
  "Lung digested" = unique(lung_digested_ids)
  ))
plot(euler_lung, quantities = TRUE)

euler_tumor <- euler(list(
  "Tumor chunks" = unique(tumor_chunks_ids),
  "Tumor digested" = unique(tumor_digested_ids)
))
plot(euler_tumor, quantities = TRUE)

#Intersections------------------------------------------------------------------
##lung====
lung_int <- enframe(
                    list(
                      "Intersection" = intersect(unique(lung_chunks_ids),
                                                 unique(lung_digested_ids)),
                      "Chunks only" = setdiff(unique(lung_chunks_ids),
                                              unique(lung_digested_ids)),
                      "Digested only" = setdiff(unique(lung_digested_ids),
                                                unique(lung_chunks_ids))
                      ),
                    "Set", "Ids"
)
lung_int <- unnest_longer(lung_int, Ids)
lung_int$AA.length <- str_length(lung_int$Ids)

ggplot(lung_int, aes(AA.length))+
  geom_bar()+
  facet_wrap(~Set)+
  geom_vline(xintercept = 9, linetype = 2)+
  labs(title = "Lung immunopeptidome repertoire")

##tumor=====
tumor_int <- enframe(
                    list(
                    "Intersection" = intersect(unique(tumor_chunks_ids),
                                               unique(tumor_digested_ids)),
                    "Chunks only" = setdiff(unique(tumor_chunks_ids),
                                            unique(tumor_digested_ids)),
                    "Digested only" = setdiff(unique(tumor_digested_ids),
                                              unique(tumor_chunks_ids))
                  ),
                  "Set", "Ids"
)
tumor_int <- as.data.frame(unnest_longer(tumor_int, Ids))
tumor_int$AA.length <- str_length(tumor_int$Ids)

ggplot(tumor_int, aes(AA.length))+
  geom_bar()+
  facet_wrap(~Set)+
  geom_vline(xintercept = 9, linetype = 2)+
  labs(title = "Tumor immunopeptidome repertoire")

#Normalization------------------------------------------------------------------
dfWLogN <- quantile_normalize_dm(as.matrix(dfWLog[5:8]))
#dfWLogN <- normalize_sample_medians_dm(as.matrix(dfWLog[5:8]))
dfWLogN <- bind_cols(dfWLog[1:4], dfWLogN)
boxplot(dfWLogN[5:8])
boxplot(dfWLog[5:8])

#Heatmap for clustering---------------------------------------------------------
#filter to at least 2 values in a row
hmT <- dfWLogN
hmT <- replace(dfWLogN, is.na(dfWLogN), 0)
hmT$index <- paste(dfWLogN$Genes, 1:nrow(dfWLogN), sep = ".")
#hmT <- drop_na(dfWLogN)
#Standardize features (proteins)
hmTbig <- drop_na(hmT)
hmTbig <- rowScales(hmTbig[5:8])
Heatmap(as.matrix(drop_na(hmTbig)),
        show_row_dend = FALSE, 
        heatmap_legend_param = list(title = "Row z-score"),
        clustering_distance_columns = "pearson", clustering_method_columns = "ward"
)

#FIND COLLAGEN PEPTIDES IN HM AND COMPARE
hmTsmall <- filter(hmT, grepl("COL", hmT$index))
rownames(hmTsmall) <- hmTsmall$index
#Standardize features (proteins)
hmTsmall <- rowScales(hmTsmall[5:8])
hmTsmall <- drop_na(hmTsmall)
Heatmap(as.matrix(hmTsmall),
        show_row_dend = FALSE, 
        heatmap_legend_param = list(title = "Row z-score"),
        clustering_distance_columns = "pearson", clustering_method_columns = "ward",
        clustering_distance_rows = "pearson", 
        row_names_gp = gpar(fontsize = 7)
)
#Correlation--------------------------------------------------------------------
corr <- cor(as.matrix(drop_na(dfWLogN[5:8])), method = "pearson")
Heatmap(corr)

#PCA----------------------------------------------------------------------------
#filter out NAs
dfWF <- drop_na(dfWLogN, 5:8)
#column to rownames for pca function
dfWF <- column_to_rownames(hmT, "index")
#pca
p <- pca(dfWF[5:8])
biplot(p, lab = colnames(dfWF[5:8]), labSize = 5, showLoadings = TRUE)
screeplot(p)
pairsplot(p, components = getComponents(p, seq_len(4)), 
          pointSize = 3, 
          lab = colnames(dfWF[5:8]),
          labSize = 3)
plotloadings(p, components = getComponents(p, seq_len(4)), drawConnectors = TRUE)
#Linear regression model--------------------------------------------------------
t <- ggplot(dfWLogN, aes(`Tumor chunks`, `Tumor digested`))+
  geom_point(alpha=1/10)+
  stat_smooth(method = "lm")+
  stat_cor(method = "pearson")

l <- ggplot(dfWLogN, aes(`Lung chunks`, `Lung digested`))+
  geom_point(alpha=1/10)+
  stat_smooth(method = "lm")

c <- ggplot(dfWLogN, aes(`Tumor chunks`, `Lung chunks`))+
  geom_point(alpha=1/10)+
  stat_smooth(method = "lm")

d <- ggplot(dfWLogN, aes(`Tumor digested`, `Lung digested`))+
  geom_point(alpha=1/10)+
  stat_smooth(method = "lm")

ggarrange(t, l, c, d)

