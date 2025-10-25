## RNA seq analysis of chemokine stimulation of DC 
#### Samples seperated by chemokine and stimulation for time point 24 h

# script to perform differential gene expression analysis using DESeq2 package

### install missing packages
# Note: It's good practice to install all packages at once or check if they are installed
# before loading them.

BiocManager::install("org.Mm.eg.db")
install.packages("xlsx")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("xlsx", force = TRUE)
BiocManager::install("apeglm")
BiocManager::install("KEGGREST", force = T)
BiocManager::install("clusterProfiler", force = T)
BiocManager::install("GSVA", force = T)
BiocManager::install("GSEABase")
BiocManager::install("DESeq2")
BiocManager::install("biomaRt")
install.packages("VennDiagram") 
BiocManager::install("ashr")
install.packages("stringr")

# load libraries

library(DESeq2)
library(tidyverse)
library(AnnotationDbi)
library(ggplot2)
library(xlsx)
library(org.Mm.eg.db)
library(biomaRt)
library(EnsDb.Mmusculus.v79)
library(RColorBrewer)
library(apeglm)
library(pheatmap)
library(clusterProfiler)
library(KEGGREST)
library(GSVA)
library(GSEABase)
library(VennDiagram) 
library(ashr)
library(edgeR)
library(msigdbr)
library(dplyr)
library(stringr)
library(tidyr)



###########################---------------Start of analysis --------------------------------

##Set directory

setwd("C:/Users/nicol/Desktop/analysis")


## Step 1:  prepare count data -------------------------------------------------------------

## read in counts data 
counts_data <- as.matrix(read.csv('counts.csv', sep = ";", row.names = 1))


##Subset for time point 24h 
counts_data <- counts_data[, c("DE97NGSUKBR136749", "DE70NGSUKBR136750", "DE43NGSUKBR136751",
                               "DE16NGSUKBR136752", "DE86NGSUKBR136753", "DE59NGSUKBR136754",
                               "DE32NGSUKBR136755", "DE05NGSUKBR136756", "DE75NGSUKBR136757", 
                               "DE48NGSUKBR136758", "DE21NGSUKBR136759", "DE91NGSUKBR136760",
                               "DE64NGSUKBR136761", "DE37NGSUKBR136762", "DE10NGSUKBR136763",
                               "DE80NGSUKBR136764", "DE53NGSUKBR136765", "DE26NGSUKBR136766")]


head(counts_data)

## read in sample info for time point 24h 

colData <- read.csv('sample_info_24h.csv', sep = ";", row.names = 1)
head(colData)
str(colData)

### making sure that row names in colData matches to colum names in counts-data
all(colnames(counts_data) %in% rownames(colData))

##are they in the same order?
all(colnames(counts_data) == rownames(colData))

#for reorganization
counts_data <- counts_data[, rownames(colData)]
all(colnames(counts_data) == rownames(colData))


# Set factor levels and create the new `group` factor
colData$stimulation_s <- factor(colData$stimulation_s, levels = c("CpG_IFN", "CpG", "IFN"))
colData$chemokine_s <- factor(colData$chemokine_s, levels = c("no", "yes"))

# Create a combined group factor
colData$group <- factor(paste(colData$stimulation_s, colData$chemokine_s, sep = "_"))
table(colData$group) # Check the new group levels
model.matrix(~group, data = colData)




##Step 2: construct a DESeqDataSet object--------------------------------------------


# Use the new `group` factor in the design formula
dds <- DESeqDataSetFromMatrix(countData = round(counts_data),
                              colData = colData,
                              design = ~group)


dds

### pre-filtering: removing rows with low gene counts 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds


### Run DESeq----------------------------------------------------------------------
dds <- DESeq(dds)
colData(dds)

### Take a look at the resultsNames----------
resultsNames(dds)

sizeFactors(dds)

##############################################################################
################################################################################

########## Check QC-Parameters 

### Boxplots for Raw Counts ###
# --- Get the counts (Raw and Normalized) and prepare them for plotting ---
# Add Group information from colData to the counts data frame
counts_df_raw <- as.data.frame(counts(dds, normalized = FALSE)) %>%
  tibble::rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "count") %>%
  left_join(as.data.frame(colData(dds)) %>% 
              tibble::rownames_to_column(var = "sample"), by = "sample")

counts_df_normalized <- as.data.frame(counts(dds, normalized = TRUE)) %>%
  tibble::rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "normalized_count") %>%
  left_join(as.data.frame(colData(dds)) %>% 
              tibble::rownames_to_column(var = "sample"), by = "sample")

# --- Step 1: Create Boxplot for RAW Counts ---
p_raw <- ggplot(counts_df_raw, aes(x = group, y = log2(count + 1))) +
  geom_boxplot() +
  labs(
    title = "Raw Counts Boxplot",
    x = "Group",
    y = "log2(Count + 1)"
  ) +
  # Using coord_cartesian() to set the y-axis boundaries without cutting data
  coord_cartesian(ylim = c(0, 20)) + 
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16)
  )

# --- Step 2: Create Boxplot for NORMALIZED Counts ---
p_normalized <- ggplot(counts_df_normalized, aes(x = group, y = log2(normalized_count + 1))) +
  geom_boxplot() +
  labs(
    title = "Normalised Counts Boxplot",
    x = "Group",
    y = "log2(Normalised Count + 1)"
  ) +
  # Using coord_cartesian() to set the y-axis boundaries without cutting data
  coord_cartesian(ylim = c(0, 20)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# --- Step 3: Arrange and Save Plots using patchwork ---
# This package lets you combine multiple ggplot objects easily.
library(patchwork)

combined_plot <- p_raw + p_normalized

ggsave("combined_boxplots_24h.png", plot = combined_plot, width = 10, height = 6)

print("Combined boxplots saved as combined_boxplots.png")


# Calculate medians for the raw counts
raw_medians <- counts_df_raw %>%
  group_by(group) %>%
  summarise(median_log2_count = median(log2(count + 1)))

# Calculate medians for the normalized counts
normalized_medians <- counts_df_normalized %>%
  group_by(group) %>%
  summarise(median_log2_normalized_count = median(log2(normalized_count + 1)))

# Print the results
print("Raw Counts Medians:")
print(raw_medians)

print("Normalized Counts Medians:")
print(normalized_medians)

# Calculate the range of the medians to quantify the difference
range_raw <- max(raw_medians$median_log2_count) - min(raw_medians$median_log2_count)
range_normalized <- max(normalized_medians$median_log2_normalized_count) - min(normalized_medians$median_log2_normalized_count)

print(paste("Range of Raw Medians:", round(range_raw, 4)))
print(paste("Range of Normalized Medians:", round(range_normalized, 4)))






##### Dispersion plot

png("DESeq2_dispersion_plot_24h.png", width = 400, height = 300)
plotDispEsts(dds, cex.lab = 1.5, cex.axis = 1.2)
dev.off()


##############################################################################
##############################################################################

########### Sample-to-sample distance heatmap 

# --- Step 1: Transform the data ---
# This is crucial for a distance heatmap. We'll use the Variance Stabilizing Transformation (VST).
# VST stabilizes the variance across the mean expression levels, making the data suitable for clustering.
vst_data <- vst(dds)
sample_dist_matrix <- dist(t(assay(vst_data)))

# --- Step 2: Create a data frame for heatmap annotations ---
# This will add your experimental groups as colored labels on the heatmap.
sample_info <- as.data.frame(colData(dds))
row.names(sample_info) <- colnames(vst_data)

# --- Step 3: Generate and save the heatmap ---
png("sample_distance_heatmap_24h.png", width = 800, height = 600)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dist_matrix,
  clustering_distance_cols = sample_dist_matrix,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_col = sample_info,
  main = "Sample-to-Sample Distance Heatmap"
)
dev.off()




################################################################################################
################################################################################################


##### Step 2: Investigate all data to get an overview -------------------------------------


#### PCA plot to investigate all data ------------------------------
## variance stabilizing transformation for all data 

##biased approach, if many of genes have large differences in counts due to the experimental 
###design, set blind=F for downstream analysis
#vsd <- vst(dds, blind= F) ###-> better use unbaised approach if no big different


# Use the combined group for PCA 
vsd <- vst(dds, blind=T)  
pcaData <- plotPCA(vsd, intgroup ="group", returnData = T)

# Extract variance explained
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot PCA
PCA <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_discrete(
    labels = c("CpG+IFN: - chemokine", "CpG+IFN: + chemokine", "CpG: - chemokine", "CpG: + chemokine", "IFN: - chemokine", "IFN: + chemokine")
  ) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Group"  # Legend title
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightgrey", color = NA),  # Light grey background
    legend.title = element_text(size = 14),                          # Increase legend title size
    legend.text = element_text(size = 12),                           # Increase legend labels size
    axis.title.x = element_text(size = 16),                          # Increase x-axis title size
    axis.title.y = element_text(size = 16),                          # Increase y-axis title size
    axis.text.x = element_text(size = 12),                           # Increase x-axis text size
    axis.text.y = element_text(size = 12),                           # Increase y-axis text size
    plot.title = element_text(size = 18, face = "bold"),              # Bold and larger title
    plot.background = element_rect(fill = "white", color = NA)    # Non-transparent plot area
  )

# Save the final plot
ggsave("PCA_plot_24h.png", plot = PCA, width = 8, height = 4, dpi = 300)


######### Gene expression of key genes associated with DC maturation ----------------------------------

# plot counts
count_CD11c <- plotCounts(dds, gene = "ENSMUSG00000030789", main = "CD11c", intgroup= c("chemokine_s","stimulation_s"), returnData = T) ## CD11c
count_MHCII <- plotCounts(dds, gene = "ENSMUSG00000079547", main = "MHCII", intgroup=c("chemokine_s","stimulation_s"), returnData = T) ## MHCII
count_CD86 <- plotCounts(dds, gene = "ENSMUSG00000022901", main = "CD86", intgroup=c("chemokine_s","stimulation_s"), returnData = T) ## CD86
count_CD80 <- plotCounts(dds, gene = "ENSMUSG00000075122", main = "CD80", intgroup=c("chemokine_s","stimulation_s"), returnData = T) ## CD80
count_CCR7 <- plotCounts(dds, gene = "ENSMUSG00000037944",main = "CCR7", intgroup=c("chemokine_s","stimulation_s"), returnData = T) ## CCR7


head(count_CD11c)

####Figure plotcount_CD11c 

count_CD11c <- as.data.frame(count_CD11c)

# Replace 'CpG_IFN' with 'CpG+IFN' in the factor levels
count_CD11c$stimulation_s <- factor(
  count_CD11c$stimulation_s,
  levels = c("CpG_IFN", "CpG", "IFN"),  # Original levels
  labels = c("CpG+IFN", "CpG", "IFN"))   # Renamed levels


CD11c <- ggplot(count_CD11c, aes(x = stimulation_s, y = count, fill = chemokine_s)) +
  geom_jitter(shape = 21, color = "black", size = 3, position = position_dodge(width = 0.8)) +
  geom_boxplot(alpha = 0.5,
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               width = 0.6
  ) +
  theme_minimal() +
  labs(
    title = "Normalised gene expression CD11c",
    x = "Stimulation",
    y = "Normalised Counts",
    fill = "CCL21"  # Change legend name
  ) +
  theme(
    legend.title = element_text(size = 14),           # Increase legend title size
    legend.text = element_text(size = 12),            # Increase legend text size
    axis.title.x = element_text(size = 14),           # Increase x-axis title size
    axis.title.y = element_text(size = 14),           # Increase y-axis title size
    axis.text = element_text(size = 12),              # Increase x/y-axis label text size
    plot.title = element_text(size = 16, face = "bold"), # Increase plot title size
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),  # Non-transparent background
    plot.background = element_rect(fill = "white", color = NA)    # Non-transparent plot area
  )+
  # Add this function to control the y-axis
  scale_y_continuous(
    limits = c(50, 500),      # Set the minimum and maximum y-axis values
    breaks = seq(50, 500, by = 100)  # Set the step width for the ticks
  )

# Save the final plot
ggsave("CD11c_plot_24h.png", plot = CD11c, width = 6, height = 4, dpi = 300)


####Figure plotcount_MHCII 

count_MHCII <- as.data.frame(count_MHCII)

# Replace 'CpG_IFN' with 'CpG+IFN' in the factor levels
count_MHCII$stimulation_s <- factor(
  count_MHCII$stimulation_s,
  levels = c("CpG_IFN", "CpG", "IFN"),  # Original levels
  labels = c("CpG+IFN", "CpG", "IFN"))   # Renamed levels


MHCII <- ggplot(count_MHCII, aes(x = stimulation_s, y = count, fill = chemokine_s)) +
  geom_jitter(shape = 21, color = "black", size = 3, position = position_dodge(width = 0.8)) +
  geom_boxplot(alpha = 0.5,
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               width = 0.6
  ) +
  theme_minimal() +
  labs(
    title = "Normalised gene expression MHCII",
    x = "Stimulation",
    y = "Normalised Counts",
    fill = "CCL21"  # Change legend name
  ) +
  theme(
    legend.title = element_text(size = 14),           # Increase legend title size
    legend.text = element_text(size = 12),            # Increase legend text size
    axis.title.x = element_text(size = 14),           # Increase x-axis title size
    axis.title.y = element_text(size = 14),           # Increase y-axis title size
    axis.text = element_text(size = 12),              # Increase x/y-axis label text size
    plot.title = element_text(size = 16, face = "bold"), # Increase plot title size
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),  # Non-transparent background
    plot.background = element_rect(fill = "white", color = NA)    # Non-transparent plot area
  )+
  # Add this function to control the y-axis
  scale_y_continuous(
    limits = c(300, 2100),      # Set the minimum and maximum y-axis values
    breaks = seq(300, 2100, by = 200)  # Set the step width for the ticks
  )

# Save the final plot
ggsave("MHCII_plot_24h.png", plot = MHCII, width = 6, height = 4, dpi = 300)



####Figure plotcount_CD80

count_CD80 <- as.data.frame(count_CD80)

# Replace 'CpG_IFN' with 'CpG+IFN' in the factor levels
count_CD80$stimulation_s <- factor(
  count_CD80$stimulation_s,
  levels = c("CpG_IFN", "CpG", "IFN"),  # Original levels
  labels = c("CpG+IFN", "CpG", "IFN"))   # Renamed levels


CD80 <- ggplot(count_CD80, aes(x = stimulation_s, y = count, fill = chemokine_s)) +
  geom_jitter(shape = 21, color = "black", size = 3, position = position_dodge(width = 0.8)) +
  geom_boxplot(alpha = 0.5,
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               width = 0.6
  ) +
  theme_minimal() +
  labs(
    title = "Normalised gene expression CD80",
    x = "Stimulation",
    y = "Normalised Counts",
    fill = "CCL21"  # Change legend name
  ) +
  theme(
    legend.title = element_text(size = 14),           # Increase legend title size
    legend.text = element_text(size = 12),            # Increase legend text size
    axis.title.x = element_text(size = 14),           # Increase x-axis title size
    axis.title.y = element_text(size = 14),           # Increase y-axis title size
    axis.text = element_text(size = 12),              # Increase x/y-axis label text size
    plot.title = element_text(size = 16, face = "bold"), # Increase plot title size
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),  # Non-transparent background
    plot.background = element_rect(fill = "white", color = NA)    # Non-transparent plot area
  )+
  # Add this function to control the y-axis
  scale_y_continuous(
    limits = c(0, 1000),      # Set the minimum and maximum y-axis values
    breaks = seq(0, 1000, by = 200)  # Set the step width for the ticks
  )


# Save the final plot
ggsave("CD80_plot_24h.png", plot = CD80, width = 6, height = 4, dpi = 300)



####Figure plotcount_CD86 

count_CD86 <- as.data.frame(count_CD86)

# Replace 'CpG_IFN' with 'CpG+IFN' in the factor levels
count_CD86$stimulation_s <- factor(
  count_CD86$stimulation_s,
  levels = c("CpG_IFN", "CpG", "IFN"),  # Original levels
  labels = c("CpG+IFN", "CpG", "IFN"))   # Renamed levels


CD86 <- ggplot(count_CD86, aes(x = stimulation_s, y = count, fill = chemokine_s)) +
  geom_jitter(shape = 21, color = "black", size = 3, position = position_dodge(width = 0.8)) +
  geom_boxplot(alpha = 0.5,
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               width = 0.6
  ) +
  theme_minimal() +
  labs(
    title = "Normalised gene expression CD86",
    x = "Stimulation",
    y = "Normalised Counts",
    fill = "CCL21"  # Change legend name
  ) +
  theme(
    legend.title = element_text(size = 14),           # Increase legend title size
    legend.text = element_text(size = 12),            # Increase legend text size
    axis.title.x = element_text(size = 14),           # Increase x-axis title size
    axis.title.y = element_text(size = 14),           # Increase y-axis title size
    axis.text = element_text(size = 12),              # Increase x/y-axis label text size
    plot.title = element_text(size = 16, face = "bold"), # Increase plot title size
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),  # Non-transparent background
    plot.background = element_rect(fill = "white", color = NA)    # Non-transparent plot area
  )+
  # Add this function to control the y-axis
  scale_y_continuous(
    limits = c(100, 1800),      # Set the minimum and maximum y-axis values
    breaks = seq(100, 1800, by = 400)  # Set the step width for the ticks
  )


# Save the final plot
ggsave("CD86_plot_24h.png", plot = CD86, width = 6, height = 4, dpi = 300)



####Figure plotcount_CCR7 

count_CCR7 <- as.data.frame(count_CCR7)

# Replace 'CpG_IFN' with 'CpG+IFN' in the factor levels
count_CCR7$stimulation_s <- factor(
  count_CCR7$stimulation_s,
  levels = c("CpG_IFN", "CpG", "IFN"),  # Original levels
  labels = c("CpG+IFN", "CpG", "IFN"))   # Renamed levels


CCR7 <- ggplot(count_CCR7, aes(x = stimulation_s, y = count, fill = chemokine_s)) +
  geom_jitter(shape = 21, color = "black", size = 3, position = position_dodge(width = 0.8)) +
  geom_boxplot(alpha = 0.5,
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               width = 0.6
  ) +
  theme_minimal() +
  labs(
    title = "Normalised gene expression CCR7",
    x = "Stimulation",
    y = "Normalised Counts",
    fill = "CCL21"  # Change legend name
  ) +
  theme(
    legend.title = element_text(size = 14),           # Increase legend title size
    legend.text = element_text(size = 12),            # Increase legend text size
    axis.title.x = element_text(size = 14),           # Increase x-axis title size
    axis.title.y = element_text(size = 14),           # Increase y-axis title size
    axis.text = element_text(size = 12),              # Increase x/y-axis label text size
    plot.title = element_text(size = 16, face = "bold"), # Increase plot title size
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),  # Non-transparent background
    plot.background = element_rect(fill = "white", color = NA)    # Non-transparent plot area
  )+
  # Add this function to control the y-axis
  scale_y_continuous(
    limits = c(2000, 12200),      # Set the minimum and maximum y-axis values
    breaks = seq(200, 12200, by = 2000)  # Set the step width for the ticks
  )

# Save the final plot
ggsave("CCR7_plot_24h.png", plot = CCR7, width = 6, height = 4, dpi = 300)





###################################################################################################
###################################################################################################

### Stp3: More specific look into DESeq Data for each maturation ----------------------------------------

# Step 1: Extract results for each group and asign gene symbols to Ensembl numbers --------------------------------------------- 

# Use the contrast argument with the new 'group' factor
# 1. CpG_yes vs CpG_no
results_cpg <- results(dds, contrast = c("group", "CpG_yes", "CpG_no"))
summary(results_cpg)
head(results_cpg)


# 2. IFN_yes vs IFN_no
results_ifn <- results(dds, contrast = c("group", "IFN_yes", "IFN_no"))
summary(results_ifn)
head(results_ifn)

# 3. CpG_IFN_yes vs CpG_IFN_no
results_cpg_ifn <- results(dds, contrast = c("group", "CpG_IFN_yes", "CpG_IFN_no"))
summary(results_cpg_ifn)
head(results_cpg_ifn)






# Add gene symbols to each result object (as before)
results_cpg$symbol <- mapIds(org.Mm.eg.db, 
                             keys = rownames(results_cpg), 
                             column = "SYMBOL", 
                             keytype = "ENSEMBL", 
                             multiVals = "first")

results_ifn$symbol <- mapIds(org.Mm.eg.db, 
                             keys = rownames(results_ifn), 
                             column = "SYMBOL", 
                             keytype = "ENSEMBL", 
                             multiVals = "first")


results_cpg_ifn$symbol <- mapIds(org.Mm.eg.db, 
                                 keys = rownames(results_cpg_ifn), 
                                 column = "SYMBOL", 
                                 keytype = "ENSEMBL", 
                                 multiVals = "first")



##### save tables containing all data for the respective condition
write.csv(results_cpg, file = "CpG-Chemokine-no-yes_24h_result-all.csv", row.names = F)
write.csv(results_ifn, file = "IFN-Chemokine-no-yes_24h_result-all.csv", row.names = F)
write.csv(results_cpg_ifn, file = "CpG-IFN-Chemokine-no-yes_24h_result-all.csv", row.names = F)


##################### Create histogram of p-values to get an overview of distribution of padj-values--------------------------------
# For CpG (using `results_cpg` with the new contrast)
png("Histogramm-CpG_24h-padj.png", width = 650, height = 450)
par(mar = c(5.1, 6.1, 5.1, 2.1))
hist(results_cpg$padj, breaks = seq(0,1, length=21), col= "grey", border = "black",
     xlab = "Adjusted p-value (padj)", ylab = "Frequency", cex.lab = 2,
     ylim= c(0,15000), main = "Frequencies of padj-values - CpG", cex.main = 2.5, xaxt = "n", cex.axis = 1.5)
axis(side = 1, at = seq(0, 1, by = 0.1,), cex.axis = 1.5)
dev.off()

# For IFN (using `results_ifn` with the new contrast)
png("Histogramm-IFN_24h-padj.png", width = 650, height = 450)
par(mar = c(5.1, 6.1, 5.1, 2.1))
hist(results_ifn$padj, breaks = seq(0,1, length=21), col= "grey", border = "black",
     xlab = "Adjusted p-value (padj)", ylab = "Frequency", cex.lab = 2,
     ylim= c(0,15000), main = "Frequencies of padj-values - IFN", cex.main =2.5, xaxt = "n", cex.axis = 1.5)
axis(side = 1, at = seq(0, 1, by = 0.1), cex.axis = 1.5)
dev.off()

# For CpG+IFN (using `results_cpg_ifn` with the new contrast)
png("Histogramm-CpG_IFN_24h-padj.png", width = 650, height = 450)
par(mar = c(5.1, 6.1, 5.1, 2.1))
hist(results_cpg_ifn$padj, breaks = seq(0,1, length=21), col= "grey", border = "black",
     xlab = "Adjusted p-value (padj)", ylab = "Frequency", cex.lab = 2,
     ylim= c(0,15000), main = "Frequencies of padj-values - CpG+IFN", cex.main =2.5, xaxt = "n", cex.axis = 1.5)
axis(side = 1, at = seq(0, 1, by = 0.1), cex.axis = 1.5)
dev.off()


####################################################################################################
####################################################################################################

###### Extract significant genes for each group you want to compare -----------------------------
###### Corrected Script for Heatmap Generation ---------------------------

# Step 1: Generate heatmaps -----------------------------------------------

# Filter Significant Results (function remains the same)
filter_significant <- function(results_obj, padj_threshold = 0.1, log2fc_threshold = 1) {
  results_obj <- results_obj[!is.na(results_obj$symbol), ]  # remove rows with NA symbols
  results_obj <- results_obj[!is.na(results_obj$padj), ] # Remove rows with NA in padj
  significant_results <- results_obj[
    results_obj$padj < padj_threshold & abs(results_obj$log2FoldChange) > log2fc_threshold, 
  ]
  rownames(significant_results) <- significant_results$symbol
  return(significant_results)
}

# Apply the filter to each set of results (from the new contrasts)
significant_cpg <- filter_significant(results_cpg)
significant_ifn <- filter_significant(results_ifn)
significant_cpg_ifn <- filter_significant(results_cpg_ifn)

# Check dimensions
dim(significant_cpg)
dim(significant_ifn)
dim(significant_cpg_ifn)

# Save significant genes in csv files (with unique names)
write.csv(significant_cpg, file = "CpG-chemokine_no_yes_sig_24h.csv", row.names = F)
write.csv(significant_ifn, file = "IFN-chemokine_no_yes_sig_24h.csv", row.names = F)
write.csv(significant_cpg_ifn, file = "CpG-IFN-chemokine_no_yes_sig_24h.csv", row.names = F)

# Step 3: Extract Normalized Counts and Fix NA Rownames
### Final Corrected Script for Heatmap Generation

normalized_counts <- counts(dds, normalized = TRUE)
rownames(normalized_counts) <- sub("\\..*", "", rownames(normalized_counts))

gene_symbols <- mapIds(
  org.Mm.eg.db,
  keys = rownames(normalized_counts),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)
# Remove rows where a gene symbol could not be found
na_rows <- is.na(gene_symbols)
normalized_counts <- normalized_counts[!na_rows, ]
gene_symbols <- gene_symbols[!na_rows]
rownames(normalized_counts) <- gene_symbols


# --- Step 3: Extract and Filter for All Significant Genes ---


# Get all significant genes based on the filtered results table
all_significant_cpg_genes <- rownames(significant_cpg)
all_significant_ifn_genes <- rownames(significant_ifn)
all_significant_cpg_ifn_genes <- rownames(significant_cpg_ifn)

# Subset the normalized counts matrix to include only these genes
log_counts_cpg <- log2(normalized_counts[all_significant_cpg_genes, colData(dds)$group %in% c("CpG_yes", "CpG_no")] + 1)
log_counts_ifn <- log2(normalized_counts[all_significant_ifn_genes, colData(dds)$group %in% c("IFN_yes", "IFN_no")] + 1)
log_counts_cpg_ifn <- log2(normalized_counts[all_significant_cpg_ifn_genes, colData(dds)$group %in% c("CpG_IFN_yes", "CpG_IFN_no")] + 1)

# Prepare annotation data frame
annot_cpg <- as.data.frame(colData(dds)[colnames(log_counts_cpg), c("group"), drop = FALSE])
annot_ifn <- as.data.frame(colData(dds)[colnames(log_counts_ifn), c("group"), drop = FALSE])
annot_cpg_ifn <- as.data.frame(colData(dds)[colnames(log_counts_cpg_ifn), c("group"), drop = FALSE])

# Drop unused factor levels
annot_cpg$group <- droplevels(annot_cpg$group)
annot_ifn$group <- droplevels(annot_ifn$group)
annot_cpg_ifn$group <- droplevels(annot_cpg_ifn$group)


# Now, rename the column and the levels
annot_cpg <- annot_cpg %>% dplyr::rename(chemokine = group)
levels(annot_cpg$chemokine) <- c("-", "+")
annot_ifn <- annot_ifn %>% dplyr::rename(chemokine = group)
levels(annot_ifn$chemokine) <- c("-", "+")
annot_cpg_ifn <- annot_cpg_ifn %>% dplyr::rename(chemokine = group)
levels(annot_cpg_ifn$chemokine) <- c("-", "+")


# Reorder columns by group for a clean heatmap layout
reordered_cpg <- log_counts_cpg[, order(annot_cpg$chemokine)]
annot_cpg <- annot_cpg[order(annot_cpg$chemokine), , drop = FALSE]

reordered_ifn <- log_counts_ifn[, order(annot_ifn$chemokine)]
annot_ifn <- annot_ifn[order(annot_ifn$chemokine), , drop = FALSE]

reordered_cpg_ifn <- log_counts_cpg_ifn[, order(annot_cpg_ifn$chemokine)]
annot_cpg_ifn <- annot_cpg_ifn[order(annot_cpg_ifn$chemokine), , drop = FALSE]


# Generate heatmap
# Note: The height of the plot is increased to accommodate all 91 genes
png("Heatmap_all_significant_CpG.png", width = 300, height = 1200) 
pheatmap(
  mat = reordered_cpg,
  scale = "row",
  annotation_col = annot_cpg,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "All Significant Genes - CpG",
  show_colnames = FALSE,
  border_color = "darkgrey"
)
dev.off()


png("Heatmap_all_significant_IFN.png", width = 300, height = 1500) 
pheatmap(
  mat = reordered_ifn,
  scale = "row",
  annotation_col = annot_ifn,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "All Significant Genes - IFN",
  show_colnames = FALSE,
  border_color = "darkgrey"
)
dev.off()

png("Heatmap_all_significant_CpG_IFN.png", width = 300, height = 1700) 
pheatmap(
  mat = reordered_cpg_ifn,
  scale = "row",
  annotation_col = annot_cpg_ifn,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "All Significant Genes - CpG+IFN",
  show_colnames = FALSE,
  border_color = "darkgrey"
)
dev.off()



# Step 3: Extract and Filter for Top 50 Genes
# --- This is the key adapted section ---

# Subset data for each contrast, now keeping only the top 50 significant genes
top50_cpg_genes <- head(rownames(significant_cpg), 50)
top50_ifn_genes <- head(rownames(significant_ifn), 50)
top50_cpg_ifn_genes <- head(rownames(significant_cpg_ifn), 50)



# Subset the normalized counts matrix to include only these 5 genes
# The 'normalized_counts' object is the correct one to use at this point in your script
log_counts_cpg <- log2(normalized_counts[top50_cpg_genes, colData(dds)$group %in% c("CpG_yes", "CpG_no")] + 1)
log_counts_ifn <- log2(normalized_counts[top50_ifn_genes, colData(dds)$group %in% c("IFN_yes", "IFN_no")] + 1)
log_counts_cpg_ifn <- log2(normalized_counts[top50_cpg_ifn_genes, colData(dds)$group %in% c("CpG_IFN_yes", "CpG_IFN_no")] + 1)

# CpG Annotation
annot_cpg <- as.data.frame(colData(dds)[colnames(log_counts_cpg), "group", drop = FALSE])
annot_cpg$group <- droplevels(annot_cpg$group)
annot_cpg <- annot_cpg %>% dplyr::rename(chemokine = group)
levels(annot_cpg$chemokine) <- c("-", "+")

# IFN Annotation
annot_ifn <- as.data.frame(colData(dds)[colnames(log_counts_ifn), "group", drop = FALSE])
annot_ifn$group <- droplevels(annot_ifn$group)
annot_ifn <- annot_ifn %>% dplyr::rename(chemokine = group)
levels(annot_ifn$chemokine) <- c("-", "+")

# CpG+IFN Annotation
annot_cpg_ifn <- as.data.frame(colData(dds)[colnames(log_counts_cpg_ifn), "group", drop = FALSE])
annot_cpg_ifn$group <- droplevels(annot_cpg_ifn$group)
annot_cpg_ifn <- annot_cpg_ifn %>% dplyr::rename(chemokine = group)
levels(annot_cpg_ifn$chemokine) <- c("-", "+")



# Reorder columns by group for a clean heatmap layout (as before)
reordered_cpg <- log_counts_cpg[, order(annot_cpg$chemokine)]
annot_cpg <- annot_cpg[order(annot_cpg$chemokine), , drop = FALSE]

reordered_ifn <- log_counts_ifn[, order(annot_ifn$chemokine)]
annot_ifn <- annot_ifn[order(annot_ifn$chemokine), , drop = FALSE]

reordered_cpg_ifn <- log_counts_cpg_ifn[, order(annot_cpg_ifn$chemokine)]
annot_cpg_ifn <- annot_cpg_ifn[order(annot_cpg_ifn$chemokine), , drop = FALSE]



# Generate heatmap for CpG contrast with only the top 5 genes
png("Heatmap_top50_24h_CpG.png", width = 300, height = 500) 
pheatmap(
  mat = reordered_cpg,
  scale = "row",
  annotation_col = annot_cpg,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Top 50 Genes - CpG",
  show_colnames = FALSE,
  border_color = "darkgrey"
)
dev.off()


png("Heatmap_top50_24h_IFN.png", width = 300, height = 500) 
pheatmap(
  mat = reordered_ifn,
  scale = "row",
  annotation_col = annot_ifn,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Top 50 Genes - IFN",
  show_colnames = FALSE,
  border_color = "darkgrey"
)
dev.off()



png("Heatmap_top50_24h_CpG+IFN.png", width = 300, height = 500) 
pheatmap(
  mat = reordered_cpg_ifn,
  scale = "row",
  annotation_col = annot_cpg_ifn,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Top 50 Genes - CpG+IFN",
  show_colnames = FALSE,
  border_color = "darkgrey"
)
dev.off()



##############################



################################################################################################
##################################################################################################

##Step 2:
####### Investigate the results/distribution using MA plots------------------------------------------


##### MA-plot without background correction
##### shows log2 fold changes attributable to a given variable over the mean o normalized counts for all
### the samples in the dataset. points will be colored red, if the adjusted p value is less than 0.1



##### MA plot for CpG
png("plotMA-CpG_24h.png", width = 600, height = 400) 
DESeq2::plotMA(results_cpg, ylim= c(-12,8), alpha = 0.05, cex = 0.8, colSig = "red")
dev.off()

##### MA-plot with background correction 

# Get the results object with the contrast you want to analyze
res_cpg <- results(dds, contrast=c("group", "CpG_yes", "CpG_no"))

# Now use lfcShrink() to create the shrunken object, using the 'res' object you just created.
#    The shrunken object 'resLFC_cpg' will be used for the MA plot.
resLFC_cpg <- lfcShrink(dds, contrast=c("group", "CpG_yes", "CpG_no"), res=res_cpg, type="ashr")

# Use the corrected object to create the plot
png("plotMA-CpG_24h-background-corr.png", width = 600, height = 400)
par(mar = c(5, 6, 4, 2) + 0.1)
DESeq2::plotMA(resLFC_cpg, ylim= c(-8,3), alpha= 0.05, cex = 1, colSig = "red",cex.lab = 2, cex.axis = 1.5, main = "MA Plot: CpG ", cex.main = 2)
dev.off()



##### MA plot for IFN

png("plotMA-IFN_24h.png", width = 600, height = 400) 
DESeq2::plotMA(results_ifn, ylim= c(-8,8), alpha = 0.05, cex = 0.8, colSig = "red")
dev.off()


##### MA-plot with background correction 

# Get the results object with the contrast you want to analyze
res_IFN <- results(dds, contrast=c("group", "IFN_yes", "IFN_no"))

# Now use lfcShrink() to create the shrunken object, using the 'res' object you just created.
#    The shrunken object 'resLFC_cpg' will be used for the MA plot.
resLFC_IFN <- lfcShrink(dds, contrast=c("group", "IFN_yes", "IFN_no"), res=res_IFN, type="ashr")

# Use the corrected object to create the plot
png("plotMA-IFN_24h-background-corr.png", width = 600, height = 400)
par(mar = c(5, 6, 4, 2) + 0.1)
DESeq2::plotMA(resLFC_IFN, ylim= c(-8,3), alpha= 0.05, cex = 1, colSig = "red", cex.lab = 2, cex.axis = 1.5, main = "MA Plot: IFN", cex.main = 2)
dev.off()



######## MA plot for CpG+IFN

png("plotMA-CpG_IFN_24h.png", width = 600, height = 400) 
DESeq2::plotMA(results_cpg_ifn, ylim= c(-14,10), alpha = 0.05, cex = 0.8, colSig = "red")
dev.off()

##### MA-plot with background correction 

# Get the results object with the contrast you want to analyze
res_CpG_IFN <- results(dds, contrast=c("group", "CpG_IFN_yes", "CpG_IFN_no"))

# Now use lfcShrink() to create the shrunken object, using the 'res' object you just created.
#    The shrunken object 'resLFC_cpg' will be used for the MA plot.
resLFC_CpG_IFN <- lfcShrink(dds, contrast=c("group", "CpG_IFN_yes", "CpG_IFN_no"), res=res_CpG_IFN, type="ashr")

# Use the corrected object to create the plot
png("plotMA-CpG-IFN_24h-background-corr.png", width = 600, height = 400)
par(mar = c(5, 6, 4, 2) + 0.1)
DESeq2::plotMA(resLFC_CpG_IFN, ylim= c(-8,3), alpha= 0.05, cex = 1, colSig = "red", cex.lab = 2, cex.axis = 1.5, main = "MA Plot: CpG+IFN", cex.main = 2)
dev.off()


############################################################################################
############################################################################################


###### Visualise the enriched genes as Volcano plots



##plot value for CpG

png("Volcano-CpG_24h-padj.png", width = 630, height = 420)
par(mar = c(5, 5, 5, 13))  # Increase right margin (last value)
plot(results_cpg$log2FoldChange, -log10(results_cpg$padj),  
     xlab = "log2 Fold Change",
     ylab = "-log10(Adjusted P-value)",
     pch = 20, cex = 1.2, 
     cex.lab = 2,
     cex.main = 2.5,
     cex.axis = 1.5,
     xlim = c(-30, 40),  # Set x-axis range from -40 to 40
     ylim = c(0, 14)     # Set y-axis range from 0 to 12
)

with(subset(results_cpg, padj < 0.1 & abs(log2FoldChange) >= 1), 
     points(log2FoldChange, -log10(padj),  
            pch = 20, 
            col = c("#00BFFF", "red")[(sign(log2FoldChange) + 3)/2], 
            cex=1.5))

legend("topright", 
       inset = c(-0.45, 0), # This value is key. Adjust it to fine-tune position.
       title = paste("Padj < ", 0.1, sep = ""),
       legend = c("downregulated", "upregulated"),
       pch = 20, col = c("#00BFFF", "red"),
       xpd = TRUE, # This is crucial to allow drawing outside the plot
       cex = 1.5
)
# Use mtext() to add the title and push it to the right
mtext("CpG: no Chemokine vs Chemokine", side = 3, adj = 0.2, line = 2, cex = 2.5, font = 2)

dev.off()


##plot value for IFN

png("Volcano-IFN_24h-padj.png", width = 630, height = 420)
par(mar = c(5, 5, 5, 13))  # Increase right margin (last value)
plot(results_ifn$log2FoldChange, -log10(results_ifn$padj),   
     xlab = "log2 Fold Change",
     ylab = "-log10(Adjusted P-value)",
     pch = 20, cex = 1.2, 
     cex.lab = 2,
     cex.main = 2.5,
     cex.axis = 1.5,
     xlim = c(-30, 40),  # Set x-axis range from -40 to 40
     ylim = c(0, 14)     # Set y-axis range from 0 to 12
)

with(subset(results_ifn, padj < 0.1 & abs(log2FoldChange) >= 1), 
     points(log2FoldChange, -log10(padj),  
            pch = 20,
            col = c("#00BFFF", "red")[(sign(log2FoldChange) + 3)/2], 
            cex=1.5))

legend("topright", inset = c(-0.45, 0),  # Adjust position (inset moves it outside)
       title = paste("Padj < ", 0.1, sep = ""),
       legend = c("downregulated", "upregulated"),
       pch = 20, col = c("#00BFFF", "red"),
       xpd = TRUE,  # Allow drawing outside the plot area
       cex = 1.5
)
# Use mtext() to add the title and push it to the right
mtext("IFN: no Chemokine vs Chemokine", side = 3, adj = 0.2, line = 2, cex = 2.5, font = 2)

dev.off()



##plot value for CpG+IFN

png("Volcano-CpG-IFN_24h-padj.png", width = 630, height = 420)

par(mar = c(5, 5, 5, 13))  # Increase right margin (last value)
plot(results_cpg_ifn$log2FoldChange, -log10(results_cpg_ifn$padj),   
     xlab = "log2 Fold Change",
     ylab = "-log10(Adjusted P-value)",
     pch = 20, cex = 1.2, 
     cex.lab = 2,
     cex.main = 2.5,
     cex.axis = 1.5,
     xlim = c(-30, 40),  # Set x-axis range from -40 to 40
     ylim = c(0, 14)     # Set y-axis range from 0 to 12
)

with(subset(results_cpg_ifn, padj < 0.1 & abs(log2FoldChange) >= 1), 
     points(log2FoldChange, -log10(padj),  
            pch = 20,
            col = c("#00BFFF", "red")[(sign(log2FoldChange) + 3)/2], 
            cex=1.5))

legend("topright", inset = c(-0.45, 0),  # Adjust position (inset moves it outside)
       title = paste("Padj < ", 0.1, sep = ""),
       legend = c("downregulated", "upregulated"),
       pch = 20, col = c("#00BFFF", "red"),
       xpd = TRUE,  # Allow drawing outside the plot area
       cex = 1.5
)
# Use mtext() to add the title and push it to the right
mtext("CpG+IFN: no Chemokine vs Chemokine", side = 3, adj = 0.2, line = 2, cex = 2.5, font = 2)

dev.off()




################################################################################################
################################################################################################

###################### GO Analysis ---------------------------------------------------------------


#### IMPORTANT: Re-establish `results_cpg_original_ensembl`
# For GO analysis, we need the original Ensembl IDs before converting rownames to symbols.
# Assuming 'dds' is your DESeqDataSet object after DESeq() run.
# This results object will keep Ensembl IDs as rownames
results_cpg_ensembl <- results(dds, contrast = c("group", "CpG_yes", "CpG_no"))

# Define the universe of genes (all genes considered in DESeq2, in Ensembl IDs)
# Make sure to remove any NA Ensembl IDs if they somehow exist at this stage
universe_genes_ensembl <- rownames(results_cpg_ensembl)
universe_genes_ensembl <- universe_genes_ensembl[!is.na(universe_genes_ensembl)]

# Define padj and log2FoldChange thresholds
padj_thresh <- 0.1
log2fc_thresh <- 1

# --- Function to perform and plot GO enrichment ---
perform_and_plot_go <- function(gene_list, ontology_type, direction_label, file_suffix, results_obj_ensembl, universe_genes_ensembl) {
  
  # Ensure the gene list passed for enrichment is in Ensembl IDs
  # This step uses the 'results_obj_ensembl' to map symbols back to Ensembl if needed
  # or directly extracts Ensembl IDs if the initial filtering was done on Ensembl IDs
  
  # For your current setup where results_cpg (used for top50 heatmap) has symbols
  # and results_cpg_ensembl has Ensembl, we need to get Ensembl IDs of significant genes
  
  # Subset the results_obj_ensembl (which has Ensembl IDs as rownames)
  # based on the gene symbols found in the 'gene_list' (from results_cpg which had symbols)
  
  # First, ensure results_obj_ensembl has the 'symbol' column
  if (!"symbol" %in% colnames(results_obj_ensembl)) {
    results_obj_ensembl$symbol <- mapIds(org.Mm.eg.db, keys = rownames(results_obj_ensembl), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  }
  
  # Filter significant genes from the Ensembl-ID-based results object
  sig_genes_ensembl <- results_obj_ensembl[!is.na(results_obj_ensembl$padj) & 
                                             !is.na(results_obj_ensembl$symbol) &
                                             results_obj_ensembl$padj < padj_thresh & 
                                             abs(results_obj_ensembl$log2FoldChange) >= log2fc_thresh, ]
  
  if (direction_label == "up regulated") {
    gene_ids_for_enrichment <- rownames(sig_genes_ensembl[sig_genes_ensembl$log2FoldChange >= log2fc_thresh, ])
  } else if (direction_label == "down regulated") {
    gene_ids_for_enrichment <- rownames(sig_genes_ensembl[sig_genes_ensembl$log2FoldChange <= -log2fc_thresh, ])
  } else {
    stop("Invalid direction_label provided.")
  }
  
  if (length(gene_ids_for_enrichment) == 0) {
    message(paste("No", direction_label, "genes found for", ontology_type, ". Skipping enrichment."))
    return(NULL)
  }
  
  go_results <- enrichGO(
    gene = gene_ids_for_enrichment,
    universe = universe_genes_ensembl, # Use the global universe of Ensembl IDs
    OrgDb = org.Mm.eg.db,
    keyType = "ENSEMBL", # Keep as ENSEMBL as universe and gene IDs are Ensembl
    ont = ontology_type,
    pAdjustMethod = "BH", # Benjamini-Hochberg for FDR correction
    qvalueCutoff = padj_thresh,
    readable = TRUE # Converts Ensembl IDs to gene symbols in the result table
  )
  
  if (is.null(go_results) || nrow(go_results@result) == 0) {
    message(paste("No enriched GO terms found for", direction_label, "genes in", ontology_type, "."))
    return(NULL)
  }
  
  # Save result as csv file
  file_name <- paste0("GO-term_CpG_24h_", direction_label, "-", ontology_type, file_suffix, ".csv")
  write.csv(go_results@result, file = file_name, row.names = F)
  message(paste("Saved GO results to", file_name))
  
  go_results@result$Description <- str_wrap(go_results@result$Description, width = 40)
  
  # Plot result
  plot_title <- paste("GO Enrichment of", direction_label, "genes -", ontology_type)
  plot_file_name <- paste0("GO-dotplot_CpG_24h_", direction_label, "-", ontology_type, file_suffix, ".png")
  
  png(plot_file_name, width = 500, height = 250)
  p <- dotplot(go_results, showCategory = 20, font.size = 12, label_format = 70) +
    scale_size_continuous(range = c(1, 7)) +
    theme_minimal() +
    ggtitle(plot_title) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size
      axis.title = element_text(size = 16),              # Increase axis label size
      axis.text = element_text(size = 14),               # Increase axis number size
      legend.text = element_text(size = 12),             # Increase legend text size
      legend.title = element_text(size = 14, face = "bold") # Increase legend title size
    )
  print(p) # Explicitly print plot when inside a function
  dev.off()
  message(paste("Saved GO dotplot to", plot_file_name))
  
  return(go_results)
}

# --- Run GO analysis for positive log2FoldChange (up-regulated genes) ---
message("--- Analyzing up-regulated genes ---")
# Biological Processes (BP)
GO_sig_cpg_pos_bp <- perform_and_plot_go(
  gene_list = NULL, # Placeholder, as gene filtering is done inside the function
  ontology_type = "BP",
  direction_label = "up regulated",
  file_suffix = "_BP",
  results_obj_ensembl = results_cpg_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Molecular Function (MF)
GO_sig_cpg_pos_mf <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "MF",
  direction_label = "up regulated",
  file_suffix = "_MF",
  results_obj_ensembl = results_cpg_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Cellular Component (CC)
GO_sig_cpg_pos_cc <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "CC",
  direction_label = "up regulated",
  file_suffix = "_CC",
  results_obj_ensembl = results_cpg_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)


# --- Run GO analysis for negative log2FoldChange (down-regulated genes) ---
message("--- Analyzing down-regulated genes ---")
# Biological Processes (BP)
GO_sig_cpg_neg_bp <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "BP",
  direction_label = "down regulated",
  file_suffix = "_BP",
  results_obj_ensembl = results_cpg_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Molecular Function (MF)
GO_sig_cpg_neg_mf <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "MF",
  direction_label = "down regulated",
  file_suffix = "_MF",
  results_obj_ensembl = results_cpg_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Cellular Component (CC)
GO_sig_cpg_neg_cc <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "CC",
  direction_label = "down regulated",
  file_suffix = "_CC",
  results_obj_ensembl = results_cpg_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)


############### For IFN ########### 

#### IMPORTANT: Re-establish `results_ifn_original_ensembl`
# For GO analysis, we need the original Ensembl IDs before converting rownames to symbols.
# Assuming 'dds' is your DESeqDataSet object after DESeq() run.
# This results object will keep Ensembl IDs as rownames
results_ifn_ensembl <- results(dds, contrast = c("group", "IFN_yes", "IFN_no"))

# Define the universe of genes (all genes considered in DESeq2, in Ensembl IDs)
# Make sure to remove any NA Ensembl IDs if they somehow exist at this stage
universe_genes_ensembl <- rownames(results_ifn_ensembl)
universe_genes_ensembl <- universe_genes_ensembl[!is.na(universe_genes_ensembl)]

# Define padj and log2FoldChange thresholds
padj_thresh <- 0.1
log2fc_thresh <- 1

# --- Function to perform and plot GO enrichment ---
perform_and_plot_go <- function(gene_list, ontology_type, direction_label, file_suffix, results_obj_ensembl, universe_genes_ensembl) {
  
  # Ensure the gene list passed for enrichment is in Ensembl IDs
  # This step uses the 'results_obj_ensembl' to map symbols back to Ensembl if needed
  # or directly extracts Ensembl IDs if the initial filtering was done on Ensembl IDs
  
  # For your current setup where results_ifn (used for top50 heatmap) has symbols
  # and results_ifn_ensembl has Ensembl, we need to get Ensembl IDs of significant genes
  
  # Subset the results_obj_ensembl (which has Ensembl IDs as rownames)
  # based on the gene symbols found in the 'gene_list' (from results_ifn which had symbols)
  
  # First, ensure results_obj_ensembl has the 'symbol' column
  if (!"symbol" %in% colnames(results_obj_ensembl)) {
    results_obj_ensembl$symbol <- mapIds(org.Mm.eg.db, keys = rownames(results_obj_ensembl), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  }
  
  # Filter significant genes from the Ensembl-ID-based results object
  sig_genes_ensembl <- results_obj_ensembl[!is.na(results_obj_ensembl$padj) & 
                                             !is.na(results_obj_ensembl$symbol) &
                                             results_obj_ensembl$padj < padj_thresh & 
                                             abs(results_obj_ensembl$log2FoldChange) >= log2fc_thresh, ]
  
  if (direction_label == "up regulated") {
    gene_ids_for_enrichment <- rownames(sig_genes_ensembl[sig_genes_ensembl$log2FoldChange >= log2fc_thresh, ])
  } else if (direction_label == "down regulated") {
    gene_ids_for_enrichment <- rownames(sig_genes_ensembl[sig_genes_ensembl$log2FoldChange <= -log2fc_thresh, ])
  } else {
    stop("Invalid direction_label provided.")
  }
  
  if (length(gene_ids_for_enrichment) == 0) {
    message(paste("No", direction_label, "genes found for", ontology_type, ". Skipping enrichment."))
    return(NULL)
  }
  
  go_results <- enrichGO(
    gene = gene_ids_for_enrichment,
    universe = universe_genes_ensembl, # Use the global universe of Ensembl IDs
    OrgDb = org.Mm.eg.db,
    keyType = "ENSEMBL", # Keep as ENSEMBL as universe and gene IDs are Ensembl
    ont = ontology_type,
    pAdjustMethod = "BH", # Benjamini-Hochberg for FDR correction
    qvalueCutoff = padj_thresh,
    readable = TRUE # Converts Ensembl IDs to gene symbols in the result table
  )
  
  if (is.null(go_results) || nrow(go_results@result) == 0) {
    message(paste("No enriched GO terms found for", direction_label, "genes in", ontology_type, "."))
    return(NULL)
  }
  
  # Save result as csv file
  file_name <- paste0("GO-term_IFN_24h_", direction_label, "-", ontology_type, file_suffix, ".csv")
  write.csv(go_results@result, file = file_name, row.names = F)
  message(paste("Saved GO results to", file_name))
  
  go_results@result$Description <- str_wrap(go_results@result$Description, width = 40)
  
  # Plot result
  plot_title <- paste("GO Enrichment of", direction_label, "genes -", ontology_type)
  plot_file_name <- paste0("GO-dotplot_IFN_24h_", direction_label, "-", ontology_type, file_suffix, ".png")
  
  png(plot_file_name, width = 500, height = 250)
  p <- dotplot(go_results, showCategory = 20, font.size = 12, label_format = 20) +
    scale_size_continuous(range = c(1, 7)) +
    theme_minimal() +
    ggtitle(plot_title) +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),  # Increase title size
  axis.title = element_text(size = 16),              # Increase axis label size
  axis.text = element_text(size = 14),               # Increase axis number size
  legend.text = element_text(size = 12),             # Increase legend text size
  legend.title = element_text(size = 14, face = "bold") # Increase legend title size
  )

  print(p) # Explicitly print plot when inside a function
  dev.off()
  message(paste("Saved GO dotplot to", plot_file_name))
  
  return(go_results)
}

# --- Run GO analysis for positive log2FoldChange (up-regulated genes) ---
message("--- Analyzing up-regulated genes ---")
# Biological Processes (BP)
GO_sig_ifn_pos_bp <- perform_and_plot_go(
  gene_list = NULL, # Placeholder, as gene filtering is done inside the function
  ontology_type = "BP",
  direction_label = "up regulated",
  file_suffix = "_BP",
  results_obj_ensembl = results_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Molecular Function (MF)
GO_sig_ifn_pos_mf <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "MF",
  direction_label = "up regulated",
  file_suffix = "_MF",
  results_obj_ensembl = results_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Cellular Component (CC)
GO_sig_ifn_pos_cc <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "CC",
  direction_label = "up regulated",
  file_suffix = "_CC",
  results_obj_ensembl = results_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)


# --- Run GO analysis for negative log2FoldChange (down-regulated genes) ---
message("--- Analyzing down-regulated genes ---")
# Biological Processes (BP)
GO_sig_ifn_neg_bp <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "BP",
  direction_label = "down regulated",
  file_suffix = "_BP",
  results_obj_ensembl = results_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Molecular Function (MF)
GO_sig_ifn_neg_mf <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "MF",
  direction_label = "down regulated",
  file_suffix = "_MF",
  results_obj_ensembl = results_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Cellular Component (CC)
GO_sig_ifn_neg_cc <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "CC",
  direction_label = "down regulated",
  file_suffix = "_CC",
  results_obj_ensembl = results_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)



############### For CpG+IFN ########### 

#### IMPORTANT: Re-establish `results_ifn_original_ensembl`
# For GO analysis, we need the original Ensembl IDs before converting rownames to symbols.
# Assuming 'dds' is your DESeqDataSet object after DESeq() run.
# This results object will keep Ensembl IDs as rownames
results_cpg_ifn_ensembl <- results(dds, contrast = c("group", "CpG_IFN_yes", "CpG_IFN_no"))

# Define the universe of genes (all genes considered in DESeq2, in Ensembl IDs)
# Make sure to remove any NA Ensembl IDs if they somehow exist at this stage
universe_genes_ensembl <- rownames(results_cpg_ifn_ensembl)
universe_genes_ensembl <- universe_genes_ensembl[!is.na(universe_genes_ensembl)]

# Define padj and log2FoldChange thresholds
padj_thresh <- 0.1
log2fc_thresh <- 1

# --- Function to perform and plot GO enrichment ---
perform_and_plot_go <- function(gene_list, ontology_type, direction_label, file_suffix, results_obj_ensembl, universe_genes_ensembl) {
  
  # Ensure the gene list passed for enrichment is in Ensembl IDs
  # This step uses the 'results_obj_ensembl' to map symbols back to Ensembl if needed
  # or directly extracts Ensembl IDs if the initial filtering was done on Ensembl IDs
  
  # For your current setup where results_cpg_ifn (used for top50 heatmap) has symbols
  # and results_cpg_ifn_ensembl has Ensembl, we need to get Ensembl IDs of significant genes
  
  # Subset the results_obj_ensembl (which has Ensembl IDs as rownames)
  # based on the gene symbols found in the 'gene_list' (from results_cpg_ifn which had symbols)
  
  # First, ensure results_obj_ensembl has the 'symbol' column
  if (!"symbol" %in% colnames(results_obj_ensembl)) {
    results_obj_ensembl$symbol <- mapIds(org.Mm.eg.db, keys = rownames(results_obj_ensembl), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  }
  
  # Filter significant genes from the Ensembl-ID-based results object
  sig_genes_ensembl <- results_obj_ensembl[!is.na(results_obj_ensembl$padj) & 
                                             !is.na(results_obj_ensembl$symbol) &
                                             results_obj_ensembl$padj < padj_thresh & 
                                             abs(results_obj_ensembl$log2FoldChange) >= log2fc_thresh, ]
  
  if (direction_label == "up regulated") {
    gene_ids_for_enrichment <- rownames(sig_genes_ensembl[sig_genes_ensembl$log2FoldChange >= log2fc_thresh, ])
  } else if (direction_label == "down regulated") {
    gene_ids_for_enrichment <- rownames(sig_genes_ensembl[sig_genes_ensembl$log2FoldChange <= -log2fc_thresh, ])
  } else {
    stop("Invalid direction_label provided.")
  }
  
  if (length(gene_ids_for_enrichment) == 0) {
    message(paste("No", direction_label, "genes found for", ontology_type, ". Skipping enrichment."))
    return(NULL)
  }
  
  go_results <- enrichGO(
    gene = gene_ids_for_enrichment,
    universe = universe_genes_ensembl, # Use the global universe of Ensembl IDs
    OrgDb = org.Mm.eg.db,
    keyType = "ENSEMBL", # Keep as ENSEMBL as universe and gene IDs are Ensembl
    ont = ontology_type,
    pAdjustMethod = "BH", # Benjamini-Hochberg for FDR correction
    qvalueCutoff = padj_thresh,
    readable = TRUE # Converts Ensembl IDs to gene symbols in the result table
  )
  
  if (is.null(go_results) || nrow(go_results@result) == 0) {
    message(paste("No enriched GO terms found for", direction_label, "genes in", ontology_type, "."))
    return(NULL)
  }
  
  # Save result as csv file
  file_name <- paste0("GO-term_Cpg-IFN_24h_", direction_label, "-", ontology_type, file_suffix, ".csv")
  write.csv(go_results@result, file = file_name, row.names = F)
  message(paste("Saved GO results to", file_name))
  
  go_results@result$Description <- str_wrap(go_results@result$Description, width = 40)
  
  # Plot result
  plot_title <- paste("GO Enrichment of", direction_label, "genes -", ontology_type)
  plot_file_name <- paste0("GO-dotplot_CpG-IFN_24h_", direction_label, "-", ontology_type, file_suffix, ".png")
  
  png(plot_file_name, width = 800, height = 600)
  p <- dotplot(go_results, showCategory = 20, font.size = 10, label_format = 70) +
    scale_size_continuous(range = c(1, 7)) +
    theme_minimal() +
    ggtitle(plot_title)
  print(p) # Explicitly print plot when inside a function
  dev.off()
  message(paste("Saved GO dotplot to", plot_file_name))
  
  return(go_results)
}

# --- Run GO analysis for positive log2FoldChange (up-regulated genes) ---
message("--- Analyzing up-regulated genes ---")
# Biological Processes (BP)
GO_sig_cpg_ifn_pos_bp <- perform_and_plot_go(
  gene_list = NULL, # Placeholder, as gene filtering is done inside the function
  ontology_type = "BP",
  direction_label = "up regulated",
  file_suffix = "_BP",
  results_obj_ensembl = results_cpg_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Molecular Function (MF)
GO_sig_cpg_ifn_pos_mf <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "MF",
  direction_label = "up regulated",
  file_suffix = "_MF",
  results_obj_ensembl = results_cpg_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Cellular Component (CC)
GO_sig_cpg_ifn_pos_cc <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "CC",
  direction_label = "up regulated",
  file_suffix = "_CC",
  results_obj_ensembl = results_cpg_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)


# --- Run GO analysis for negative log2FoldChange (down-regulated genes) ---
message("--- Analyzing down-regulated genes ---")
# Biological Processes (BP)
GO_sig_cpg_ifn_neg_bp <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "BP",
  direction_label = "down regulated",
  file_suffix = "_BP",
  results_obj_ensembl = results_cpg_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Molecular Function (MF)
GO_sig_cpg_ifn_neg_mf <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "MF",
  direction_label = "down regulated",
  file_suffix = "_MF",
  results_obj_ensembl = results_cpg_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

# Cellular Component (CC)
GO_sig_cpg_ifn_neg_cc <- perform_and_plot_go(
  gene_list = NULL,
  ontology_type = "CC",
  direction_label = "down regulated",
  file_suffix = "_CC",
  results_obj_ensembl = results_cpg_ifn_ensembl,
  universe_genes_ensembl = universe_genes_ensembl
)

###############################################################################
#################################################################################


###----------------- GSEA analysis----------------------------


################################### CpG #############################
# --- Step 1: Prepare a pre-ranked list of genes with unique IDs ---

# Get the full results object with no filters applied
results_cpg_full <- results(dds, contrast = c("group", "CpG_yes", "CpG_no"))

# Convert results to a data frame and add an Ensembl ID column
results_df <- as.data.frame(results_cpg_full)
results_df$ensembl_id <- rownames(results_df)

# Map Ensembl IDs to Entrez IDs
results_df$entrez <- mapIds(
  org.Mm.eg.db,
  keys = results_df$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove rows with NA values (no Entrez ID or other missing data)
results_df <- na.omit(results_df)

# Handle duplicate Entrez IDs by keeping the one with the largest absolute log2FoldChange.
# This ensures each Entrez ID is unique and represents the most significant change.
results_df <- results_df[order(abs(results_df$log2FoldChange), decreasing = TRUE), ]
results_df <- results_df[!duplicated(results_df$entrez), ]

# Create the final ranked list with unique Entrez IDs
gene_list_entrez_ranked <- results_df$log2FoldChange
names(gene_list_entrez_ranked) <- results_df$entrez

# Sort the gene list in descending order for GSEA
gene_list_entrez_ranked <- sort(gene_list_entrez_ranked, decreasing = TRUE)

# --- Step 2: Load Gene Sets for GSEA ---

# Get the mouse GO gene sets from the C5 category
go_terms_df <- msigdbr(species = "Mus musculus", category = "C5")
term2gene_go <- go_terms_df %>% dplyr::select(gs_name, entrez_gene)

# --- Step 3: Run the GSEA Analysis ---
gsea_results <- GSEA(
  geneList = gene_list_entrez_ranked,
  TERM2GENE = term2gene_go,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# --- Step 4: Save and Visualize Results with ggplot2 ---
# Save the full results table
write.csv(gsea_results@result, file = "GSEA_GO_results_CpG_24h.csv", row.names = FALSE)
print("Full GO GSEA results table saved to 'GSEA_GO_results_CpG_24h.csv'")

if (nrow(gsea_results@result) > 0) {
  
  # Get data frame and wrap long descriptions
  plot_data <- gsea_results@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  # Dotplot of the top enriched pathways with ggplot2
  png("GSEA_GO_dotplot_CpG.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of GO Terms",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p1)
  dev.off()
  
  # Detailed enrichment plot for the top pathway
  png("GSEA_GO_enrichment_plot_CpG.png", width = 800, height = 600)
  p2 <- gseaplot2(gsea_results, geneSetID = 1, title = gsea_results$Description[1])
  print(p2)
  dev.off()
  
  print("GO GSEA analysis complete and plots saved.")
  
} else {
  print("No significantly enriched pathways found.")
}


####### Another GSEA category -------- 
# Get the mouse Hallmark gene sets
# --- Step 1 & 2: Prepare Gene List and Load Gene Sets ---
hallmark_df <- msigdbr(species = "Mus musculus", category = "H")

# Correct the gene set names by replacing underscores with spaces
hallmark_df$gs_name <- stringr::str_replace_all(hallmark_df$gs_name, "_", " ")
term2gene_hallmark <- hallmark_df %>% dplyr::select(gs_name, entrez_gene)

gsea_results_hallmark <- GSEA(
  geneList = gene_list_entrez_ranked,
  TERM2GENE = term2gene_hallmark,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# --- Step 4: Save and Visualize Results with ggplot2 ---
# Save the full results table
write.csv(gsea_results_hallmark@result, file = "GSEA_hallmark_results_CpG_24h.csv", row.names = FALSE)
print("Full GSEA results table saved to 'GSEA_hallmark_results_CpG_24h.csv'")

if (nrow(gsea_results_hallmark@result) > 0) {
  plot_data <- gsea_results_hallmark@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_hallmark_dotplot_CpG.png", width = 550, height = 300)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Hallmark Terms",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      # Increase the size of the plot title
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      # Increase the size of the x-axis title
      axis.title.x = element_text(size = 16, face = "bold"),
      # Increase the size of the y-axis title
      axis.title.y = element_text(size = 16, face = "bold"),
      # Increase the size of the axis tick labels
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      # Move legend to the right of the plot
      legend.position = "right"
    )
  
  print(p1)
  dev.off()
  
  png("GSEA_hallmark_enrichment_plot_CpG.png", width = 500, height = 400)
  p2 <- gseaplot(gsea_results_hallmark, geneSetID = 1, title = gsea_results_hallmark$Description[1])
  print(p2)
  dev.off()
  
  print("GSEA analysis complete and plots saved.")
  
} else {
  print("No significantly enriched pathways found.")
}


# Find the ID of the top enriched pathway from the results table
top_pathway_id <- gsea_results_hallmark@result$ID[1]

# Use the `gsInfo` function to get the leading edge genes for that ID
# The `gsInfo` function is part of the `enrichplot` package
# The "core_enrichment" column in your results table is the key
leading_edge_genes <- gsea_results_hallmark@result$core_enrichment[gsea_results_hallmark@result$ID == top_pathway_id]

# The genes are stored as a string, separated by slashes.
# We need to split this string to get a list of individual gene IDs.
gene_list <- unlist(strsplit(leading_edge_genes, "/"))

# Use mapIds to convert the Entrez IDs to gene symbols
# The `org.Mm.eg.db` package is for Mus musculus (mouse)
gene_symbols <- mapIds(
  x = org.Mm.eg.db,
  keys = gene_list,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)
# Print the list of gene IDs
print(gene_symbols)


# You can also get the number of leading edge genes
num_leading_edge <- length(gene_list)
print(paste("Number of leading edge genes:", num_leading_edge))



################ Immunologic gene set enrichment ###############
########################################################

# Get the mouse Immunologic gene sets from the C7 category

# --- Step 1 & 2: Prepare Gene List and Load Gene Sets ---
c7_df <- msigdbr(species = "Mus musculus", category = "C7")

# Correct the gene set names by replacing underscores with spaces
c7_df$gs_name <- stringr::str_replace_all(c7_df$gs_name, "_", " ")
term2gene_c7 <- c7_df %>% dplyr::select(gs_name, entrez_gene)

# Run the GSEA analysis with the new gene sets
gsea_results_c7 <- GSEA(
  geneList = gene_list_entrez_ranked,
  TERM2GENE = term2gene_c7,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# --- Step 4: Save and Visualize Results with ggplot2 ---
# Save the full results table
write.csv(gsea_results_c7@result, file = "GSEA_C7_results_CpG_24h.csv", row.names = FALSE)
print("Full GSEA C7 results table saved to 'GSEA_C7_results_CpG_24h.csv'")

if (nrow(gsea_results_c7@result) > 0) {
  plot_data <- gsea_results_c7@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_C7_dotplot_CpG.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Immunologic Signatures",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p1)
  dev.off()
  
  print("C7 GSEA analysis complete and plots saved.")
  
} else {
  print("No significantly enriched C7 pathways found.")
}


################ Canocical pathways ############
########################################################

# Get the mouse gene sets for canonical pathways (C2)
c2_df <- msigdbr(species = "Mus musculus", category = "C2")

# Filter for apoptotic pathways (you can add more keywords as needed)
apoptotic_pathways <- c2_df %>% 
  filter(grepl("APOPTOSIS|CELL_DEATH", gs_name, ignore.case = TRUE)) %>% 
  dplyr::select(gs_name, entrez_gene)

# Check how many pathways you have
print(paste("Found", length(unique(apoptotic_pathways$gs_name)), "apoptotic pathways."))


# --- Step 2: Prepare your ranked gene list ---
# Assuming 'res_cpg' is your full DESeq2 results object for the CpG comparison
# You need a ranked list of ALL genes, not just the significant ones
# You did this in your previous scripts, so this part should be familiar
res_cpg_shrink <- lfcShrink(dds, contrast = c("group", "CpG_yes", "CpG_no"), res = res_cpg, type = "ashr")

# Get a ranked list of all genes based on log2FoldChange
# The list should have Entrez IDs as names
res_cpg_shrink_df <- as.data.frame(res_cpg_shrink)
res_cpg_shrink_df$ensembl_id <- rownames(res_cpg_shrink_df)
res_cpg_shrink_df$entrez <- mapIds(
  org.Mm.eg.db,
  keys = res_cpg_shrink_df$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove rows with NA values
res_cpg_shrink_df <- na.omit(res_cpg_shrink_df)
res_cpg_shrink_df <- res_cpg_shrink_df[order(abs(res_cpg_shrink_df$log2FoldChange), decreasing = TRUE), ]
res_cpg_shrink_df <- res_cpg_shrink_df[!duplicated(res_cpg_shrink_df$entrez), ]

# Create the final ranked list
gene_list_entrez_ranked_cpg <- res_cpg_shrink_df$log2FoldChange
names(gene_list_entrez_ranked_cpg) <- res_cpg_shrink_df$entrez
gene_list_entrez_ranked_cpg <- sort(gene_list_entrez_ranked_cpg, decreasing = TRUE)

# --- Step 3: Run the GSEA Analysis for apoptotic pathways ---
gsea_results_apoptosis <- GSEA(
  geneList = gene_list_entrez_ranked_cpg,
  TERM2GENE = apoptotic_pathways,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# --- Step 4: Save and Visualize Results with ggplot2 ---
# Save the full results table
write.csv(gsea_results_apoptosis@result, file = "GSEA_apoptosis_results_CpG_24h.csv", row.names = FALSE)
print("Full GSEA Apoptosis results table saved to 'GSEA_apoptosis_results_CpG_24h.csv'")

if (nrow(gsea_results_apoptosis@result) > 0) {
  plot_data <- gsea_results_apoptosis@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_apoptosis_dotplot_CpG.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Apoptotic Pathways",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p1)
  dev.off()
  
  print("Apoptosis GSEA analysis complete and plots saved.")
  
} else {
  print("No significantly enriched apoptotic pathways found.")
}


######## transcript factor enrichment analysis
##################################################################

res_cpg_shrink <- lfcShrink(dds, contrast=c("group", "CpG_yes", "CpG_no"), res=res_cpg, type="ashr")

# Assuming 'res_cpg_shrink' is your shrunk DESeq2 result for CpG
significant_degs <- subset(res_cpg_shrink, padj < 0.05 & abs(log2FoldChange) > 1)

# Get the Ensembl IDs of your significant DEGs
degs_ensembl <- rownames(significant_degs)

# Map Ensembl IDs to Entrez IDs
degs_entrez <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = degs_ensembl,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove any genes that didn't map to an Entrez ID
degs_entrez <- degs_entrez[!is.na(degs_entrez)]

# Get only the Transcription Factor Targets (TFT) from the C3 collectio
# Get the full C3 collection (Regulatory Target Motifs)
c3_df <- msigdbr(species = "Mus musculus", category = "C3")

# Now, filter the data frame to get only the Transcription Factor Targets (TFT)
tft_df <- c3_df %>% filter(gs_subcat == "TFT")


# Run the enrichment analysis with the filtered TFT gene sets
tf_enrichment <- enricher(
  gene = degs_entrez,
  TERM2GENE = tft_df[, c("gs_name", "entrez_gene")],
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

# --- Step 3: Save and Visualize Results with ggplot2 ---
# Save the full results table
write.csv(tf_enrichment@result, file = "TFT_results_CpG_24h.csv", row.names = FALSE)
print("Full TFT enrichment results table saved to 'TFT_results_CpG_24h.csv'")

if (nrow(tf_enrichment@result) > 0) {
  plot_data <- tf_enrichment@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("TFT_dotplot_CpG.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = Count,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "Transcription Factor Enrichment",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p1)
  dev.off()
  
  print("TFT enrichment analysis complete and plots saved.")
  
} else {
  print("No significantly enriched TFT pathways found.")
}


#################################### IFN #############################
# --- Step 1: Prepare a pre-ranked list of genes with unique IDs ---

# Get the full results object with no filters applied
results_ifn_full <- results(dds, contrast = c("group", "IFN_yes", "IFN_no"))

# Convert results to a data frame and add an Ensembl ID column
results_df_ifn <- as.data.frame(results_ifn_full)
results_df_ifn$ensembl_id <- rownames(results_df_ifn)

# Map Ensembl IDs to Entrez IDs
results_df_ifn$entrez <- mapIds(
  org.Mm.eg.db,
  keys = results_df_ifn$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove rows with NA values (no Entrez ID or other missing data)
results_df_ifn <- na.omit(results_df_ifn)

# Handle duplicate Entrez IDs by keeping the one with the largest absolute log2FoldChange.
# This ensures each Entrez ID is unique and represents the most significant change.
results_df_ifn <- results_df_ifn[order(abs(results_df_ifn$log2FoldChange), decreasing = TRUE), ]
results_df_ifn <- results_df_ifn[!duplicated(results_df_ifn$entrez), ]

# Create the final ranked list with unique Entrez IDs
gene_list_entrez_ranked_ifn <- results_df_ifn$log2FoldChange
names(gene_list_entrez_ranked_ifn) <- results_df_ifn$entrez

# Sort the gene list in descending order for GSEA
gene_list_entrez_ranked_ifn <- sort(gene_list_entrez_ranked_ifn, decreasing = TRUE)

# --- GO Analysis (C5) ---
print("Starting GO GSEA Analysis...")
go_terms_df <- msigdbr(species = "Mus musculus", category = "C5")
term2gene_go <- go_terms_df %>% dplyr::select(gs_name, entrez_gene)

gsea_results_ifn <- GSEA(
  geneList = gene_list_entrez_ranked_ifn,
  TERM2GENE = term2gene_go,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(gsea_results_ifn@result, file = "GSEA_GO_results_IFN_24h.csv", row.names = FALSE)
print("Full GO GSEA results table saved to 'GSEA_GO_results_IFN_24h.csv'")

if (nrow(gsea_results_ifn@result) > 0) {
  plot_data <- gsea_results_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_GO_dotplot_IFN.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of GO Terms",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  png("GSEA_GO_enrichment_plot_ifn.png", width = 800, height = 600)
  p2 <- gseaplot2(gsea_results_ifn, geneSetID = 1, title = gsea_results_ifn$Description[1])
  print(p2)
  dev.off()
  
  print("GO GSEA analysis complete and plots saved.")
} else {
  print("No significantly enriched GO pathways found.")
}

---
  
  ### Hallmark GSEA Analysis (H)
  
  print("Starting Hallmark GSEA Analysis...")
hallmark_df <- msigdbr(species = "Mus musculus", category = "H")
hallmark_df$gs_name <- stringr::str_replace_all(hallmark_df$gs_name, "_", " ")
term2gene_hallmark <- hallmark_df %>% dplyr::select(gs_name, entrez_gene)

gsea_results_hallmark_ifn <- GSEA(
  geneList = gene_list_entrez_ranked_ifn,
  TERM2GENE = term2gene_hallmark,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(gsea_results_hallmark_ifn@result, file = "GSEA_hallmark_results_IFN_24h.csv", row.names = FALSE)
print("Full GSEA Hallmark results table saved to 'GSEA_hallmark_results_IFN_24h.csv'")

if (nrow(gsea_results_hallmark_ifn@result) > 0) {
  plot_data <- gsea_results_hallmark_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_hallmark_dotplot_IFN.png", width = 550, height = 300)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Hallmark Terms",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  print(p1)
  dev.off()
  
  png("GSEA_hallmark_enrichment_plot_IFN.png", width = 500, height = 400)
  p2 <- gseaplot(gsea_results_hallmark_ifn, geneSetID = 1, title = gsea_results_hallmark_ifn$Description[1])
  print(p2)
  dev.off()
  
  print("Hallmark GSEA analysis complete and plots saved.")
} else {
  print("No significantly enriched Hallmark pathways found.")
}

---
  
  ### Immunologic GSEA Analysis (C7)
  
  print("Starting Immunologic GSEA Analysis...")
c7_df <- msigdbr(species = "Mus musculus", category = "C7")
c7_df$gs_name <- stringr::str_replace_all(c7_df$gs_name, "_", " ")
term2gene_c7 <- c7_df %>% dplyr::select(gs_name, entrez_gene)

gsea_results_c7_ifn <- GSEA(
  geneList = gene_list_entrez_ranked_ifn,
  TERM2GENE = term2gene_c7,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(gsea_results_c7_ifn@result, file = "GSEA_C7_results_IFN_24h.csv", row.names = FALSE)
print("Full GSEA C7 results table saved to 'GSEA_C7_results_IFN_24h.csv'")

if (nrow(gsea_results_c7_ifn@result) > 0) {
  plot_data <- gsea_results_c7_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_C7_dotplot_IFN_24h.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Immunologic Signatures",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  print("C7 GSEA analysis complete and plots saved.")
} else {
  print("No significantly enriched C7 pathways found.")
}

---
  
  ### Apoptotic Canonical Pathway GSEA Analysis (C2)
  
  print("Starting Apoptosis GSEA Analysis...")
c2_df <- msigdbr(species = "Mus musculus", category = "C2")
apoptotic_pathways_ifn <- c2_df %>%
  filter(grepl("APOPTOSIS|CELL_DEATH", gs_name, ignore.case = TRUE)) %>%
  dplyr::select(gs_name, entrez_gene)

print(paste("Found", length(unique(apoptotic_pathways_ifn$gs_name)), "apoptotic pathways."))

# Note: This is an extra step to get the full ranked list from shrunk data
# It's not strictly necessary if you already have the ranked list above,
# but it's good practice for reproducibility if the list was prepared differently.
res_ifn_shrink <- lfcShrink(dds, contrast = c("group", "IFN_yes", "IFN_no"), res = res_IFN, type = "ashr")

res_ifn_shrink_df <- as.data.frame(res_ifn_shrink)
res_ifn_shrink_df$ensembl_id <- rownames(res_ifn_shrink_df)
res_ifn_shrink_df$entrez <- mapIds(
  org.Mm.eg.db,
  keys = res_ifn_shrink_df$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res_ifn_shrink_df <- na.omit(res_ifn_shrink_df)
res_ifn_shrink_df <- res_ifn_shrink_df[order(abs(res_ifn_shrink_df$log2FoldChange), decreasing = TRUE), ]
res_ifn_shrink_df <- res_ifn_shrink_df[!duplicated(res_ifn_shrink_df$entrez), ]

gene_list_entrez_ranked_ifn <- res_ifn_shrink_df$log2FoldChange
names(gene_list_entrez_ranked_ifn) <- res_ifn_shrink_df$entrez
gene_list_entrez_ranked_ifn <- sort(gene_list_entrez_ranked_ifn, decreasing = TRUE)

gsea_results_apoptosis_ifn <- GSEA(
  geneList = gene_list_entrez_ranked_ifn,
  TERM2GENE = apoptotic_pathways_ifn,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(gsea_results_apoptosis_ifn@result, file = "GSEA_apoptosis_results_IFN_24h.csv", row.names = FALSE)
print("Full GSEA Apoptosis results table saved to 'GSEA_apoptosis_results_IFN_24h.csv'")

if (nrow(gsea_results_apoptosis_ifn@result) > 0) {
  plot_data <- gsea_results_apoptosis_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_apoptosis_dotplot_IFN.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Apoptotic Pathways",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  print("Apoptosis GSEA analysis complete and plots saved.")
} else {
  print("No significantly enriched apoptotic pathways found.")
}

---
  
  ### Transcription Factor Enrichment Analysis (C3 TFT)
  
  print("Starting TFT Enrichment Analysis...")
res_ifn_shrink <- lfcShrink(dds, contrast=c("group", "IFN_yes", "IFN_no"), res=res_IFN, type="ashr")
significant_degs_ifn <- subset(res_ifn_shrink, padj < 0.05 & abs(log2FoldChange) > 1)

degs_ensembl_ifn <- rownames(significant_degs_ifn)
degs_entrez_ifn <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = degs_ensembl_ifn,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
degs_entrez_ifn <- degs_entrez_ifn[!is.na(degs_entrez_ifn)]

c3_df <- msigdbr(species = "Mus musculus", category = "C3")
tft_df <- c3_df %>% filter(gs_subcat == "TFT")

tf_enrichment_ifn <- enricher(
  gene = degs_entrez_ifn,
  TERM2GENE = tft_df[, c("gs_name", "entrez_gene")],
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(tf_enrichment_ifn@result, file = "TFT_results_IFN_24h.csv", row.names = FALSE)
print("Full TFT enrichment results table saved to 'TFT_results_IFN_24h.csv'")

if (nrow(tf_enrichment_ifn@result) > 0) {
  plot_data <- tf_enrichment_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("TFT_dotplot_IFN.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = Count,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "Transcription Factor Enrichment",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  print("TFT enrichment analysis complete and plots saved.")
} else {
  print("No significantly enriched TFT pathways found.")
}


#################################### CpG+IFN #############################
# --- Step 1: Prepare a pre-ranked list of genes with unique IDs ---

# Get the full results object with no filters applied
results_cpg_ifn_full <- results(dds, contrast = c("group", "CpG_IFN_yes", "CpG_IFN_no"))

# Convert results to a data frame and add an Ensembl ID column
results_df_cpg_ifn <- as.data.frame(results_cpg_ifn_full)
results_df_cpg_ifn$ensembl_id <- rownames(results_df_cpg_ifn)

# Map Ensembl IDs to Entrez IDs
results_df_cpg_ifn$entrez <- mapIds(
  org.Mm.eg.db,
  keys = results_df_cpg_ifn$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove rows with NA values (no Entrez ID or other missing data)
results_df_cpg_ifn <- na.omit(results_df_cpg_ifn)

# Handle duplicate Entrez IDs by keeping the one with the largest absolute log2FoldChange.
# This ensures each Entrez ID is unique and represents the most significant change.
results_df_cpg_ifn <- results_df_cpg_ifn[order(abs(results_df_cpg_ifn$log2FoldChange), decreasing = TRUE), ]
results_df_cpg_ifn <- results_df_cpg_ifn[!duplicated(results_df_cpg_ifn$entrez), ]

# Create the final ranked list with unique Entrez IDs
gene_list_entrez_ranked_cpg_ifn <- results_df_cpg_ifn$log2FoldChange
names(gene_list_entrez_ranked_cpg_ifn) <- results_df_cpg_ifn$entrez

# Sort the gene list in descending order for GSEA
gene_list_entrez_ranked_cpg_ifn <- sort(gene_list_entrez_ranked_cpg_ifn, decreasing = TRUE)

# --- GO Analysis (C5) ---
print("Starting GO GSEA Analysis...")
go_terms_df <- msigdbr(species = "Mus musculus", category = "C5")
term2gene_go <- go_terms_df %>% dplyr::select(gs_name, entrez_gene)

gsea_results_cpg_ifn <- GSEA(
  geneList = gene_list_entrez_ranked_cpg_ifn,
  TERM2GENE = term2gene_go,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(gsea_results_cpg_ifn@result, file = "GSEA_GO_results_CpG_IFN_24h.csv", row.names = FALSE)
print("Full GO GSEA results table saved to 'GSEA_GO_results_CpG_IFN_24h.csv'")

if (nrow(gsea_results_cpg_ifn@result) > 0) {
  plot_data <- gsea_results_cpg_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_GO_dotplot_CpG_IFN.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of GO Terms",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  png("GSEA_GO_enrichment_plot_cpg_ifn.png", width = 800, height = 600)
  p2 <- gseaplot2(gsea_results_cpg_ifn, geneSetID = 1, title = gsea_results_cpg_ifn$Description[1])
  print(p2)
  dev.off()
  
  print("GO GSEA analysis complete and plots saved.")
} else {
  print("No significantly enriched GO pathways found.")
}

---
  
  ### Hallmark GSEA Analysis (H)
  
  print("Starting Hallmark GSEA Analysis...")
hallmark_df <- msigdbr(species = "Mus musculus", category = "H")
hallmark_df$gs_name <- stringr::str_replace_all(hallmark_df$gs_name, "_", " ")
term2gene_hallmark <- hallmark_df %>% dplyr::select(gs_name, entrez_gene)

gsea_results_hallmark_cpg_ifn <- GSEA(
  geneList = gene_list_entrez_ranked_cpg_ifn,
  TERM2GENE = term2gene_hallmark,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(gsea_results_hallmark_cpg_ifn@result, file = "GSEA_hallmark_results_CpG_IFN_24h.csv", row.names = FALSE)
print("Full GSEA Hallmark results table saved to 'GSEA_hallmark_results_CpG_IFN_24h.csv'")

if (nrow(gsea_results_hallmark_cpg_ifn@result) > 0) {
  plot_data <- gsea_results_hallmark_cpg_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_hallmark_dotplot_CPG_IFN.png", width = 550, height = 300)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Hallmark Terms",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    )
  print(p1)
  dev.off()
  
  png("GSEA_hallmark_enrichment_plot_CpG_IFN.png", width = 500, height = 400)
  p2 <- gseaplot(gsea_results_hallmark_cpg_ifn, geneSetID = 1, title = gsea_results_hallmark_cpg_ifn$Description[1])
  print(p2)
  dev.off()
  
  print("Hallmark GSEA analysis complete and plots saved.")
} else {
  print("No significantly enriched Hallmark pathways found.")
}



---
  
  ### Immunologic GSEA Analysis (C7)
  
  print("Starting Immunologic GSEA Analysis...")
c7_df <- msigdbr(species = "Mus musculus", category = "C7")
c7_df$gs_name <- stringr::str_replace_all(c7_df$gs_name, "_", " ")
term2gene_c7 <- c7_df %>% dplyr::select(gs_name, entrez_gene)

gsea_results_c7_cpg_ifn <- GSEA(
  geneList = gene_list_entrez_ranked_cpg_ifn,
  TERM2GENE = term2gene_c7,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(gsea_results_c7_cpg_ifn@result, file = "GSEA_C7_results_CpG_IFN_24h.csv", row.names = FALSE)
print("Full GSEA C7 results table saved to 'GSEA_C7_results_CpG_IFN_24h.csv'")

if (nrow(gsea_results_c7_cpg_ifn@result) > 0) {
  plot_data <- gsea_results_c7_cpg_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_C7_dotplot_CpG_IFN_24h.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Immunologic Signatures",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  print("C7 GSEA analysis complete and plots saved.")
} else {
  print("No significantly enriched C7 pathways found.")
}

---
  
  ### Apoptotic Canonical Pathway GSEA Analysis (C2)
  
  print("Starting Apoptosis GSEA Analysis...")
c2_df <- msigdbr(species = "Mus musculus", category = "C2")
apoptotic_pathways_cpg_ifn <- c2_df %>%
  filter(grepl("APOPTOSIS|CELL_DEATH", gs_name, ignore.case = TRUE)) %>%
  dplyr::select(gs_name, entrez_gene)

print(paste("Found", length(unique(apoptotic_pathways_cpg_ifn$gs_name)), "apoptotic pathways."))

# Note: This is an extra step to get the full ranked list from shrunk data
# It's not strictly necessary if you already have the ranked list above,
# but it's good practice for reproducibility if the list was prepared differently.
res_cpg_ifn_shrink <- lfcShrink(dds, contrast = c("group", "CpG_IFN_yes", "CpG_IFN_no"), res = res_CpG_IFN, type = "ashr")

res_cpg_ifn_shrink_df <- as.data.frame(res_cpg_ifn_shrink)
res_cpg_ifn_shrink_df$ensembl_id <- rownames(res_cpg_ifn_shrink_df)
res_cpg_ifn_shrink_df$entrez <- mapIds(
  org.Mm.eg.db,
  keys = res_cpg_ifn_shrink_df$ensembl_id,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
res_cpg_ifn_shrink_df <- na.omit(res_cpg_ifn_shrink_df)
res_cpg_ifn_shrink_df <- res_cpg_ifn_shrink_df[order(abs(res_cpg_ifn_shrink_df$log2FoldChange), decreasing = TRUE), ]
res_cpg_ifn_shrink_df <- res_cpg_ifn_shrink_df[!duplicated(res_cpg_ifn_shrink_df$entrez), ]

gene_list_entrez_ranked_cpg_ifn <- res_cpg_ifn_shrink_df$log2FoldChange
names(gene_list_entrez_ranked_cpg_ifn) <- res_cpg_ifn_shrink_df$entrez
gene_list_entrez_ranked_cpg_ifn <- sort(gene_list_entrez_ranked_cpg_ifn, decreasing = TRUE)

gsea_results_apoptosis_cpg_ifn <- GSEA(
  geneList = gene_list_entrez_ranked_cpg_ifn,
  TERM2GENE = apoptotic_pathways_cpg_ifn,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(gsea_results_apoptosis_cpg_ifn@result, file = "GSEA_apoptosis_results_CpG_IFN_24h.csv", row.names = FALSE)
print("Full GSEA Apoptosis results table saved to 'GSEA_apoptosis_results_CpG_IFN_24h.csv'")

if (nrow(gsea_results_apoptosis_cpg_ifn@result) > 0) {
  plot_data <- gsea_results_apoptosis_cpg_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("GSEA_apoptosis_dotplot_CpG_IFN.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = setSize,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "GSEA of Apoptotic Pathways",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  print("Apoptosis GSEA analysis complete and plots saved.")
} else {
  print("No significantly enriched apoptotic pathways found.")
}

---
  
  ### Transcription Factor Enrichment Analysis (C3 TFT)
  
  print("Starting TFT Enrichment Analysis...")
res_cpg_ifn_shrink <- lfcShrink(dds, contrast=c("group", "CpG_IFN_yes", "CpG_IFN_no"), res=res_CpG_IFN, type="ashr")
significant_degs_cpg_ifn <- subset(res_cpg_ifn_shrink, padj < 0.05 & abs(log2FoldChange) > 1)

degs_ensembl_cpg_ifn <- rownames(significant_degs_cpg_ifn)
degs_entrez_cpg_ifn <- AnnotationDbi::mapIds(
  org.Mm.eg.db,
  keys = degs_ensembl_cpg_ifn,
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)
degs_entrez_cpg_ifn <- degs_entrez_cpg_ifn[!is.na(degs_entrez_cpg_ifn)]

c3_df <- msigdbr(species = "Mus musculus", category = "C3")
tft_df <- c3_df %>% filter(gs_subcat == "TFT")

tf_enrichment_cpg_ifn <- enricher(
  gene = degs_entrez_cpg_ifn,
  TERM2GENE = tft_df[, c("gs_name", "entrez_gene")],
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH"
)

write.csv(tf_enrichment_cpg_ifn@result, file = "TFT_results_CpG_IFN_24h.csv", row.names = FALSE)
print("Full TFT enrichment results table saved to 'TFT_results_CpG_IFN_24h.csv'")

if (nrow(tf_enrichment_cpg_ifn@result) > 0) {
  plot_data <- tf_enrichment_cpg_ifn@result
  plot_data$Description <- stringr::str_wrap(plot_data$Description, width = 40)
  
  png("TFT_dotplot_CpG_IFN.png", width = 600, height = 500)
  p1 <- ggplot(plot_data, aes(
    x = -log10(p.adjust),
    y = reorder(Description, -log10(p.adjust)),
    size = Count,
    color = p.adjust
  )) +
    geom_point() +
    scale_color_viridis_c(direction = -1) +
    labs(
      title = "Transcription Factor Enrichment",
      x = "-log10(Adjusted P-value)",
      y = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  print(p1)
  dev.off()
  
  print("TFT enrichment analysis complete and plots saved.")
} else {
  print("No significantly enriched TFT pathways found.")
}



#########################################################################################
#########################################################################################

############---------------  GSVA analysis -------------------------------------------------------


# --- Step 1: Prepare the data for limma with Entrez IDs ---

# Get the DESeq2 normalized counts
dge <- DGEList(counts = counts(dds))
dge <- calcNormFactors(dge)

# Define the design matrix based on your sample groups
limma_design <- model.matrix(~group, data = colData(dds))

# Use 'voom' to prepare the data for limma
v <- voom(dge, limma_design, plot = TRUE)

# --- NEW: Convert the row names (Ensembl) to Entrez IDs ---

# Get the Ensembl to Entrez ID mapping for all genes in your data
ensembl_to_entrez <- mapIds(
  org.Mm.eg.db,
  keys = rownames(v),
  column = "ENTREZID",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# Remove genes with no corresponding Entrez ID
v <- v[!is.na(ensembl_to_entrez), ]
ensembl_to_entrez <- ensembl_to_entrez[!is.na(ensembl_to_entrez)]

# Handle duplicates: multiple Ensembl IDs might map to the same Entrez ID.
# This code removes the duplicate Entrez IDs, keeping the first occurrence.
v <- v[!duplicated(ensembl_to_entrez), ]
ensembl_to_entrez <- ensembl_to_entrez[!duplicated(ensembl_to_entrez)]

# Set the row names of the 'voom' object to be the Entrez IDs
rownames(v) <- ensembl_to_entrez

# --- Step 2: Load Gene Sets with Matching Entrez IDs ---

# The hallmark_list from msigdbr already uses Entrez IDs, so this step is simple.
hallmark_msigdb <- msigdbr(species = "Mus musculus", category = "H")
hallmark_list <- split(hallmark_msigdb$entrez_gene, hallmark_msigdb$gs_name)

# --- Step 3: Run the Gene Set Analysis with the corrected gene list ---

# Define the contrasts for your comparisons
my_contrasts <- makeContrasts(
  CpG_vs_no_CpG = groupCpG_yes,
  IFN_vs_no_IFN = groupIFN_yes,
  CpG_IFN_vs_no_CpG_IFN = groupCpG_IFN_yes,
  levels = limma_design
)

######## For CpG ####################

# Run 'roast' analysis for the CpG comparison
roast_results_cpg <- roast(
  v,
  index = hallmark_list,
  design = limma_design,
  contrast = my_contrasts[, "CpG_vs_no_CpG"]
)

# Run 'camera' analysis for the CpG comparison
camera_results_cpg <- camera(
  v,
  index = hallmark_list,
  design = limma_design,
  contrast = my_contrasts[, "CpG_vs_no_CpG"]
)

# View the top results from both methods
print("Top pathways from roast (CpG vs. no CpG):")
print(head(roast_results_cpg, 20))
print("Top pathways from camera (CpG vs. no CpG):")
print(head(camera_results_cpg, 20))


# Assuming you want to plot the top 20 pathways from the camera results
camera_results_top20 <- head(camera_results_cpg, 20)

sum(camera_results_cpg$FDR < 0.05)



# Create a bar plot
camera_barplot_cpg <-ggplot(camera_results_top20, aes(x = reorder(factor(rownames(camera_results_top20)), -log10(FDR)), y = -log10(FDR), fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top Enriched Hallmark Pathways: CpG - Chemokine vs. no Chemokine",
    x = "Pathway",
    y = "-log10(FDR)"
  ) +
  theme_classic()

ggsave("camera_barplot_cpg.png", plot = camera_barplot_cpg, width = 10, height = 8, dpi = 300)

# Save the full roast results to a CSV file
write.csv(roast_results_cpg, "roast_results_cpg.csv")

# Save the full camera results to a CSV file
write.csv(camera_results_cpg, "camera_results_cpg.csv")


# Filter the full results to keep only significant pathways (FDR < 0.05)
camera_results_significant <- subset(camera_results_cpg, FDR < 0.05)

# If any significant pathways were found, plot the top 20 of them
if (nrow(camera_results_significant) > 0) {
  # Order by significance and take the top 20
  camera_results_top20 <- camera_results_significant[order(camera_results_significant$FDR), ][1:20, ]
  
  # Now, create the bar plot using this filtered data frame
  camera_barplot_cpg <- ggplot(camera_results_top20, aes(x = reorder(factor(rownames(camera_results_top20)), -log10(FDR)), y = -log10(FDR), fill = Direction)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = "Top 20 Significant Hallmark Pathways: CpG - Chemokine vs. no Chemokine",
      x = "Pathway",
      y = "-log10(FDR)"
    ) +
    theme_classic()
  
  print(camera_barplot_cpg)
  
} else {
  print("No significantly enriched pathways found (FDR < 0.05).")
}

######## For IFN ####################

# Run 'roast' analysis for the IFN comparison
roast_results_ifn <- roast(
  v,
  index = hallmark_list,
  design = limma_design,
  contrast = my_contrasts[, "IFN_vs_no_IFN"]
)

# Run 'camera' analysis for the IFN comparison
camera_results_ifn <- camera(
  v,
  index = hallmark_list,
  design = limma_design,
  contrast = my_contrasts[, "IFN_vs_no_IFN"]
)

# View the top results from both methods
print("Top pathways from roast (IFN vs. no IFN):")
print(head(roast_results_ifn, 20))
print("Top pathways from camera (IFN vs. no IFN):")
print(head(camera_results_ifn, 20))


# Assuming you want to plot the top 20 pathways from the camera results
camera_results_top20 <- head(camera_results_ifn, 9)

sum(camera_results_ifn$FDR < 0.05)

# Create a bar plot
camera_barplot_ifn <- ggplot(camera_results_top20, aes(x = reorder(factor(rownames(camera_results_top20)), -log10(FDR)), y = -log10(FDR), fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top Enriched Hallmark Pathways: IFN - Chemokine vs. no Chemokine",
    x = "Pathway",
    y = "-log10(FDR)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.8),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )

ggsave("camera_barplot_ifn.png", plot = camera_barplot_ifn, width = 10, height = 8, dpi = 300)

# Save the full roast results to a CSV file
write.csv(roast_results_ifn, "roast_results_ifn.csv")

# Save the full camera results to a CSV file
write.csv(camera_results_ifn, "camera_results_ifn.csv")




######## For CpG+IFN ####################
# Run 'roast' analysis for the CpG comparison
roast_results_cpg_ifn <- roast(
  v,
  index = hallmark_list,
  design = limma_design,
  contrast = my_contrasts[, "CpG_IFN_vs_no_CpG_IFN"]
)

# Run 'camera' analysis for the CpG comparison
camera_results_cpg_ifn <- camera(
  v,
  index = hallmark_list,
  design = limma_design,
  contrast = my_contrasts[, "CpG_IFN_vs_no_CpG_IFN"]
)

# View the top results from both methods
print("Top pathways from roast (CpG+IFN vs. no CpG+IFN):")
print(head(roast_results_cpg_ifn, 20))
print("Top pathways from camera (CpG+IFN vs. no CpG+IFN):")
print(head(camera_results_cpg_ifn, 20))


# Assuming you want to plot the top 20 pathways from the camera results
camera_results_top20 <- head(camera_results_cpg_ifn, 20)

sum(camera_results_cpg_ifn$FDR < 0.05)

camera_results_top20 <- head(camera_results_cpg_ifn, 8)

# Create a bar plot
camera_barplot_cpg_ifn <-ggplot(camera_results_top20, aes(x = reorder(factor(rownames(camera_results_top20)), -log10(FDR)), y = -log10(FDR), fill = Direction)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Top Enriched Hallmark Pathways: CpG+IFN - Chemokine vs. no Chemokine",
    x = "Pathway",
    y = "-log10(FDR)"
  ) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.8),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12)
  )

ggsave("camera_barplot_cpg_ifn.png", plot = camera_barplot_cpg_ifn, width = 10, height = 8, dpi = 300)

# Save the full roast results to a CSV file
write.csv(roast_results_cpg_ifn, "roast_results_cpg_ifn.csv")

# Save the full camera results to a CSV file
write.csv(camera_results_cpg_ifn, "camera_results_cpg_ifn.csv")

