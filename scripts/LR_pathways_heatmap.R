#### Make heatmap of sending pathways

library(ComplexHeatmap)
library(stringr)
library(circlize)
library(tidyverse)
library(dplyr)
library(rstatix)
library(stats)
library(tidyr)

#read sender df 
sender.df <- read.csv('/path_to_data/sender_df.csv')

#read object with annotation 
tri_state_object <- readRDS('/path_to_data/interface_object.rds')
tri_state_object.subset <- subset(tri_state_object, subset = interface_w_border == 'Border')
length(tri_state_object$global_two_state)

#Add samples and annotation to sender df
sender.df$Samples <- str_sub(sender.df$X, -1, -1)
sender.df$annotation = tri_state_object.subset$global_two_state
sender.df <- as_tibble(sender.df) %>% group_by(Samples, annotation) %>% arrange(Samples, annotation) 

#get df to plot 
top.plot <- sender.df[,-c(1, 38, 104)]
top.plot[is.na(top.plot)] <- 0
colnames(top.plot) <- gsub('s.', '', colnames(top.plot))

#Add annotation
col.anno <- HeatmapAnnotation(df = data.frame(Samples = sender.df[,c('Samples')], annotation = sender.df$annotation), 
                              col =list(annotation = c("Active" = '#FF3C31',"Passive" = "#5B8AFD", 'Intermediate' = '#b5179e'),
                                        Samples = c("A" = '#F9766D', "B" ='#C19746', "C" ='#9AA848', "D" ='#6EB55A',
                                                    "E" ='#72BCA2', "G" ='#6DB3DA', "I" ='#789AF1', "K" ='#C17CED',"M" ='#ED6CC0')))

#Plot only selected genes 
mat <- as.matrix(t(top.plot))
selected_genes <- c("WNT", "MK", "ncWNT", "SPP1", "ANGPTL", "MIF", "SEMA3", "TNF", "PTN", "FGF", 'PDGF', 'TGFB', 'VEGF', 'VISFATIN') 
selected_indices <- which(rownames(mat) %in% selected_genes)

#Set annotation
right_anno <- rowAnnotation(
  mark = anno_mark(
    at = selected_indices, 
    labels = rownames(mat)[selected_indices],
    labels_gp = gpar(fontsize = 8)
  )
)

#make heatmap 
col_fun = colorRamp2(c(0,1), c("lightgrey", "red"))
Heatmap(mat, row_names_gp = gpar(fontsize = 5), 
        cluster_rows = T, name = 'Score', top_annotation = col.anno, 
        col = col_fun, cluster_columns = F,
        show_row_names = FALSE, 
        right_annotation = right_anno,)

Heatmap(mat, row_names_gp = gpar(fontsize = 6), 
        cluster_rows = T, name = 'Score', top_annotation = col.anno, 
        col = col_fun, cluster_columns = F,
        show_row_names = T, column_title = 'LR Signaling at TSLIs',
        column_title_gp = gpar(fontsize = 14, fontface = "bold"))




#Make hallmark heatmap
hallmark.mat <- tri_state_object.subset@meta.data[, grepl('HALLMARK_', colnames(tri_state_object.subset@meta.data))]
hallmark.mat <- as.data.frame(hallmark.mat)
hallmark.mat$annotation = tri_state_object.subset$global_two_state
hallmark.mat$Samples <- gsub('sample_', '', tri_state_object.subset$orig.ident)
hallmark.mat<- as_tibble(hallmark.mat) %>% group_by(Samples, annotation) %>% arrange(Samples, annotation) 


colnames(hallmark.mat) <- gsub("HALLMARK_", "", colnames(hallmark.mat))
col_fun = colorRamp2(c(-1,0,1), c('blue',"lightgrey", "red"))
Heatmap(as.matrix(t(hallmark.mat[,1:50])), row_names_gp = gpar(fontsize = 9), col = col_fun,
        cluster_rows = T, name = 'Score', top_annotation = col.anno, cluster_columns = F, column_names_gp = gpar(fontsize = 0), 
        column_title = 'Hallmark Pathway Scores', column_title_gp = gpar(fontface = "bold"))


#Add annotation 
col.anno <- HeatmapAnnotation(df = sender.df[,c('Samples', 'Activity')], 
                              col =list(Activity = c("Active" = '#FF3C31', 'Intermediate' = 'pink', "Passive" = "#5B8AFD"),
                                        Samples = c("A" = '#F9766D', "B" ='#C19746', "C" ='#9AA848', "D" ='#6EB55A',
                                                    "E" ='#72BCA2', "G" ='#6DB3DA', "I" ='#789AF1', "K" ='#C17CED',"M" ='#ED6CC0')))

#make heatmap 
col_fun = colorRamp2(c(0,1), c("lightgrey", "red"))
Heatmap(as.matrix(t(top.plot)), row_names_gp = gpar(fontsize = 5), 
        cluster_rows = T, name = 'Score', top_annotation = col.anno, 
        col = col_fun, cluster_columns = F)



#### 3. Perform Krustal Wallis test on LR pathways ##############################

#Format mat
analysis_df <- as.data.frame(top.plot)
analysis_df$annotation <- factor(sender.df$annotation)
feature_cols <- colnames(top.plot)

#Run Kruskal-Wallis test 
kw_results <- lapply(feature_cols, function(gene) {
  # Skip zero-variance features
  if (length(unique(analysis_df[[gene]])) > 1) {
    res <- kruskal.test(formula(paste0("`", gene, "` ~ annotation")), data = analysis_df)
    data.frame(
      feature = gene,
      statistic = res$statistic,
      pvalue = res$p.value,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      feature = gene,
      statistic = 0,
      pvalue = 1,
      stringsAsFactors = FALSE
    )
  }
}) %>% bind_rows()

#Do FDR adjustment
kw_results <- kw_results %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  filter(padj < 0.05) %>%
  arrange(padj)

#plot results of significant LR 
Heatmap(as.matrix(t(top.plot[kw_results$feature])), row_names_gp = gpar(fontsize = 5), 
        cluster_rows = T, name = 'Score', top_annotation = col.anno, 
        col = col_fun, cluster_columns = F)




### 3.2 Plot only pathways that are expressed in more than one sample
plot_matrix <- top.plot
plot_matrix$Samples <- sender.df$Samples

#Count the number of distinct samples with non-zero expression per pathway
pathway_sample_counts <- sapply(colnames(top.plot), function(gene) {
  plot_matrix %>%
    filter(.data[[gene]] > 0) %>%
    pull(Samples) %>%
    n_distinct()
})

#Filter to keep only pathways expressed in more than one sample
multi_sample_pathways <- names(pathway_sample_counts[pathway_sample_counts > 1])
top.plot_filtered <- top.plot[, multi_sample_pathways]

#Make heatmap
Heatmap(
  as.matrix(t(top.plot_filtered[kw_results$feature[kw_results$feature %in% colnames(top.plot_filtered)]])), 
  row_names_gp = gpar(fontsize = 7), 
  cluster_rows = TRUE, 
  name = 'Score', 
  top_annotation = col.anno, 
  col = col_fun, 
  cluster_columns = FALSE,
)



### 3.3 Plot only pathways that are higher at a particular cluster
sig_genes <- kw_results$feature

#Calculate mean score per annotation group for each LR
group_means <- analysis_df %>%
  group_by(annotation) %>%
  summarise(across(all_of(sig_genes), mean, na.rm = TRUE)) %>%
  pivot_longer(
    cols = -annotation,
    names_to = "feature",
    values_to = "mean_score"
  )

#Annotate highest group
kw_results_enriched <- group_means %>%
  group_by(feature) %>%
  slice_max(order_by = mean_score, n = 1) %>%
  rename(highest_group = annotation, max_mean_score = mean_score) %>%
  inner_join(kw_results, by = "feature") %>%
  arrange(highest_group, padj)

#Save res
write.csv(kw_results_enriched, '/path_to_data/LR_signaling_state_comparison_results.csv')

#plot one cluster
kw_results_enriched.filtered <- kw_results_enriched %>% filter(highest_group == 'Active')
Heatmap(
  as.matrix(t(top.plot_filtered[kw_results_enriched.filtered $feature[kw_results_enriched.filtered$feature %in% colnames(top.plot_filtered)]])), 
  row_names_gp = gpar(fontsize = 7), 
  cluster_rows = TRUE, 
  name = 'Score', 
  top_annotation = col.anno, 
  col = col_fun, 
  cluster_columns = FALSE,
  column_title = "LR Upregulted in 'Passive' State",
  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)




kw_results_enriched.filtered <- kw_results_enriched %>% filter(highest_group == 'Intermediate')

sender.df.filtered <- sender.df[sender.df$annotation == 'Intermediate',]
kw_results_enriched.filtered <- kw_results_enriched.filtered[,sender.df$annotation == 'Intermediate']


#### 4. Make violin plot of selected pathways 

#get data for specific 
sample_target <- 'K'
plotting.df <- read.csv('/path_to_data/receiver_df.csv')
plotting.df$Samples <- str_sub(sender.df$X, -1, -1)
plotting.df$annotation = tri_state_object.subset$global_tri_state

#filter
plotting.df_filtered <- plotting.df %>% 
  filter(Sample == sample_target)

#Select LR
target_pathway <- "r.IL1"  # Replace with any pathway column name from sender.df

#Set colors
annotation_colors <- c(
  "Active" = "#FF3C31",
  "Intermediate" = "#b5179e",
  "Passive" = "#5B8AFD"
)

#Plot violin plot
ggplot(plotting.df_filtered , aes(x = annotation, y = .data[[target_pathway]], fill = annotation)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.85) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha = 0.6) +
  scale_fill_manual(values = annotation_colors) +
  theme_minimal(base_size = 14) +
  labs(
    title = paste0("LR Signaling Score: ", gsub('r.', '', target_pathway)),
    x = "Activity State",
    y = "Signaling Score"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank()
  )
ggsave(paste0('/path_to_data/global_3_state/', sample_target, 
              '/', gsub('r.', '', target_pathway), '_violinplot.png'), units = 'in', height = 4, width = 4)


#### 5. Barplot of total sending per region ######################################

#read sender df 
sender.df <- read.csv('/path_to_data/sender_df.csv')

#read new annotation 
tri_state_object <- readRDS('/Users/lfq7/Documents/Lab/Projects/Border_transition/global_states/interface_object.rds')
tri_state_object.subset <- subset(tri_state_object, subset = interface_w_border == 'Border')
length(tri_state_object$global_two_state)

#add annotation to df
sender.df$Samples <- str_sub(sender.df$X, -1, -1)
sender.df$annotation = tri_state_object.subset$global_two_state
sender.df <- as_tibble(sender.df) %>% group_by(Samples, annotation) %>% arrange(Samples, annotation) 

#sum by sample
df_summed <- sender.df %>%
  group_by(Samples) %>%
  summarise(spot_sum = sum(c_across(where(is.numeric)), na.rm = TRUE), .groups = "drop")

#plot LR singlaing by sample
ggplot(df_summed, aes(x = Samples, y = spot_sum, fill = Samples)) +
  geom_col() +
  labs(
    title = "LR Signaling by Sample",
    x = "Sample",
    y = "Total Signaling"
  ) +
  theme_minimal()

#Plot mean Signaling at Activity states
df_summed <- sender.df %>%
  group_by(Samples, annotation) %>%
  summarise(spot_sum = mean(c_across(where(is.numeric)), na.rm = TRUE), .groups = "drop")
ggplot(df_summed, aes(x = Samples, y = spot_sum, fill = annotation)) +
  scale_fill_manual(values = c("Active" = "red", "Passive" = "#6588F5")) +
  geom_col(position = "dodge") +
  labs(
    title = "Mean Signaling of Spots in Active State",
    x = "Sample",
    y = "Mean Signaling Score",
    fill = "Active State"
  ) +
  theme_minimal()
