### Get spatial regression coefficients for Active Score

library(Seurat)
library(cluster)
library(ggplot2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(EnhancedVolcano)
library(purrr)


source('/path_to_CALISTA/ModuleA.R')

#Set parameters
letters <- c("A", "B", "C", "D", "E", "G", "I", "K", "M")
size <- c(rep(10500, 4), 5500, 3000, 4000, 5500, 3500)
names(size) <- letters

#Load object 
object.interface <- readRDS('/path_to_data/interface_object.rds')

#subset to keep only TSLI region 
object.inter <- subset(object.interface, subset = interface_w_border == 'Border')

#read sender df 
receiver.df <- read.csv('/path_to_data/receiver_df.csv')
object.inter <- AddMetaData(object.inter, receiver.df)

#get data for every sample 
coords.list <- list()
x.list <- list()
y.list <-  list()
for(sample in letters){
  #Subset sample
  object.sample <- subset(object.inter, subset = orig.ident == paste0('sample_', sample))
  
  #read clone annotation 
  input_dir <- "/Users/lfq7/Library/CloudStorage/OneDrive-RutgersUniversity/Lab/Data/Bladder_Cancer/Spatial_Transcriptomics/Seurat_objects/dat.so.infercnv."
  clone.object <- readRDS(paste0(input_dir, sample, ".rds"))
  object.sample$clone <- as.factor(clone.object@meta.data$infercnv_subcluster)
  
  #Get coords list
  coords.list[[sample]] <- GetTissueCoordinates(object.sample)
  
  #Get minor cell types
  minor.cell.types <- as.factor(object.sample$cell_types)
  minor.cell.types.mat <- model.matrix(~ minor.cell.types  - 1)
  minor.cell.types.mat <- as.data.frame(minor.cell.types.mat )
  colnames(minor.cell.types.mat) <- gsub('minor.cell.types', '', colnames(minor.cell.types.mat))
  minor.cell.types.mat <- minor.cell.types.mat[,c('Bcell', 'Tcell', 'Endothelial', 'Icaf', 'Macrophage', 'Mcaf', 'Myeloid', 'Tcell')]
  
  #Get clones
  clones <- as.factor(object.sample$clone)
  levels(clones) <- c(levels(clones), "0")
  clones[is.na(clones)] <- "0"
  clones.mat <- model.matrix(~ clones  - 1)
  clones.mat <- as.data.frame(clones.mat[,c("clones1", "clones2", "clones3", "clones4", "clones5")])
  
  #Get LR MAT
  LR.mat <- object.sample@meta.data[,grepl("^r\\.", colnames(object.sample@meta.data))]
  LR.mat <- as.data.frame(LR.mat)
  
  #Combine all into one mat
  X.mat <- cbind(minor.cell.types.mat, clones.mat, LR.mat)
  X.mat <- X.mat[, colSums(is.na(X.mat)) != nrow(X.mat)]
  x.list[[sample]] <- X.mat
  
  
  #Get cancer Hallmarks mat and get y.list
  active.score <- as.data.frame(object.sample$Active_Score1)
  colnames(active.score)[1] <- 'Active_Score'
  y.list[[sample]] <-  as.data.frame(active.score)
  
}

#Get Hallmark Coefficient
pathways.coeff <- pathways.spatial.regression(coords.list, y.list, x.list, z.list = NULL, to.analyze = NULL, cells.to.ignore = NULL)

#Analyze each sample
sample <- 'E'
pathways.coeff.mat <- pathways.coeff[[sample]]$Active_Score
df <- as.data.frame(pathways.coeff.mat)

#Make Volcano Plot 
EnhancedVolcano(
  df,
  lab = rownames(df),
  x = 'coefficients',
  y = 'p.value',
  pCutoff = 0.05,            
  FCcutoff = 0.02,           
  pointSize = 3.0,
  labSize = 4.0,
  title = 'Regression Coefficients for Activity Score',
  subtitle = paste0('Sample ', sample),
  xlab = 'Coefficient',
  ylab = bquote(~-Log[10] ~ italic(P)),
  legendPosition = 'right',
  xlim = c(-0.4, 0.6),        # Limits the X-axis range
  ylim = c(0, 15),
  selectLab = rownames(df)[df$p.value < 0.05 | abs(df$coefficients) > 0.02]
)


#extract coefficients and save
df_combined <- map_dfr(pathways.coeff, ~ .x$Active_Score, .id = "Sample") %>%
  mutate(Feature = unlist(lapply(pathways.coeff, function(x) rownames(x$Active_Score))))
write.csv(df_combined, '/path_to_data/regression_coeffs.csv')
