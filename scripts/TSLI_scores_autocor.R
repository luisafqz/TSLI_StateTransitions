#### Get autocorrelation of Activity and significant pathways

library(Seurat)
source('/path_to_data/ModuleA.R')


#set parameters
base.dir <- '/path_to_data/'
letters <- c("A", "B", "C", "D", "E", "G", "I", "K", "M")
sample.num <- 9
samples <- c("sampleA", "sampleB", "sampleC", "sampleD", "sampleE", "sampleG", "sampleI", 
             "sampleK", "sampleM")
size <- c(rep(10500, 4), 5500, 3000, 4000, 5500, 3500)
cols <- c('Tumor' = '#F5DE95', 'Stroma' = "#DABCDF", 'Interface' = 'grey', 
          "Active" = '#FF3C31', "Passive" = "#5B8AFD", 'Intermediate' = '#b5179e')


#read object
object.interface <- readRDS('/path_to_data/interface_object.rds')

#subset to keep only TSLI region 
object.tsli <- subset(object.interface, subset = interface_w_border == 'Border')

#readlists 
MDA.list <- readRDS('/path_to_data/global_3_state/MDA_list.rds')
MDI.list <- readRDS('/path_to_data/global_3_state/MDI_list.rds')

# Get matching MDA and MDI
match.sig <- intersect(MDI.list, MDA.list)
match.sig <- gsub('EMT', 'EPITHELIAL_MESENCHYMAL_TRANSITION', match.sig)
match.sig <- paste0('HALLMARK_', match.sig)

#subset by sample 
coords.list <- list()
scores.list <- list()
for(i in letters){
  object.sample <- subset(object.tsli, subset = orig.ident == paste0('sample_', i))
  coords.list[[i]] <- GetTissueCoordinates(object.sample)
  clusters <- as.factor((object.sample$global_tri_state))
  X.mat <- model.matrix(~ clusters  - 1)
  scores.list[[i]] <- cbind(object.sample@meta.data[c(match.sig, 'Active_Score1')], X.mat)
}

#Get autocor of pathways
autocor.res <- scores.autocor(coords.list, scores.list)
autocor.res.pvalues <- scores.autocor(coords.list, scores.list, return.pvalue = T)

#save results 
write.csv(autocor.res, '/path_to_data/supplementary_tables/tsli_autocor_res.csv')
write.csv(autocor.res.pvalues, '/path_to_data/supplementary_tables/tsli_autocor_pvals_res.csv')

#plot autocor
Heatmap(autocor.res, name = "Moran's I", column_names_gp = gpar(fontsize = 15),
        column_title = 'Autocorrelation in TSLIs', column_title_gp = gpar(fontsize = 15, fontface = 'bold'))

