library(SPOTlight)
library(fda)
library(magrittr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
setwd("D:\\study\\undergraduate\\spolight\\code\\fda_ST\\fun_spotlight")
source("fun_spotlight.r")

#parameter
n_syn_mixt = 1000
data_name = "cerebellum_singlecell"
marker_gene_access = 1
hvg = 0
num_cl_n = 20
num_repeat = 5

#load data

if (data_name == "se_Quartz") {
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\SPOTlight_deconvolution_analysis-master\\analysis\\tool_benchmarking"
  )
  data = readRDS("se_quartz.rds")
  clust_vr = "nnet2"
  table(data$nnet2)
}
if (data_name == "PDAC-A_itai_processed") {
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\SPOTlight_deconvolution_analysis-master\\analysis\\pancreas_PDAC\\data"
  )
  data = readRDS("PDAC-A_itai_processed.rds")
  clust_vr = "annotation"
  table(data$clust_vr)
}
if (data_name == "PDAC-B_itai_processed") {
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\SPOTlight_deconvolution_analysis-master\\analysis\\pancreas_PDAC\\data"
  )
  data = readRDS("PDAC-B_itai_processed.rds")
  clust_vr = "annotation"
  table(data$clust_vr)
}
if (data_name == "cerebellum_singlecell_sn") {
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\data_RCTD\\RCTD_data\\mouse cerebellum"
  )
  data = readRDS("scRefSubsampled1000_cerebellum_singlenucleus.rds")
  clust_vr = "liger_ident_coarse"
  table(data$clust_vr)
}
if (data_name == "cerebellum_singlecell") {
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\data_RCTD\\RCTD_data\\mouse cerebellum"
  )
  data = readRDS("1000cellsSubsampled_cerebellum_singlecell.rds")
  clust_vr = "liger_ident_coarse"
  table(data$liger_ident_coarse)
  ##cell_type = names(table(data$liger_ident_coarse))[which(table(data$liger_ident_coarse) >
  ##                                                          200)]
  #data=data[,data$liger_ident_coarse%in%cell_type]
  
  data_counts = data@assays$RNA@counts
  #data_counts <- data_counts[!duplicated(rownames(data_counts)), ]
  #data_counts=data_counts[sample(nrow(data_counts),round(0.5*nrow(data_counts))),]
  cell_type_factor = data$liger_ident_coarse
  nnet2 = data.frame(cell_type_factor)
  colnames(nnet2) = "nnet2"
  rownames(nnet2) = colnames(data_counts)
  
  data = CreateSeuratObject(counts = data_counts,
                            assay = "RNA",
                            meta.data = nnet2)
  #data = data[, data$nnet2 %in% cell_type]
  clust_vr = "nnet2"
}
#test_spot=FALSE
#generate st data
for (i in 1:num_repeat) {
  print(sprintf("Generating %s synthetic test mixtures", n_syn_mixt))
  set.seed(i)

  test_spot_ls <- test_spot_fun(se_obj = data,
                                clust_vr = clust_vr,
                                n = n_syn_mixt)
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\fda_ST\\parameter_test\\cl_n\\st_data"
  )
  saveRDS(test_spot_ls, file = sprintf("test_spot_ls%s.rds", i))
  
}

#marker gene
if (marker_gene_access == 0) {
  se_sc <- Seurat::SCTransform(object = data,
                               assay = "RNA")
  # se_sc <- Seurat::SCTransform(object = data,
  #                              assay = "RNA",
  #                              verbose = FALSE)
  Seurat::Idents(se_sc) <- se_sc$nnet2
  marker_genes <- Seurat::FindAllMarkers(
    object = se_sc,
    assay = "SCT",
    slot = "data",
    min.pct = 0,
    only.pos = TRUE,
    logfc.threshold = 0
  )
  
  marker_genes %>% dplyr::count(cluster)
  
  marker_genes_filt <- marker_genes %>%
    dplyr::filter(pct.1 > 0.9)
  
  marker_genes_filt %>% dplyr::count(cluster)
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\fda_ST\\parameter_test\\cl_n\\markergene"
  )
  saveRDS(marker_genes, file = sprintf("%s_markergene.rds", data_name))
} else{
  setwd(
    "D:\\study\\undergraduate\\spolight\\code\\fda_ST\\parameter_test\\cl_n\\markergene"
  )
  marker_genes = readRDS(sprintf("%s_markergene.rds", data_name))
  
  marker_genes %>% dplyr::count(cluster)
  
  marker_genes_filt <- marker_genes %>%
    dplyr::filter(pct.1 > 0.9)
  
  marker_genes_filt = marker_genes %>%
    dplyr::filter(pct.2<0.1)
  
  marker_genes_filt %>% dplyr::count(cluster)
}

data_performance = data.frame(matrix(0, nrow = (4 * num_cl_n), ncol = 3))
colnames(data_performance) = c("cl_n", "value", "metric")
data_performance[, 3] = rep(c("accuracy", "sensitivity", "specificity", "F1"), num_cl_n)

data_performance_JSD = data.frame(matrix(0, nrow = (3 * num_repeat * num_cl_n), ncol =
                                           3))
colnames(data_performance_JSD) = c("cl_n", "value", "metric")
data_performance_JSD[, 3] = rep(c("0.25", "0.50", "0.75"), num_cl_n * num_repeat)

for (i in 1:num_cl_n) {
  cl_n = i * 10
  data_performance[(1 + (i - 1) * 4):(4 + (i - 1) * 4), 1] = rep(cl_n, 4)
  data_performance_JSD[(1 + (i - 1) * 3 * num_repeat):(3 * num_repeat +
                                                         (i - 1) * 3 * num_repeat), 1] = rep(cl_n, 3 * num_repeat)
  
  accuracy_sum = 0
  sensitivity_sum = 0
  specificity_sum = 0
  precision_sum = 0
  recall_sum = 0
  F1_sum = 0
  
  for (j in 1:num_repeat) {
    setwd(
      "D:\\study\\undergraduate\\spolight\\code\\fda_ST\\parameter_test\\cl_n\\st_data"
    )
    test_spot_ls = readRDS(file = sprintf("test_spot_ls%s.rds", j))
    test_spot_counts <- as.matrix(test_spot_ls[[1]])
    colnames(test_spot_counts) <-
      paste("mixt", 1:ncol(test_spot_counts), sep = "_")
    test_spot_metadata <- test_spot_ls[[2]]
    #test_spot_metadata=test_spot_metadata/rowSums(test_spot_metadata)
    
    se_sc_down <- downsample_se_obj(
      se_obj = data,
      clust_vr = clust_vr,
      cluster_markers = marker_genes,
      cl_n = cl_n,
      hvg = hvg
    )
    
    print("Deconvolute synthetic spots")
    
    geno = as.matrix(se_sc_down@assays$RNA@counts)
    ST.matrix = test_spot_counts
    ST.matrix = ST.matrix[rownames(geno), ]
    cell.type.factor = se_sc_down$nnet2
    
    decon_mtrx = fun_spotlight(
      geno = geno,
      ST.matrix = ST.matrix,
      cell.type.factor = cell.type.factor,
      nkonts = 150,
      ct_mode = "median",
      marker_gene_order = "1",
      cluster_marker = marker_genes
    )
    
    decon_mtrx = t(decon_mtrx)
    
    fun_spotlight_deconv <-
      decon_mtrx[, colnames(decon_mtrx) != "ress_ss"]
    ct_cols <- colnames(fun_spotlight_deconv)
    spatial_decon_syn <-
      test_synthetic_performance(
        test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
        spot_composition_mtrx = fun_spotlight_deconv[, ct_cols]
      )
    data_performance_JSD[1 + (j - 1) * 3 + (i - 1) * 3 * num_repeat, 2] =
      spatial_decon_syn$JSD[1]
    data_performance_JSD[2 + (j - 1) * 3 + (i - 1) * 3 * num_repeat, 2] =
      spatial_decon_syn$JSD[2]
    data_performance_JSD[3 + (j - 1) * 3 + (i - 1) * 3 * num_repeat, 2] =
      spatial_decon_syn$JSD[3]
    
    tp = spatial_decon_syn$TP
    tn = spatial_decon_syn$TN
    fp = spatial_decon_syn$FP
    fn = spatial_decon_syn$FN
    
    accuracy <- round((tp + tn) / (tp + tn + fp + fn), 2)
    sensitivity <- round(tp / (tp + fn), 2)
    specificity <- round(tn / (tn + fp), 2)
    precision <- round(tp / (tp + fp), 2)
    recall <- round(tp / (tp + fn), 2)
    F1 <- round(2 * ((precision * recall) / (precision + recall)), 2)
    
    accuracy_sum = accuracy_sum + accuracy
    sensitivity_sum = sensitivity_sum + sensitivity
    specificity_sum = specificity_sum + specificity
    precision_sum = precision_sum + precision
    recall_sum = recall_sum + recall
    F1_sum = F1_sum + F1
    
  }
  
  accuracy = accuracy_sum / num_repeat
  sensitivity = sensitivity_sum / num_repeat
  specificity = specificity_sum / num_repeat
  precision = precision_sum / num_repeat
  recall = recall_sum / num_repeat
  F1 = F1_sum / num_repeat
  
  data_performance[1 + (i - 1) * 4, 2] = accuracy
  data_performance[2 + (i - 1) * 4, 2] = sensitivity
  data_performance[3 + (i - 1) * 4, 2] = specificity
  data_performance[4 + (i - 1) * 4, 2] = F1
}
cl_n_list=seq(1,20)*10
colourCount=20
#plot
performance_plt <- data_performance %>%
  dplyr::mutate(
    cl_n = factor(x = cl_n,
                  levels = cl_n_list)
  ) %>%
  ggplot() +
  geom_point(aes(x = cl_n,
                 y = value,
                 color = cl_n),
             size = 5,
             alpha = 0.9) +
  ylim(0.25,1)+
  facet_wrap(. ~ metric) +
  labs(title = "Deconvolution parameter benchmarking(cell number)",
       x = "cl_n/num_cell_per_celltype",
       y = "Metric Value") +
  theme_classic()  +
  scale_fill_manual(values = colorRampPalette(brewer.pal(20, "Set3"))(colourCount))+
  theme(
    strip.text = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15)
  )

performance_plt
colnames(data_performance_JSD)=c("cl_n","value","quantile")
quan_list=c("0.25","0.50","0.75")
#scale_color_brewer(palette = "Set3")
# Tool = factor(x = Tool,
#               levels = tech_list)

performance_JSD_plt <- data_performance_JSD %>%
  dplyr::mutate(
    quantile=factor(x=quantile,levels = quan_list)
  ) %>%
  ggplot() +
  geom_point(aes(x = cl_n,
                 y = value,
                 color = quantile),
             size = 1,
             alpha = 0.9) +
  ylim(0,0.5)+
  labs(title = "Deconvolution parameter benchmarking(cell number)",
       x = "cl_n",
       y = "JSD") +
  theme_classic()  +
  scale_fill_manual(values = colorRampPalette(brewer.pal(13, "Set3"))(colourCount))+
  theme(
    strip.text = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(size = 12, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15)
  )

performance_JSD_plt
#facet_wrap(. ~ metric)


