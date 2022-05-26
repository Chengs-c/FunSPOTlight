library(SPOTlight)
library(fda)
library(magrittr)
setwd("D:\\study\\undergraduate\\spolight\\code\\fda_ST\\fun_spotlight")
source("fun_spotlight.r")

spatial_decon_syn_assessment_fun_spotlight_fun=function(se_sc,
                                                        clust_vr,
                                                        n_syn_mixt = 1000,
                                                        cl_n = 100,
                                                        hvg = 0,
                                                        ntop = NULL,
                                                        min_cont = 0.01,
                                                        assay = "RNA",
                                                        slot = "counts"){
  
  
  ##########################
  ### Generate test data ###
  ##########################
  print(sprintf("Generating %s synthetic test mixtures", n_syn_mixt))
  test_spot_ls <- test_spot_fun(se_obj = se_sc,
                                clust_vr = clust_vr,
                                n = n_syn_mixt)
  
  test_spot_counts <- as.matrix(test_spot_ls[[1]])
  colnames(test_spot_counts) <- paste("mixt", 1:ncol(test_spot_counts), sep = "_")
  test_spot_metadata <- test_spot_ls[[2]]
  
  ####################
  ### Downsampling ###
  ####################
  # Downsample number of genes and number of samples
  print("Downsampling genes and cells")
  
  se_sc <- Seurat::SCTransform(object = se_sc,
                                   assay = "RNA",
                                   verbose = FALSE)
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
  
  
  
  # Downsample scRNAseq to select gene set and number of cells to train the model
  se_sc_down <- downsample_se_obj(se_obj = se_sc,
                                  clust_vr = clust_vr,
                                  cluster_markers = marker_genes,
                                  cl_n = cl_n,
                                  hvg = hvg)
  
  
  print("Deconvolute synthetic spots")
  
  geno = as.matrix(se_sc_down@assays$RNA@counts)
  ST.matrix = test_spot_counts
  ST.matrix = ST.matrix[rownames(geno),]
  cell.type.factor = se_sc_down$nnet2
  
  decon_mtrx = fun_spotlight(geno = geno,
                             ST.matrix = ST.matrix,
                             cell.type.factor = cell.type.factor)
  
  decon_mtrx = t(decon_mtrx)
  
  fun_spotlight_deconv <-
    decon_mtrx[, colnames(decon_mtrx) != "ress_ss"]
  ct_cols <- colnames(fun_spotlight_deconv)
  raw_statistics_ls <- test_synthetic_performance(test_spots_metadata_mtrx = as.matrix(test_spot_metadata[, ct_cols]),
                                                  spot_composition_mtrx = fun_spotlight_deconv[, ct_cols])
  
  return(list("fun_spotlight_deconv" = fun_spotlight_deconv,
              "stats" = raw_statistics_ls,
              "test_spot"=test_spot_metadata))
}
  
  
  
  
