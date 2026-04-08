rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(dplyr)
library(Seurat)
library(patchwork)
library(liger) # install_github('MacoskoLab/liger')
library(SeuratWrappers) # remotes::install_github('satijalab/seurat-wrappers')


# create count matrices

    D56_iCtrl = Read10X(data.dir=paste0(FASTQ_folder,"D56_iCtrl/D56_iCtrl/outs/filtered_feature_bc_matrix"))
    D56_Mut = Read10X(data.dir=paste0(FASTQ_folder,"D56_Mut/D56_Mut/outs/filtered_feature_bc_matrix"))
    D70_iCtrl = Read10X(data.dir=paste0(FASTQ_folder,"D70_iCtrl/D70_iCtrl/outs/filtered_feature_bc_matrix"))
    D70_Mut = Read10X(data.dir=paste0(FASTQ_folder,"D70_Mut/D70_Mut/outs/filtered_feature_bc_matrix"))
    D100_iCtrl = Read10X(data.dir=paste0(FASTQ_folder,"D100_iCtrl/D100_iCtrl/outs/filtered_feature_bc_matrix"))
    D100_Mut = Read10X(data.dir=paste0(FASTQ_folder,"D100_Mut/D100_Mut/outs/filtered_feature_bc_matrix"))

# merge to one seurat object

    datExpr <- merge(D56_iCtrl, y = c(D56_Mut, D70_iCtrl, D70_Mut, D100_iCtrl, D100_Mut),
                          add.cell.ids = c("D56_iCtrl", "D56_Mut", "D70_iCtrl", "D70_Mut", "D100_iCtrl", "D100_Mut"),
                          project = "all_combined")
    datExpr[["percent.mt"]] <- PercentageFeatureSet(datExpr, pattern = "^MT-")
    rb.genes <- rownames(datExpr)[grep("^RP[SL]",rownames(datExpr))]
    datExpr[["percent.ribo"]] <- PercentageFeatureSet(datExpr, features = rb.genes)

# filtering 

    datMeta=as.data.frame(as.matrix(datExpr@meta.data))
    MaxGene=mean(as.numeric(datMeta$nFeature_RNA)) + c(3*sd(as.numeric(datMeta$nFeature_RNA)))

    datExpr_unfilt = datExpr
    datExpr_unfilt
        # An object of class Seurat 
        # 27379 features across 55870 samples within 1 assay 
        # Active assay: RNA (27379 features, 0 variable features)


    datExpr <- subset(datExpr, subset = nFeature_RNA > 500 & nFeature_RNA < MaxGene & percent.mt < 10)
    datExpr
        # An object of class Seurat 
        # 27379 features across 49942 samples within 1 assay 
        # Active assay: RNA (27379 features, 0 variable features)

# normalization, batch correction, and clustering

    datExpr <- NormalizeData(datExpr)
    datExpr <- FindVariableFeatures(datExpr)
    datExpr <- ScaleData(datExpr, split.by = "orig.ident", do.center = FALSE)
    datExpr <- RunOptimizeALS(datExpr, k = 20, lambda = 5, split.by = "orig.ident")
        # Finished in 5.501731 mins, 30 iterations.
        # Max iterations set: 30.
        # Final objective delta: 7.629359e-06.
        # Best results with seed 1.
        # Warning: No columnames present in cell embeddings, setting to 'riNMF_1:20'
    datExpr <- RunQuantileNorm(datExpr, split.by = "orig.ident")
    datExpr <- FindNeighbors(datExpr, reduction = "iNMF", dims = 1:20)
    datExpr <- FindClusters(datExpr, resolution = 0.3) # 14 clusters
    datExpr <- RunUMAP(datExpr, dims = 1:ncol(datExpr[["iNMF"]]), reduction = "iNMF")

# rename clusters

    Idents(datExpr) <- "seurat_clusters"
    new.cluster.ids <- c("CFuPNs","Inh.Neurons","oRGs","CPNs","RGCs","Cycling_Prog.",
                         "Unknown","IPCs","CPNs","CFuPNs","Inh.Neurons","Unknown","Cajal-Retzius","CFuPNs")

    names(new.cluster.ids) <- levels(datExpr)
    datExpr <- RenameIdents(datExpr, new.cluster.ids)
    datExpr[["new.cluster.ids"]] <- Idents(datExpr)

    table(Idents(datExpr),datExpr$Genotype)
        # iCtrl Mut
        # CFuPNs        4718 5612
        # Inh.Neurons   3411 3475
        # oRGs          3076 2204
        # CPNs          3877 4682
        # RGCs          3493 1363
        # Cycling_Prog. 2269 1736
        # Unknown       3295 1929
        # IPCs          2212 1600
        # Cajal-Retzius  584  406

# avg expression for heatmap visualization

    datExpravg = AverageExpression(datExpr, return.seurat = T, add.ident = 'Genotype') 

# Differential Expression and GO analysis

library(gprofiler2)

    Idents(datExpr) <- "new.cluster.ids"
    datExpr$cluster_genotype <- paste(Idents(datExpr), datExpr$Genotype, sep = "_")
    Idents(datExpr) <- "cluster_genotype"
    clusterNum=data.frame(table(datExpr$new.cluster.ids))

    for (i in 1:length(clusterNum$Var1)) {
          print(paste0("Starting cluster ",clusterNum$Var1[i]))
          mymarkers <- FindMarkers(datExpr, ident.1 = paste0(clusterNum$Var1[i],"_Rett"), ident.2 = paste0(clusterNum$Var1[i],"_Ctrl"), verbose = TRUE)
          mymarkers = mymarkers[mymarkers$p_val_adj<0.05,]
          write.csv(mymarkers,file=paste0(folder,"_",clusterNum$Var1[i],"_DEGs_Rett_vs_Ctrl.csv"))

          up.list=mymarkers[mymarkers$avg_logFC > 0,]
          up.query=rownames(up.list)
          print(length(up.query))
          up.results = gost(up.query, organism = "hsapiens",ordered_query = TRUE,
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                          measure_underrepresentation = FALSE, evcodes = FALSE,
                          user_threshold = 0.05, correction_method = c("bonferroni"),
                          user_threshold = 0.05, correction_method = c("bonferroni"), custom_bg = rownames(datExpr),
                          numeric_ns = "", sources = c("GO:BP","GO:MF","GO:CC","HP"), as_short_link = TRUE)
          print(up.results)

          down.list=mymarkers[mymarkers$avg_logFC < 0,]
          down.list=down.list[order(down.list$avg_logFC,decreasing=F),]
          down.query=rownames(down.list)
          print(length(down.query))
          down.results = gost(down.query, organism = "hsapiens",ordered_query = TRUE,
                            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                            measure_underrepresentation = FALSE, evcodes = FALSE,
                            user_threshold = 0.05, correction_method = c("bonferroni"), custom_bg = rownames(datExpr),
                            numeric_ns = "", sources = c("GO:BP","GO:MF","GO:CC","HP"), as_short_link = TRUE)
          print(down.results)
    }
 
# subset and recluster inhibitory neuron cluster

mycells=WhichCells(datExpr, idents = "Inh.Neurons")
datExprInh = subset(datExpr, cells=mycells)
Idents(datExprInh) <- "orig.ident"
datExprInh <- FindNeighbors(datExprInh, reduction = "iNMF", dims = 1:20)
datExprInh <- FindClusters(datExprInh, resolution = 0.3) # 7 clusters
datExprInh <- RunUMAP(datExprInh, dims = 1:ncol(datExprInh[["iNMF"]]), reduction = "iNMF")

Inhavgexp = AverageExpression(datExprInh, return.seurat = T, add.ident = 'Genotype') # avg for heatmap visualization


# scDC: Single cell differential composition analysis https://sydneybiox.github.io/scDC/index.html

# install.packages(c("boot", "class", "cli", "codetools", "digest", "glue", "IRkernel", 
#                    "jsonlite", "KernSmooth", "lattice", "MASS", "Matrix", "mgcv", "nlme", "nnet", "pbdZMQ", 
#                    "pillar", "rlang", "spatial", "survival", "vctrs", "DescTools", "lme4", "reshape2", "ggridges", "lme4", "mice"))
# BiocManager::install(c("BiocNeighbors","scater","scran"))
# devtools::install_github("taiyunkim/scClustBench")
# devtools::install_github("SydneyBioX/scDC")

library(scDC)

# for all clusters

    cellTypes = datExpr$new.cluster.ids # cluster assignments 
    subject = datExpr$Genotype # genotype assignments

    res_percentile = scDC_noClustering(cellTypes, subject, calCI = TRUE, calCI_method = "percentile", verbose = TRUE)

    res_percentile <- res$results

# for inhibitory cluster (reclustered)
  
  # select cells from inhibitory cluster
      mycells=WhichCells(datExpr, idents = "Inh.Neurons") 
      datExprInh = subset(datExpr, cells=mycells)
      Idents(datExprInh) <- "orig.ident"
  
  # recluster 
      datExprInh <- FindNeighbors(datExprInh, reduction = "iNMF", dims = 1:20)
      datExprInh <- FindClusters(datExprInh, resolution = 0.3)
      datExprInh <- RunUMAP(datExprInh, dims = 1:ncol(datExprInh[["iNMF"]]), reduction = "iNMF")
  
  cellTypesInh = datExprInh$seurat_clusters # cluster assignments 
  subjectInh = datExprInh$Genotype # genotype assignments

  res_percentile_Inh = scDC_noClustering(cellTypesInh, subjectInh, calCI = TRUE, calCI_method = "percentile", verbose = TRUE)

  res_percentile_Inh <- res_percentile_Inh$results
