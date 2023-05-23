<<<<<<< HEAD
library(scrattch.hicat)
library(Matrix)
library(matrixStats)
library(CYANO)
library(ROGUE)
library(tibble)
library(ktplots)
library(beeswarm)
library(WGCNA)
library(edgeR)   
library(feather)
library(dendextend)
library(ggplot2)
library(dplyr)
library(Seurat)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)
library(knitr)
library(parallel)
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(hdWGCNA)

set.seed(12345)

file_path_root = "/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/mouseALM/data1/"
data = read.csv(paste0(file_path_root, "data.csv"),row.names=1)
data_meta = readRDS(paste0(file_path_root,"classificiation_resault.rds"))
new_data = data[,rownames(data_meta)]
mouse_data <- CreateSeuratObject(counts = new_data, min.cells = 0, min.features = 0, project = "example")
mouse_data <- AddMetaData(mouse_data, data_meta)
Idents(mouse_data) = mouse_data$cluster_label
mouse_data <- NormalizeData(object = mouse_data)
mouse_data <- FindVariableFeatures(object = mouse_data, selection.method = 'vst', nfeatures = 2000)
mouse_data <- ScaleData(mouse_data)

mouse_data <- RunPCA(mouse_data, features = VariableFeatures(object = mouse_data),npcs=100)
seurat_obj <- SetupForWGCNA(
  mouse_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("class_label","subclass_label","cluster_label"), # specify the columns in seurat_obj@meta.data to group by
  k = 3, # nearest-neighbors parameter
  max_shared = 2, # maximum number of shared cells between two metacells
  ident.group = 'cluster_label',
  min_cells = 4# set the Idents of the metacell seurat object
)
metacell_obj <- GetMetacellObject(seurat_obj)
length(unique(metacell_obj$cells_merged))
new_hdwgcna_data1 = metacell_obj@assays$RNA
new_hdwgcna_meta1 = metacell_obj@meta.data
new_hdwgcna_data_new1 = as.matrix(new_hdwgcna_data1@data)
saveRDS(new_hdwgcna_data_new1, paste0(file_path_root,"wgcna_pse_data.rds"))
saveRDS(new_hdwgcna_meta1, paste0(file_path_root, "wgcna_pse_meta.rds"))

source("/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/get_file_name_function_for_hdwgcna.R")
cell_number=10
anno_meta = readRDS(paste0(file_path_root, "classificiation_resault.rds"))
result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",TRUE,paste0(file_path_root,"wgcna_pse_data.rds"),paste0(file_path_root,"wgcna_pse_meta.rds"))
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "wgcna_pse_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "wgcna_pse_roc_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "wgcna_pse_roc_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "wgcna_pse_roc.csv"),
                   paste0(file_path_root, "wgcna_pse_roc_old.csv"),paste0(file_path_root, "wgcna_pse_roc_raw.csv"))




set.seed(12345)

file_path_root = "/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/mouseALM/data2/"
data = read.csv(paste0(file_path_root, "data.csv"),row.names=1)
data_meta = readRDS(paste0(file_path_root,"classificiation_resault.rds"))
new_data = data[,rownames(data_meta)]
mouse_data <- CreateSeuratObject(counts = new_data, min.cells = 0, min.features = 0, project = "example")
mouse_data <- AddMetaData(mouse_data, data_meta)
Idents(mouse_data) = mouse_data$cluster_label
mouse_data <- NormalizeData(object = mouse_data)
mouse_data <- FindVariableFeatures(object = mouse_data, selection.method = 'vst', nfeatures = 2000)
mouse_data <- ScaleData(mouse_data)

mouse_data <- RunPCA(mouse_data, features = VariableFeatures(object = mouse_data),npcs=100)
seurat_obj <- SetupForWGCNA(
  mouse_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("class_label","subclass_label","cluster_label"), # specify the columns in seurat_obj@meta.data to group by
  k = 3, # nearest-neighbors parameter
  max_shared = 2, # maximum number of shared cells between two metacells
  ident.group = 'cluster_label',
  min_cells = 4# set the Idents of the metacell seurat object
)
metacell_obj <- GetMetacellObject(seurat_obj)
length(unique(metacell_obj$cells_merged))
new_hdwgcna_data1 = metacell_obj@assays$RNA
new_hdwgcna_meta1 = metacell_obj@meta.data
new_hdwgcna_data_new1 = as.matrix(new_hdwgcna_data1@data)
saveRDS(new_hdwgcna_data_new1, paste0(file_path_root,"wgcna_pse_data.rds"))
saveRDS(new_hdwgcna_meta1, paste0(file_path_root, "wgcna_pse_meta.rds"))

source("/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/get_file_name_function_for_hdwgcna.R")
cell_number=10
anno_meta = readRDS(paste0(file_path_root, "classificiation_resault.rds"))
result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",TRUE,paste0(file_path_root,"wgcna_pse_data.rds"),paste0(file_path_root,"wgcna_pse_meta.rds"))
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "wgcna_pse_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "wgcna_pse_roc_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "wgcna_pse_roc_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "wgcna_pse_roc.csv"),
                   paste0(file_path_root, "wgcna_pse_roc_old.csv"),paste0(file_path_root, "wgcna_pse_roc_raw.csv"))

=======
library(scrattch.hicat)
library(Matrix)
library(matrixStats)
library(CYANO)
library(ROGUE)
library(tibble)
library(ktplots)
library(beeswarm)
library(WGCNA)
library(edgeR)   
library(feather)
library(dendextend)
library(ggplot2)
library(dplyr)
library(Seurat)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)
library(knitr)
library(parallel)
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(hdWGCNA)

set.seed(12345)

file_path_root = "/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/mouseALM/data1/"
data = read.csv(paste0(file_path_root, "data.csv"),row.names=1)
data_meta = readRDS(paste0(file_path_root,"classificiation_resault.rds"))
new_data = data[,rownames(data_meta)]
mouse_data <- CreateSeuratObject(counts = new_data, min.cells = 0, min.features = 0, project = "example")
mouse_data <- AddMetaData(mouse_data, data_meta)
Idents(mouse_data) = mouse_data$cluster_label
mouse_data <- NormalizeData(object = mouse_data)
mouse_data <- FindVariableFeatures(object = mouse_data, selection.method = 'vst', nfeatures = 2000)
mouse_data <- ScaleData(mouse_data)

mouse_data <- RunPCA(mouse_data, features = VariableFeatures(object = mouse_data),npcs=100)
seurat_obj <- SetupForWGCNA(
  mouse_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("class_label","subclass_label","cluster_label"), # specify the columns in seurat_obj@meta.data to group by
  k = 3, # nearest-neighbors parameter
  max_shared = 2, # maximum number of shared cells between two metacells
  ident.group = 'cluster_label',
  min_cells = 4# set the Idents of the metacell seurat object
)
metacell_obj <- GetMetacellObject(seurat_obj)
length(unique(metacell_obj$cells_merged))
new_hdwgcna_data1 = metacell_obj@assays$RNA
new_hdwgcna_meta1 = metacell_obj@meta.data
new_hdwgcna_data_new1 = as.matrix(new_hdwgcna_data1@data)
saveRDS(new_hdwgcna_data_new1, paste0(file_path_root,"wgcna_pse_data.rds"))
saveRDS(new_hdwgcna_meta1, paste0(file_path_root, "wgcna_pse_meta.rds"))

source("/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/get_file_name_function_for_hdwgcna.R")
cell_number=10
anno_meta = readRDS(paste0(file_path_root, "classificiation_resault.rds"))
result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",TRUE,paste0(file_path_root,"wgcna_pse_data.rds"),paste0(file_path_root,"wgcna_pse_meta.rds"))
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "wgcna_pse_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "wgcna_pse_roc_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "wgcna_pse_roc_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "wgcna_pse_roc.csv"),
                   paste0(file_path_root, "wgcna_pse_roc_old.csv"),paste0(file_path_root, "wgcna_pse_roc_raw.csv"))




set.seed(12345)

file_path_root = "/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/mouseALM/data2/"
data = read.csv(paste0(file_path_root, "data.csv"),row.names=1)
data_meta = readRDS(paste0(file_path_root,"classificiation_resault.rds"))
new_data = data[,rownames(data_meta)]
mouse_data <- CreateSeuratObject(counts = new_data, min.cells = 0, min.features = 0, project = "example")
mouse_data <- AddMetaData(mouse_data, data_meta)
Idents(mouse_data) = mouse_data$cluster_label
mouse_data <- NormalizeData(object = mouse_data)
mouse_data <- FindVariableFeatures(object = mouse_data, selection.method = 'vst', nfeatures = 2000)
mouse_data <- ScaleData(mouse_data)

mouse_data <- RunPCA(mouse_data, features = VariableFeatures(object = mouse_data),npcs=100)
seurat_obj <- SetupForWGCNA(
  mouse_data,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("class_label","subclass_label","cluster_label"), # specify the columns in seurat_obj@meta.data to group by
  k = 3, # nearest-neighbors parameter
  max_shared = 2, # maximum number of shared cells between two metacells
  ident.group = 'cluster_label',
  min_cells = 4# set the Idents of the metacell seurat object
)
metacell_obj <- GetMetacellObject(seurat_obj)
length(unique(metacell_obj$cells_merged))
new_hdwgcna_data1 = metacell_obj@assays$RNA
new_hdwgcna_meta1 = metacell_obj@meta.data
new_hdwgcna_data_new1 = as.matrix(new_hdwgcna_data1@data)
saveRDS(new_hdwgcna_data_new1, paste0(file_path_root,"wgcna_pse_data.rds"))
saveRDS(new_hdwgcna_meta1, paste0(file_path_root, "wgcna_pse_meta.rds"))

source("/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/get_file_name_function_for_hdwgcna.R")
cell_number=10
anno_meta = readRDS(paste0(file_path_root, "classificiation_resault.rds"))
result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",TRUE,paste0(file_path_root,"wgcna_pse_data.rds"),paste0(file_path_root,"wgcna_pse_meta.rds"))
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "wgcna_pse_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "wgcna_pse_roc_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "wgcna_pse_roc_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "wgcna_pse_roc.csv"),
                   paste0(file_path_root, "wgcna_pse_roc_old.csv"),paste0(file_path_root, "wgcna_pse_roc_raw.csv"))

>>>>>>> 3e2f6a6 (add large file)
