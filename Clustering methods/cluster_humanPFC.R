<<<<<<< HEAD
#consensus 1
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
library(CYANO)
file_path_root = "/root/datasets/humanPFC/"
de.param <- de_param(padj.th=0.01,
                     q1.th       = 0.3, 
                     q.diff.th   = 0.7, 
                     de.score.th = 100,
                     min.cells=20)


data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"data1/run_cluster"))
result_1 <- RunConsensus1(norm.dat,
                          de.param = de.param,
                          override = TRUE,
                          output_dir = paste0(file_path_root,"data1/run_cluster"),
                          mc.cores = 1)
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
saveRDS(result_1, paste0(file_path_root,"data1/consesus_1_result_1.rds")) 
saveRDS(group_meta, paste0(file_path_root,"data1/group_meta_1_new.rds")) 
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_consensus1.rds"))



data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"data2/run_cluster"))
result_1 <- RunConsensus1(norm.dat,
                          de.param = de.param,
                          override = TRUE,
                          output_dir = paste0(file_path_root,"data2/run_cluster"),
                          mc.cores = 1)
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
saveRDS(result_1, paste0(file_path_root,"data2/consesus_1_result_2.rds")) 
saveRDS(group_meta, paste0(file_path_root,"data2/group_meta_2_new.rds")) 
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_consensus1.rds"))






#Seurat
library(dplyr)
library(Seurat)
library(patchwork)

data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
pbmc <- CreateSeuratObject(counts = data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:14)
pbmc <- FindClusters(pbmc, resolution = 10)
write.csv(pbmc@meta.data, paste0(file_path_root, "data1/seurat_result_10.csv"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_seurat.rds"))




data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
pbmc <- CreateSeuratObject(counts = data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:14)
pbmc <- FindClusters(pbmc, resolution = 10)
write.csv(pbmc@meta.data, paste0(file_path_root, "data2/seurat_result_10.csv"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_seurat.rds"))






#monocle3
library(monocle3)
library(dplyr) # imported for some downstream data manipulation
data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
pd_1 = data.frame(colnames(data))
rownames(pd_1) = colnames(data)
fd_1 = data.frame(rownames(data))
rownames(fd_1) = rownames(data)
colnames(fd_1)=c("gene_short_name")
pd <- new("AnnotatedDataFrame", data = pd_1)
fd <- new("AnnotatedDataFrame", data = fd_1)
# cds <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
#                phenoData = pd,
#                featureData = fd,
#                lowerDetectionLimit = 0.5,
#                expressionFamily = negbinomial.size())
cds <- new_cell_data_set(expression_data=as.matrix(data),
                         cell_metadata=pd_1,
                         gene_metadata=fd_1)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, umap.fast_sgd=TRUE,preprocess_method = 'PCA')
cds <- cluster_cells(cds, k=10,cluster_method="louvain")
saveRDS(cds, paste0(file_path_root, "data1/monocle3_10.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_monocle.rds"))





data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
pd_1 = data.frame(colnames(data))
rownames(pd_1) = colnames(data)
fd_1 = data.frame(rownames(data))
rownames(fd_1) = rownames(data)
colnames(fd_1)=c("gene_short_name")
pd <- new("AnnotatedDataFrame", data = pd_1)
fd <- new("AnnotatedDataFrame", data = fd_1)
# cds <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
#                phenoData = pd,
#                featureData = fd,
#                lowerDetectionLimit = 0.5,
#                expressionFamily = negbinomial.size())
cds <- new_cell_data_set(expression_data=as.matrix(data),
                         cell_metadata=pd_1,
                         gene_metadata=fd_1)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, umap.fast_sgd=TRUE,preprocess_method = 'PCA')
cds <- cluster_cells(cds, k=10,cluster_method="louvain")
saveRDS(cds, paste0(file_path_root, "data2/monocle3_10.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_monocle.rds"))







#SC3
library(SingleCellExperiment)
library(SC3)
library(scater)

data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
sce  <- SingleCellExperiment(assays = list(counts = data,
                                           logcounts = log2(as.matrix(data) + 1)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3(sce, ks = 100:150,n_cores=20)
sce <- sc3_run_svm(sce, ks = 100:150)
saveRDS(sce, paste0(file_path_root, "data1/SC3.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_sc3.rds"))


data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
sce  <- SingleCellExperiment(assays = list(counts = data,
                                           logcounts = log2(as.matrix(data) + 1)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3(sce, ks = 100:150,n_cores=20)
sce <- sc3_run_svm(sce, ks = 100:150)
saveRDS(sce, paste0(file_path_root, "data2/SC3.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_sc3.rds"))

#scCCESS-kmeans
library(scCCESS)
library(SingleCellExperiment)
library(reticulate)
library(tensorflow)
#use_virtualenv("r-reticulate")
data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
dat <- SingleCellExperiment(assays = list(counts = data))
dat <- counts(dat)

dat=prefilter(dat)
k = estimate_k(dat,
               seed = 1, 
               cluster_func = function(x,centers) { 
                 set.seed(42);
                 kmeans(x, centers)
               },
               criteria_method = "NMI",
               krange = 60:100, ensemble_sizes = 10,
               cores = 10
)
cluster = ensemble_cluster(dat, 
                           seed = 1, 
                           cluster_func = function(x) {
                             set.seed(1)
                             kmeans(x, centers = k$ngroups)
                           }, 
                           cores = 10, 
                           genes_as_rows = T, 
                           ensemble_sizes = 10, 
                           verbose = 0, 
                           scale = F, 
                           batch_size = 64
)
saveRDS(cluster,paste0(file_path_root, "data1/scCCESS_kmeans.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_kmeans.rds"))

data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
dat <- SingleCellExperiment(assays = list(counts = data))
dat <- counts(dat)

dat=prefilter(dat)
k = estimate_k(dat,
               seed = 1, 
               cluster_func = function(x,centers) { 
                 set.seed(42);
                 kmeans(x, centers)
               },
               criteria_method = "NMI",
               krange = 60:100, ensemble_sizes = 10,
               cores = 10
)

cluster = ensemble_cluster(dat, 
                           seed = 1, 
                           cluster_func = function(x) {
                             set.seed(1)
                             kmeans(x, centers = k$ngroups)
                           }, 
                           cores = 10, 
                           genes_as_rows = T, 
                           ensemble_sizes = 10, 
                           verbose = 0, 
                           scale = F, 
                           batch_size = 64
)
saveRDS(cluster,paste0(file_path_root, "data2/scCCESS_kmeans.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_kmeans.rds"))


#scCCESS-SIMLR
library(scCCESS)
library(SIMLR)

data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
dat <- SingleCellExperiment(assays = list(counts = data))
dat <- counts(dat)

k = estimate_k(dat,
               seed = 1, 
               cluster_func = function(x,centers) {
                 set.seed(42);
                 SIMLR_Large_Scale(t(x), c=centers,kk=15)
               },
               criteria_method = "NMI",
               krange = 60:100, ensemble_sizes = 10,
               cores = 10
)


cluster = ensemble_cluster(dat, 
                           seed = 1, 
                           cluster_func = function(x) {
                             set.seed(1)
                             SIMLR_Large_Scale(t(x), c=k$ngroups,kk=15)
                           }, 
                           cores = 10, 
                           genes_as_rows = T, 
                           ensemble_sizes = 10, 
                           verbose = 0, 
                           scale = F, 
                           batch_size = 64
)
saveRDS(cluster,paste0(file_path_root, "data1/scCCESS_simlr.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_simlr.rds"))


data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
dat <- SingleCellExperiment(assays = list(counts = data))
dat <- counts(dat)

k = estimate_k(dat, 
               seed = 1, 
               cluster_func = function(x,centers) {
                 set.seed(42);
                 SIMLR_Large_Scale(t(x), c=centers,kk=15)
               },
               criteria_method = "NMI",
               krange = 60:100, ensemble_sizes = 10,
               cores = 10
)


cluster = ensemble_cluster(dat, 
                           seed = 1, 
                           cluster_func = function(x) {
                             set.seed(1)
                             SIMLR_Large_Scale(t(x), c=k$ngroups,kk=15)
                           }, 
                           cores = 10, 
                           genes_as_rows = T, 
                           ensemble_sizes = 10, 
                           verbose = 0, 
                           scale = F, 
                           batch_size = 64
)
saveRDS(cluster,paste0(file_path_root, "data2/scCCESS_simlr.rds"))

end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_simlr.rds"))








#consensus 100
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



data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
#data = data1_new
# data = aggregate(data[,-1], by = list(data[,1]), mean)s
# rownames(data) = data[,1]
# data = data[,-1]
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"run_cluster2"))
result_1 <- run_consensus_clust(norm.dat,
                                niter = 100,
                                de.param = de.param,
                                dim.method = "pca",
                                output_dir = paste0(file_path_root,"run_cluster2"),
                                mc.cores = 10)



group_meta = data.frame(names(result_1$cl.result$cl),result_1$cl.result$cl)
colnames(group_meta) = c("cell","group")

saveRDS(result_1, paste0(file_path_root,"data1/consesus_100_result_1.rds")) 
saveRDS(group_meta, paste0(file_path_root,"data1/group_meta_100.rds")) 
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_consensus100.rds"))



data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
#data = data1_new
# data = aggregate(data[,-1], by = list(data[,1]), mean)s
# rownames(data) = data[,1]
# data = data[,-1]
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"run_cluster2"))
result_1 <- run_consensus_clust(norm.dat,
                                niter = 100,
                                de.param = de.param,
                                dim.method = "pca",
                                output_dir = paste0(file_path_root,"run_cluster2"),
                                mc.cores = 10)



group_meta = data.frame(names(result_1$cl.result$cl),result_1$cl.result$cl)
colnames(group_meta) = c("cell","group")

saveRDS(result_1, paste0(file_path_root,"data2/consesus_100_result_2.rds")) 
saveRDS(group_meta, paste0(file_path_root,"data2/group_meta_100.rds")) 
end <- Sys.time()
runningtime <- end-start
=======
#consensus 1
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
library(CYANO)
file_path_root = "/root/datasets/humanPFC/"
de.param <- de_param(padj.th=0.01,
                     q1.th       = 0.3, 
                     q.diff.th   = 0.7, 
                     de.score.th = 100,
                     min.cells=20)


data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"data1/run_cluster"))
result_1 <- RunConsensus1(norm.dat,
                          de.param = de.param,
                          override = TRUE,
                          output_dir = paste0(file_path_root,"data1/run_cluster"),
                          mc.cores = 1)
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
saveRDS(result_1, paste0(file_path_root,"data1/consesus_1_result_1.rds")) 
saveRDS(group_meta, paste0(file_path_root,"data1/group_meta_1_new.rds")) 
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_consensus1.rds"))



data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"data2/run_cluster"))
result_1 <- RunConsensus1(norm.dat,
                          de.param = de.param,
                          override = TRUE,
                          output_dir = paste0(file_path_root,"data2/run_cluster"),
                          mc.cores = 1)
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
saveRDS(result_1, paste0(file_path_root,"data2/consesus_1_result_2.rds")) 
saveRDS(group_meta, paste0(file_path_root,"data2/group_meta_2_new.rds")) 
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_consensus1.rds"))






#Seurat
library(dplyr)
library(Seurat)
library(patchwork)

data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
pbmc <- CreateSeuratObject(counts = data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:14)
pbmc <- FindClusters(pbmc, resolution = 10)
write.csv(pbmc@meta.data, paste0(file_path_root, "data1/seurat_result_10.csv"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_seurat.rds"))




data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
pbmc <- CreateSeuratObject(counts = data, project = "pbmc3k", min.cells = 3, min.features = 200)

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:14)
pbmc <- FindClusters(pbmc, resolution = 10)
write.csv(pbmc@meta.data, paste0(file_path_root, "data2/seurat_result_10.csv"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_seurat.rds"))






#monocle3
library(monocle3)
library(dplyr) # imported for some downstream data manipulation
data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
pd_1 = data.frame(colnames(data))
rownames(pd_1) = colnames(data)
fd_1 = data.frame(rownames(data))
rownames(fd_1) = rownames(data)
colnames(fd_1)=c("gene_short_name")
pd <- new("AnnotatedDataFrame", data = pd_1)
fd <- new("AnnotatedDataFrame", data = fd_1)
# cds <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
#                phenoData = pd,
#                featureData = fd,
#                lowerDetectionLimit = 0.5,
#                expressionFamily = negbinomial.size())
cds <- new_cell_data_set(expression_data=as.matrix(data),
                         cell_metadata=pd_1,
                         gene_metadata=fd_1)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, umap.fast_sgd=TRUE,preprocess_method = 'PCA')
cds <- cluster_cells(cds, k=10,cluster_method="louvain")
saveRDS(cds, paste0(file_path_root, "data1/monocle3_10.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_monocle.rds"))





data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
pd_1 = data.frame(colnames(data))
rownames(pd_1) = colnames(data)
fd_1 = data.frame(rownames(data))
rownames(fd_1) = rownames(data)
colnames(fd_1)=c("gene_short_name")
pd <- new("AnnotatedDataFrame", data = pd_1)
fd <- new("AnnotatedDataFrame", data = fd_1)
# cds <- newCellDataSet(as(as.matrix(data), "sparseMatrix"),
#                phenoData = pd,
#                featureData = fd,
#                lowerDetectionLimit = 0.5,
#                expressionFamily = negbinomial.size())
cds <- new_cell_data_set(expression_data=as.matrix(data),
                         cell_metadata=pd_1,
                         gene_metadata=fd_1)
cds <- preprocess_cds(cds, num_dim = 100)
cds <- reduce_dimension(cds, umap.fast_sgd=TRUE,preprocess_method = 'PCA')
cds <- cluster_cells(cds, k=10,cluster_method="louvain")
saveRDS(cds, paste0(file_path_root, "data2/monocle3_10.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_monocle.rds"))







#SC3
library(SingleCellExperiment)
library(SC3)
library(scater)

data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
sce  <- SingleCellExperiment(assays = list(counts = data,
                                           logcounts = log2(as.matrix(data) + 1)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3(sce, ks = 100:150,n_cores=20)
sce <- sc3_run_svm(sce, ks = 100:150)
saveRDS(sce, paste0(file_path_root, "data1/SC3.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_sc3.rds"))


data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
sce  <- SingleCellExperiment(assays = list(counts = data,
                                           logcounts = log2(as.matrix(data) + 1)))
rowData(sce)$feature_symbol <- rownames(sce)
sce <- sc3(sce, ks = 100:150,n_cores=20)
sce <- sc3_run_svm(sce, ks = 100:150)
saveRDS(sce, paste0(file_path_root, "data2/SC3.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_sc3.rds"))

#scCCESS-kmeans
library(scCCESS)
library(SingleCellExperiment)
library(reticulate)
library(tensorflow)
#use_virtualenv("r-reticulate")
data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
dat <- SingleCellExperiment(assays = list(counts = data))
dat <- counts(dat)

dat=prefilter(dat)
k = estimate_k(dat,
               seed = 1, 
               cluster_func = function(x,centers) { 
                 set.seed(42);
                 kmeans(x, centers)
               },
               criteria_method = "NMI",
               krange = 60:100, ensemble_sizes = 10,
               cores = 10
)
cluster = ensemble_cluster(dat, 
                           seed = 1, 
                           cluster_func = function(x) {
                             set.seed(1)
                             kmeans(x, centers = k$ngroups)
                           }, 
                           cores = 10, 
                           genes_as_rows = T, 
                           ensemble_sizes = 10, 
                           verbose = 0, 
                           scale = F, 
                           batch_size = 64
)
saveRDS(cluster,paste0(file_path_root, "data1/scCCESS_kmeans.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_kmeans.rds"))

data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
dat <- SingleCellExperiment(assays = list(counts = data))
dat <- counts(dat)

dat=prefilter(dat)
k = estimate_k(dat,
               seed = 1, 
               cluster_func = function(x,centers) { 
                 set.seed(42);
                 kmeans(x, centers)
               },
               criteria_method = "NMI",
               krange = 60:100, ensemble_sizes = 10,
               cores = 10
)

cluster = ensemble_cluster(dat, 
                           seed = 1, 
                           cluster_func = function(x) {
                             set.seed(1)
                             kmeans(x, centers = k$ngroups)
                           }, 
                           cores = 10, 
                           genes_as_rows = T, 
                           ensemble_sizes = 10, 
                           verbose = 0, 
                           scale = F, 
                           batch_size = 64
)
saveRDS(cluster,paste0(file_path_root, "data2/scCCESS_kmeans.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_kmeans.rds"))


#scCCESS-SIMLR
library(scCCESS)
library(SIMLR)

data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
dat <- SingleCellExperiment(assays = list(counts = data))
dat <- counts(dat)

k = estimate_k(dat,
               seed = 1, 
               cluster_func = function(x,centers) {
                 set.seed(42);
                 SIMLR_Large_Scale(t(x), c=centers,kk=15)
               },
               criteria_method = "NMI",
               krange = 60:100, ensemble_sizes = 10,
               cores = 10
)


cluster = ensemble_cluster(dat, 
                           seed = 1, 
                           cluster_func = function(x) {
                             set.seed(1)
                             SIMLR_Large_Scale(t(x), c=k$ngroups,kk=15)
                           }, 
                           cores = 10, 
                           genes_as_rows = T, 
                           ensemble_sizes = 10, 
                           verbose = 0, 
                           scale = F, 
                           batch_size = 64
)
saveRDS(cluster,paste0(file_path_root, "data1/scCCESS_simlr.rds"))
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_simlr.rds"))


data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
dat <- SingleCellExperiment(assays = list(counts = data))
dat <- counts(dat)

k = estimate_k(dat, 
               seed = 1, 
               cluster_func = function(x,centers) {
                 set.seed(42);
                 SIMLR_Large_Scale(t(x), c=centers,kk=15)
               },
               criteria_method = "NMI",
               krange = 60:100, ensemble_sizes = 10,
               cores = 10
)


cluster = ensemble_cluster(dat, 
                           seed = 1, 
                           cluster_func = function(x) {
                             set.seed(1)
                             SIMLR_Large_Scale(t(x), c=k$ngroups,kk=15)
                           }, 
                           cores = 10, 
                           genes_as_rows = T, 
                           ensemble_sizes = 10, 
                           verbose = 0, 
                           scale = F, 
                           batch_size = 64
)
saveRDS(cluster,paste0(file_path_root, "data2/scCCESS_simlr.rds"))

end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_simlr.rds"))








#consensus 100
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



data = read.csv(paste0(file_path_root,"data1/","data.csv"),row.names=1)
start <- Sys.time()
#data = data1_new
# data = aggregate(data[,-1], by = list(data[,1]), mean)s
# rownames(data) = data[,1]
# data = data[,-1]
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"run_cluster2"))
result_1 <- run_consensus_clust(norm.dat,
                                niter = 100,
                                de.param = de.param,
                                dim.method = "pca",
                                output_dir = paste0(file_path_root,"run_cluster2"),
                                mc.cores = 10)



group_meta = data.frame(names(result_1$cl.result$cl),result_1$cl.result$cl)
colnames(group_meta) = c("cell","group")

saveRDS(result_1, paste0(file_path_root,"data1/consesus_100_result_1.rds")) 
saveRDS(group_meta, paste0(file_path_root,"data1/group_meta_100.rds")) 
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root, "data1/","running_time_consensus100.rds"))



data = read.csv(paste0(file_path_root,"data2/","data.csv"),row.names=1)
start <- Sys.time()
#data = data1_new
# data = aggregate(data[,-1], by = list(data[,1]), mean)s
# rownames(data) = data[,1]
# data = data[,-1]
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"run_cluster2"))
result_1 <- run_consensus_clust(norm.dat,
                                niter = 100,
                                de.param = de.param,
                                dim.method = "pca",
                                output_dir = paste0(file_path_root,"run_cluster2"),
                                mc.cores = 10)



group_meta = data.frame(names(result_1$cl.result$cl),result_1$cl.result$cl)
colnames(group_meta) = c("cell","group")

saveRDS(result_1, paste0(file_path_root,"data2/consesus_100_result_2.rds")) 
saveRDS(group_meta, paste0(file_path_root,"data2/group_meta_100.rds")) 
end <- Sys.time()
runningtime <- end-start
>>>>>>> 3e2f6a6 (add large file)
saveRDS(runningtime, paste0(file_path_root, "data2/","running_time_consensus100.rds"))