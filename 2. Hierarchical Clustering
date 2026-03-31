library(factoextra)

# 1.- Select genes for clustering -----------------------------------------

# 1.1 Filter more variable genes

# 1.1.1 Pick 5000 more variable genes

top.genes  <- order(microarray_data.unique$variance, decreasing = TRUE)[1:500] 

# 1.1.2 Create df with only 5000 genes

microarray_data.fil <- microarray_data.unique[top.genes, ] 

# 1.2 Delete variance columns

microarray_data.unique$variance <- NULL

microarray_data.fil$variance <- NULL

# 1.3 Transpose so that the rows correspond to patients and the columns to genes

microarray_data.org <- t(microarray_data.fil)

# 1.4 Normalize

microarray_data.norm <- scale(microarray_data.org)


# 2.- Clustering ----------------------------------------------------------


# 2.1 Select clustering number

fviz_nbclust(microarray_data.norm, FUN = hcut, method = "silhouette")

# 2.2 Distance matrix between samples

dist_microarray.brca <- get_dist(microarray_data.norm, method = "euclidean")

# 2.3 Clustering jerárquico

hc.out_microarray.brca <- hclust(dist_microarray.brca, method = "ward.D2")


# 3.- Clustering characteristics ------------------------------------------


# 3.1 Dendrogram

plot(
  hc.out_microarray.brca,
  labels = FALSE,
  hang = -1,
  main = "CLUSTERING JERÁRQUICO METABRIC"
)

rect.hclust(hc.out_microarray.brca,
            k = 2,
            border = "purple") # bottom up approach

# 3.2 Tree of clusters

clusters <- cutree(hc.out_microarray.brca, k = 2)

# 3.2.2 Observe how many patients in each cluster

table(clusters)

# 3.3 Create object as data frame

cluster_df <- data.frame(CLUSTER = clusters)

# 3.4 Asign a column named PATIENT ID with the rownames such that we can then merge based on that column

cluster_df$"PATIENT_ID" <- rownames(microarray_data.norm)

# 3.5 Merge object so that in the metadata there is a column corresponding to that patients cluster

metadata <- merge(
  cluster_df,
  metadata,  
  by = "PATIENT_ID"
)
