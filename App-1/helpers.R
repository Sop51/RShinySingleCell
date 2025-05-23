library(Seurat)
library(ggplot2)

# create a function to generate a feature plot for a certain gene
generate_feature_plot <- function(dataset, gene){
  # code to create the plot
  plot <- FeaturePlot(object = dataset, 
                   features = c(gene),
                   order = TRUE,
                   min.cutoff = 'q10', 
                   repel = TRUE)
  return(plot)
}

# create a function to generate a UMAP with the selected metric
generate_umap_plot <- function(dataset, selected_metric){
  # code to create the plot
  plot <- DimPlot(object = dataset,
                  reduction = "umap", 
                  group.by = selected_metric,
                  repel = TRUE,
                  label = TRUE)
  return(plot)
}

# create a function to return a DE table based on selected cell type
generate_de_table <- function(dataset, cell_type, cell_type_meta) {
  de_table <- FindMarkers(dataset,
                          slot = 'data',
                          ident.1 = cell_type,
                          group.by = cell_type_meta)
  de_table <- as.data.frame(de_table)
  return(de_table)
}

# create a helper function to create a violin plot from a meta var
generate_violin_plot <- function(dataset, gene, meta_var){
    plot <- VlnPlot(dataset, features = gene, group.by = meta_var)
    return(plot)
}

# generate a umap of a specific cell type
generate_subcluster_umap <- function(dataset, meta_var){
  plot <- DimPlot(object = dataset,
                  reduction = "umap", 
                  group.by = meta_var)
  return(plot)
}

# generate subcluster feature plot
generate_subcluster_featureplot <- function(dataset, gene){
  # code to create the plot
  plot <- FeaturePlot(object = dataset, 
                      features = c(gene),
                      order = TRUE,
                      min.cutoff = 'q10', 
                      repel = TRUE)
  return(plot)
}

# function to create a cell type count summary table
cell_type_count <- function(dataset, cell_type_col){
  cell_type_counts <- table(dataset[[cell_type_col]])
  cell_type_counts_df <- as.data.frame(cell_type_counts)
  colnames(cell_type_counts_df) <- c("Cell Type", "Cell Count")
  return(cell_type_counts_df)
}

# function to create plot of gene across cell types
gene_across_cell_type <- function(dataset, gene, cell_col_meta) {
  plot <- VlnPlot(dataset, features = gene, group.by = cell_col_meta)
  return(plot)
}