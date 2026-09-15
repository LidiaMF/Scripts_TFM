#===============================================================================
# # TFM - Análisis integrado de la expresión génica y de proteínas de superficie 
# a nivel unicelular en muestras de médula ósea y sangre periférica de pacientes 
# con LLA-B (GSE230295)
#
#
# INTEGRACIÓN Y CLUSTERING:
# Normalización, selección de genes variables, escalado y PCA por librería;
# integración de las 6 librerías (RPCA) para corregir el efecto de lote; clustering
# y proyección UMAP.
#
#===============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(here)

#===============================================================================

#Cargamos datos

seu_list_filtered <- readRDS(here("2_seu_list_filtered_demux.rds"))

#Normalizamos los datos --> método "LogNormalize", normaliza las medidas de expresión 
#de cada células por el total de la expresion y lo multiplica x10000 y transforma el resultado mediante log.

seu_list_filtered <- lapply(seu_list_filtered, NormalizeData)


#Identificamos las variables con mayor variabilidad

seu_list_filtered <- lapply(seu_list_filtered, function(obj) {
  FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
})

#Representamos las 10 con mayor variabilidad

for (nm in names(seu_list_filtered)) {
  obj <- seu_list_filtered[[nm]]
  top10 <- head(VariableFeatures(obj), 10)
  
  
  # Comprobamos que no haya ningún objeto sin genes con variabilidad
  if (length(top10) == 0) {
    message(paste("No hay genes variables en", nm))
    next
  }
  
  
  # Plot base sin "log" para evitar el aviso de valores infinitos
  p_base <- VariableFeaturePlot(obj, log = FALSE) + 
    ggtitle(paste0("HVF - ", nm))
  
  
  hvf_top10 <- LabelPoints(
    plot = p_base,
    points = top10,
    repel = TRUE,
    xnudge = 0,
    ynudge = 0
  )
  
  # Guardamos 
  ggsave(
    filename = here("plots1", paste0("HVF_", nm, ".png")),
    plot = hvf_top10,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# Escalamos los datos de cada librería 

for (nm in names(seu_list_filtered)) {
  
  all.genes <- rownames(seu_list_filtered[[nm]])
  
  seu_list_filtered[[nm]] <- ScaleData(
    seu_list_filtered[[nm]],
    features = all.genes
  )
}

# Reduccción linear de la dimensionalidad (PCA)

for (nm in names(seu_list_filtered)) {
  obj <- seu_list_filtered[[nm]]
  obj <- RunPCA(
    seu_list_filtered[[nm]],
    features = VariableFeatures(object =  seu_list_filtered[[nm]]),
    npcs = 30,
    verbose = FALSE
  )
  
  seu_list_filtered[[nm]] <- obj
}




#Representamos y visualizamos el PCA 

print(seu_list_filtered[["L5"]][["pca"]], dims = 1:5, nfeatures = 10)

print(seu_list_filtered[["L6"]][["pca"]], dims = 1:5, nfeatures = 10)

print(seu_list_filtered[["L7"]][["pca"]], dims = 1:5, nfeatures = 10)

print(seu_list_filtered[["L8"]][["pca"]], dims = 1:5, nfeatures = 10)

print(seu_list_filtered[["L9"]][["pca"]], dims = 1:5, nfeatures = 10)

print(seu_list_filtered[["L10"]][["pca"]], dims = 1:5, nfeatures = 10)



# Determinar dimensionalidad de los datos con ElbowPlot

for (nm in names(seu_list_filtered)) {
  
  obj <- seu_list_filtered[[nm]]
  
  Ep <- ElbowPlot(obj, ndims = 30) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    ) +
    ggtitle(paste0("ElbowPlot - ", nm))  # Generamos el ElbowPlot de cada librería
  
  
  print(Ep)
  
  
  ggsave(
    filename = here("plots1", paste0("Elbow_", nm, ".png")),
    plot = Ep,
    width = 7,
    height = 5,
    dpi = 300
  )
} #las primeras 15 dims capturan gran parte de la varianza en todos los lotes, continuamos con 15dims


# Usamos VizDimLoadings y DimHeatmap para revisar los componentes principales


for (nm in names(seu_list_filtered)) {
  obj <- seu_list_filtered[[nm]]
  
  VizDimplot <- VizDimLoadings(obj,dims = 1:2, reduction = "pca") +
    ggtitle(paste0("PCA - ", nm))
  
  print(VizDimplot)
  
  ggsave(
    filename = here("plots1", paste0("VizDim_", nm, ".png")),
    plot = VizDimplot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  png(
    filename = here("plots1", paste0("Heatmap_", nm, "_PC1-PC5.png")),
    width = 1200,
    height = 1000,
    res = 150
  )
  
  DimHeatmap(
    obj,
    dims = 1:5,
    cells = 200,
    balanced = TRUE
  )
  
  dev.off()
  
  
}    


##INTEGRACIÓN LOTES## 
# Seleccionamos genes comunes para la integración
features <- SelectIntegrationFeatures(
  object.list = seu_list_filtered,
  nfeatures = 3000
)

# Comprobamos que se han seleccionado genes
length(features)
head(features)


# Recalculamos ScaleData y PCA usando esos genes comunes

seu_list_filtered <- lapply(seu_list_filtered, function(obj) {
  
  obj <- ScaleData(
    obj,
    features = features,
    verbose = FALSE
  )
  
  obj <- RunPCA(
    obj,
    features = features,
    npcs = 30,
    verbose = FALSE
  )
  
  return(obj)
})



#Buscamos anchors de integración usando RPCA
anchors_rpca <- FindIntegrationAnchors(
  object.list = seu_list_filtered,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:15
)


# Integramos las muestras y las escalamos 
dataset_integrado <- IntegrateData(anchorset = anchors_rpca,dims = 1:15)

DefaultAssay(dataset_integrado) <- "integrated" #Indicamos que dataset_integrado sean los datos corregidos (assay "integrated") y no los de RNA originales.
dataset_integrado <- ScaleData(dataset_integrado, verbose = FALSE)

#Volvemos a hacer PCA sobre el dataset integrado

dataset_integrado <- RunPCA(dataset_integrado, npcs = 30, verbose = FALSE)

#Buscamos vecinos para después hacer el clustering

dataset_integrado <- FindNeighbors(dataset_integrado, dims = 1:15)


#SELECCIÓN DE LA RESOLUCIÓN: probamos varias resoluciones y visualizamos con
#clustree cómo se van dividiendo los clusters al subir la resolución.
if (!requireNamespace("clustree", quietly = TRUE)) install.packages("clustree")
library(clustree)

resoluciones_a_probar <- seq(0.2, 1.2, by = 0.2)

dataset_integrado <- FindClusters(dataset_integrado, resolution = resoluciones_a_probar)

p_clustree <- clustree(dataset_integrado, prefix = "integrated_snn_res.")
p_clustree

ggsave(
  filename = here("plots1", "clustree_resoluciones.png"),
  plot = p_clustree,
  width = 12,
  height = 14,
  dpi = 300
)


#Tabla con nº de clusters resultante a cada resolución
n_clusters_por_resolucion <- sapply(resoluciones_a_probar, function(r) {
  col <- paste0("integrated_snn_res.", r)
  length(unique(dataset_integrado@meta.data[[col]]))
})

tabla_resoluciones <- data.frame(resolucion = resoluciones_a_probar, n_clusters = n_clusters_por_resolucion)
tabla_resoluciones

write.csv(tabla_resoluciones, file = here("tabla_resoluciones_clustering.csv"), row.names = FALSE)

#Clustering (con resolución 0.6)
dataset_integrado <- FindClusters(dataset_integrado, resolution = 0.6)


#Comprobación: esto debe coincidir con el nº de clusters que se obtienen con la resolución elegida(20)
length(unique(dataset_integrado$seurat_clusters))


#Usamos UMAP para visualizar el clustering
dataset_integrado <- RunUMAP(dataset_integrado, dims = 1:15)

#UMAP por clusters
Plot_UMAP <- DimPlot(dataset_integrado, reduction = "umap", label = TRUE) + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("UMAP-CLUSTERS")
Plot_UMAP

ggsave( filename = here("plots1", "UMAP_clusters.png"), plot = Plot_UMAP,  width = 8,  height = 6,  dpi = 300)

#UMAP por librería
UMAPlot_library <- DimPlot(dataset_integrado,  reduction = "umap",  group.by = "library") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("UMAP-by library")
UMAPlot_library

ggsave( filename = here("plots1", "UMAP_library.png"), plot = UMAPlot_library,  width = 8,  height = 6,  dpi = 300)

#UMAP por tejido (MO vs SP)
UMAPlot_tissue <- DimPlot(dataset_integrado,  reduction = "umap",  group.by = "tissue") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("UMAP - tissue")
UMAPlot_tissue
ggsave( filename = here("plots1", "UMAP_tissue.png"),plot = UMAPlot_tissue,width = 8,height = 6,dpi = 300)

#UMAP por momento clínico (Dx, 2d y 15d)
UMAPlot_timepoint <- DimPlot(dataset_integrado,  reduction = "umap",  group.by = "timepoint") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("UMAP - momento clínico")
UMAPlot_timepoint

ggsave( filename = here("plots1", "UMAP_timepoint.png"),plot = UMAPlot_timepoint,width = 8,height = 6,dpi = 300)


#UMAP por momento clínico pero solo con las células que sí tienen momento clínico asignado (sin NA)

dataset_con_timepoint <- subset(dataset_integrado, subset = !is.na(timepoint))

UMAPlot_timepoint_sin_NA <- DimPlot(dataset_con_timepoint, reduction = "umap", group.by = "timepoint") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("UMAP - momento clínico (solo células asignadas)")
UMAPlot_timepoint_sin_NA

ggsave( filename = here("plots1", "UMAP_timepoint_sin_NA.png"), plot = UMAPlot_timepoint_sin_NA, width = 8, height = 6, dpi = 300)

#Guardamos el objeto final tras la integración
saveRDS(dataset_integrado, file = here("3_dataset_integrado.rds"))

#Comprobamos cuantas células hay por cluster y las guardamos
table(dataset_integrado$seurat_clusters)

table(dataset_integrado$seurat_clusters, dataset_integrado$tissue)


tabla_cluster_tissue <- as.data.frame(table(dataset_integrado$seurat_clusters, dataset_integrado$tissue))

colnames(tabla_cluster_tissue) <- c("cluster", "tissue", "cell_number")

write.csv(tabla_cluster_tissue, file = here("cluster_by_tissue.csv"), row.names = FALSE)


#Calculamos la proporcion de células de cada tejido que hay en cada cluster (margin=1) para comprobar que la integración ha ido bien y se ha corregido el "batch effect"
cluster_percent_by_cluster <- prop.table(table(dataset_integrado$seurat_clusters, dataset_integrado$tissue),margin = 1) * 100
cluster_percent_by_cluster

cluster_percent_by_cluster_df <- as.data.frame(cluster_percent_by_cluster)
colnames(cluster_percent_by_cluster_df) <- c("cluster","tissue","percentage")


write.csv(cluster_percent_by_cluster_df, file = here("cluster_percent_by_cluster.csv"), row.names = FALSE)

#Calculamos la proporcion de células que hay en cada cluster dentro de cadaa tejido (margin=2)

cluster_percent_by_tissue <- prop.table(table(dataset_integrado$seurat_clusters, dataset_integrado$tissue),margin = 2) * 100
cluster_percent_by_tissue

cluster_percent_by_tissue_df <- as.data.frame(cluster_percent_by_tissue)
colnames(cluster_percent_by_tissue_df) <- c("cluster","tissue","percentage")


write.csv(cluster_percent_by_tissue_df, file = here("cluster_percent_by_tissue.csv"),row.names = FALSE)

#Calculamos la proporción de células en cada cluster dentro de cada librería 

table(dataset_integrado$seurat_clusters, dataset_integrado$library)

cluster_percent_by_library <- prop.table(table(dataset_integrado$seurat_clusters, dataset_integrado$library), margin = 1) * 100
cluster_percent_by_library

cluster_percent_by_library_df <- as.data.frame(cluster_percent_by_library)
colnames(cluster_percent_by_library_df) <- c("cluster","library","percentage")

write.csv(cluster_percent_by_library_df, file = here("cluster_percent_by_library.csv"), row.names = FALSE)

#Información de la sesión
sessionInfo()

