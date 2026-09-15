#===============================================================================
# # TFM - Análisis integrado de la expresión génica y de proteínas de superficie 
# a nivel unicelular en muestras de médula ósea y sangre periférica de pacientes 
# con LLA-B (GSE230295)
#
# Anotacion celular automática con SingleR, usando dos referencias de celldex:
# NovershternHematopoieticData y HumanPrimaryCellAtlasData . 
#
# 
#===============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(here)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("SingleR", quietly = TRUE)) BiocManager::install("SingleR")
if (!requireNamespace("celldex", quietly = TRUE)) BiocManager::install("celldex")

library(SingleR)
library(celldex)

#===============================================================================


# Cargamos los objetos seurat 
dataset_integrado <- readRDS(here("4_dataset_integrado_anotado.rds"))
#Necesitaremos los datos normalizados, el assay "RNA"
DefaultAssay(dataset_integrado) <- "RNA"


#===============================================================================
# REFERENCIA 1: NovershternHematopoieticData 
#===============================================================================
# Cargamos los datos del dataset de referencia (NovershternHematopoieticData )
ref_novershtern <- NovershternHematopoieticData()

Result_singler_novershtern <- SingleR(
  test = GetAssayData(dataset_integrado, assay = "RNA", layer = "data"),#le decimos que coja la matriz de expresión normalizada dentro del objeto seurat
  ref = ref_novershtern,
  labels = ref_novershtern$label.main #Con esto obtenmos una etiqueta/célula de categorías celulares principales 
)


#Añadimos las etiquetas que predice SingleR a los metadatos del objeto seurat
dataset_integrado$singler_novershtern <- Result_singler_novershtern$labels

#Añadimos por otra parte las etiquetas "podadas" para evitar asignaciones ambiguas(se muestras como NA)
dataset_integrado$singler_novershtern_pruned <- Result_singler_novershtern$pruned.labels

#Células con etiqueta "podada"
table(is.na(dataset_integrado$singler_novershtern_pruned))

#% Células ambiguas
pct_pruned_novershtern <- round(mean(is.na(dataset_integrado$singler_novershtern_pruned)) * 100, 2)

pct_pruned_novershtern


#% Células ambiguas por cluster
pruned_por_cluster_novershtern <-
  dataset_integrado@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_celulas = n(),
    n_pruned = sum(
      is.na(singler_novershtern_pruned)
      ),
    pct_pruned = round(
      100 * n_pruned / n_celulas,
      2
      )
    
    ) %>%
  arrange(desc(pct_pruned))
pruned_por_cluster_novershtern 


#Guardamos los datos
write.csv(pruned_por_cluster_novershtern, file = here("SingleR_Novershtern_pruned_summary.csv"), row.names = FALSE)



#Visualizamos los clusters anotados (NovershternHematopoieticData)  mediante SingleR
p_umap_novershtern <- DimPlot(
  dataset_integrado,
  reduction = "umap",
  group.by = "singler_novershtern"
) +
  ggtitle("UMAP anotado por SingleR (Novershtern)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  theme(legend.position = "right")
p_umap_novershtern

ggsave(
  filename = here("plots1", "UMAP_singleR_novershtern.png"),
  plot = p_umap_novershtern,
  width = 12,
  height = 8,
  dpi = 300
)

#Etiqueta de SingleR mayoritaria dentro de cada cluster
tabla_cluster_novershtern <- as.data.frame(table(dataset_integrado$seurat_clusters, dataset_integrado$singler_novershtern))
colnames(tabla_cluster_novershtern) <- c("cluster", "singler_label", "cell_number")

novershtern_por_cluster <- tabla_cluster_novershtern %>%
  group_by(cluster) %>%
  slice_max(order_by = cell_number, n = 1) %>%
  ungroup()
colnames(novershtern_por_cluster) <- c("cluster", "novershtern_label", "novershtern_n_cells")
novershtern_por_cluster$cluster <- as.character(novershtern_por_cluster$cluster)

#Tabla comparativa: anotación manual vs Novershtern (juntamos tablas)
tabla_markers <- read.csv(here("tabla_markers_by_cluster.csv"))
tabla_markers$cluster <- as.character(tabla_markers$cluster)

comparacion_novershtern <- merge(tabla_markers, novershtern_por_cluster, by = "cluster")
comparacion_novershtern$pct_confianza_novershtern <- round(
  comparacion_novershtern$novershtern_n_cells / comparacion_novershtern$n_celulas * 100, 1
)

comparacion_novershtern

#Categoría de concordancia al menos de linaje (Coincidente/No coincidente), decidida tras revisar caso por caso
#el linaje de la etiqueta de Novershtern frente a la anotación manual
categoria_novershtern <- c(
  "0" = "Coincidente", "1" = "No coincidente", "2" = "Coincidente", "3" = "No coincidente",
  "4" = "No coincidente", "5" = "Coincidente", "6" = "No coincidente", "7" = "Coincidente",
  "8" = "Coincidente", "9" = "Coincidente", "10" = "No coincidente", "11" = "Coincidente",
  "12" = "Coincidente", "13" = "Coincidente", "14" = "Coincidente", "15" = "Coincidente",
  "16" = "No coincidente", "17" = "No coincidente", "18" = "No coincidente", "19" = "Coincidente"
)
comparacion_novershtern$coincidencia <- unname(categoria_novershtern[comparacion_novershtern$cluster])
comparacion_novershtern <- comparacion_novershtern[order(as.numeric(comparacion_novershtern$cluster)), ]

comparacion_novershtern

write.csv(comparacion_novershtern, file = here("comparacion_singleR_novershtern.csv"), row.names = FALSE)


#===============================================================================
# REFERENCIA 2: HumanPrimaryCellAtlasData 
#===============================================================================

ref_hpca <- HumanPrimaryCellAtlasData()

Result_singler_hpca <- SingleR(
  test = GetAssayData(dataset_integrado, assay = "RNA", layer = "data"),
  ref = ref_hpca,
  labels = ref_hpca$label.main
)


#Añadimos las etiquetas que predice SingleR a los metadatos del objeto seurat
dataset_integrado$singler_hpca <- Result_singler_hpca$labels

#Añadimos por otra parte las etiquetas "podadas" para evitar asignaciones ambiguas(se muestras como NA)
dataset_integrado$singler_hpca_pruned <- Result_singler_hpca$pruned.labels

#Células con etiqueta "podada"
table(is.na(dataset_integrado$singler_hpca_pruned))

#% Células ambiguas
pct_pruned_hpca <- round(mean(is.na(dataset_integrado$singler_hpca_pruned)) * 100, 2)

pct_pruned_hpca


#% Células ambiguas por cluster
pruned_por_cluster_hpca <-
  dataset_integrado@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_celulas = n(),
    n_pruned = sum(
      is.na(singler_hpca_pruned)
    ),
    pct_pruned = round(
      100 * n_pruned / n_celulas,
      2
    )
    
  ) %>%
  arrange(desc(pct_pruned))
pruned_por_cluster_hpca 

#Guardamos los datos
write.csv(pruned_por_cluster_hpca, file = here("SingleR_HPCA_pruned_summary.csv"), row.names = FALSE)


#UMAP anotado por SingleR (HPCA)
p_umap_hpca <- DimPlot(
  dataset_integrado,
  reduction = "umap",
  group.by = "singler_hpca"
) +
  ggtitle("UMAP anotado por SingleR (HPCA)") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  theme(legend.position = "right")
p_umap_hpca

ggsave(
  filename = here("plots1", "UMAP_singleR_hpca.png"),
  plot = p_umap_hpca,
  width = 12,
  height = 8,
  dpi = 300
)

#Etiqueta de SingleR mayoritaria dentro de cada cluster
tabla_cluster_hpca <- as.data.frame(table(dataset_integrado$seurat_clusters, dataset_integrado$singler_hpca))
colnames(tabla_cluster_hpca) <- c("cluster", "singler_label", "cell_number")

hpca_por_cluster <- tabla_cluster_hpca %>%
  group_by(cluster) %>%
  slice_max(order_by = cell_number, n = 1) %>%
  ungroup()
colnames(hpca_por_cluster) <- c("cluster", "hpca_label", "hpca_n_cells")
hpca_por_cluster$cluster <- as.character(hpca_por_cluster$cluster)

#Tabla comparativa: anotación manual vs HPCA
comparacion_hpca <- merge(tabla_markers, hpca_por_cluster, by = "cluster")
comparacion_hpca$pct_confianza_hpca <- round(
  comparacion_hpca$hpca_n_cells / comparacion_hpca$n_celulas * 100, 1
)

comparacion_hpca

#Categoría de concordancia (Coincidente/No coincidente) para HPCA

categoria_hpca <- c(
  "0" = "Coincidente", "1" = "No coincidente", "2" = "Coincidente", "3" = "No coincidente",
  "4" = "Coincidente", "5" = "Coincidente", "6" = "Coincidente", "7" = "Coincidente",
  "8" = "No coincidente", "9" = "Coincidente", "10" = "No coincidente", "11" = "Coincidente",
  "12" = "Coincidente", "13" = "Coincidente", "14" = "Coincidente", "15" = "Coincidente",
  "16" = "No coincidente", "17" = "Coincidente", "18" = "No coincidente", "19" = "Coincidente"
)
comparacion_hpca$coincidencia <- unname(categoria_hpca[comparacion_hpca$cluster])
comparacion_hpca <- comparacion_hpca[order(as.numeric(comparacion_hpca$cluster)), ]

comparacion_hpca

write.csv(comparacion_hpca, file = here("comparacion_singleR_hpca.csv"), row.names = FALSE)

#Tabla comparativa anotación manual vs. SingleR
tabla_anotacion_manual_singleR <- comparacion_novershtern %>%
  select(
    cluster,
    id_asignada,
    novershtern_label,
  ) %>%
  left_join(
    comparacion_hpca %>%
      select(
        cluster,
        hpca_label
      ),
    by = "cluster"
  )

tabla_anotacion_manual_singleR 

write.csv(tabla_anotacion_manual_singleR, file = here("tabla_anotacion_manual_singleR.csv"), row.names = FALSE)

#Volvemos a dejar RNA como assay por defecto
DefaultAssay(dataset_integrado) <- "RNA"

#Guardamos los datos: dataset_integrado con las dos etiquetas de SingleR añadidas
saveRDS(dataset_integrado, file = here("5_dataset_integrado_validado.rds"))

#Información de la sesión
sessionInfo()