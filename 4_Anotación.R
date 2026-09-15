#===============================================================================
# # TFM - Análisis integrado de la expresión génica y de proteínas de superficie 
# a nivel unicelular en muestras de médula ósea y sangre periférica de pacientes 
# con LLA-B (GSE230295)
#
# ANOTACIÓN: 
# Identificación de genes marcadores por clúster (FindAllMarkers), construcción de
# un panel de marcadores de ARN y combinación con la proteína de superficie
# (ADT/CITE-seq) para anotar manualmente la identidad de cada clúster. 
#
#===============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(here)
#===============================================================================

# Cargamos los objetos seurat 
dataset_integrado <- readRDS(here("3_dataset_integrado.rds"))
adt_list <- readRDS(here("1_adt_list.rds"))

#BUSCAR GENES QUE DEFINEN CADA CLUSTER
#Ahora necesitamos trabajar con el assay "RNA", es decir, con los valores de expresión brutos
#Juntamos los datos de expresión de las distintas librerías que están en distintas capas en una.

dataset_integrado <- JoinLayers(dataset_integrado,assay = "RNA")
DefaultAssay(dataset_integrado) <- "RNA"

Idents(dataset_integrado) <- "seurat_clusters"

#Identificar los marcadores que se expresan en cada cluster con FindAllMarkers()
cluster_markers <- FindAllMarkers(dataset_integrado,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25, verbose = FALSE)
head(cluster_markers) 

write.csv(cluster_markers,file = here("cluster_markers.csv"),row.names = FALSE)

#Identificar los 10 genes que están más expresados en cada cluster, con un avg_log2FC mayor de 1 y significativos (p_val_adj < 0.05)

top10_markers <- cluster_markers %>%
  filter(avg_log2FC > 1, p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

top10_markers

write.csv(top10_markers, file = here("top10_markers_avglog2FC_mayor_1.csv"),row.names = FALSE)

  #Comprobamos los markers de los tres cluster sin identidad
  top10_sin_identidad <- top10_markers %>%
    filter(cluster %in% c("1", "3", "10"))
  
  top10_sin_identidad
  
  write.csv(top10_sin_identidad, file = here("top10_markers_clusters_sin_identidad.csv"), row.names = FALSE)
  


#Heatmap de los 10 genes más expresados
 #Al representar el Heatmap sale un aviso indicando que algunos de los 10 genes más expresados no están en la capa scale.data del assay RNA
 #Esta capa sólo tiene los 300 gnes de integración.Por ello, volvemos a escalar los genes que se van a representar en el heatmap para evitar
 #que se omitan. Luego representamos el heatmap
dataset_integrado <- ScaleData(dataset_integrado, features = top10_markers$gene, assay = "RNA")

heatmap_top10 <- DoHeatmap(dataset_integrado, features = top10_markers$gene, label = TRUE) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("Heatmap de los 10 genes más upregulados por Cluster")

heatmap_top10

ggsave(
  filename = here("plots1", "Heatmap_top10_genes_por_cluster.png"),
  plot = heatmap_top10,
  width = 14,
  height = 10,
  dpi = 300
)



#ANOTACION DE GENES
#Definimos poblaciones en función de marcadores canónicos ya descritos en la literatura
marker_genes <- c(
  # Linfocitos B y precursoras
  "CD79A", "CD79B", "MS4A1", "CD19", "VPREB1", "IGLL1", "MME",
  
  # Progenitores multipotentes / Células troncales hematopoyéticasl (HSC);
  # NPM1: marcador general de HSPC, Belderbos et al. 2026, Nat Immunol)
  "CD34", "KIT", "PROM1", "AVP", "HLF", "MECOM", "NPM1",
  
  # Progenitor linfoide temprano / pro-B vs pre-B 
  #DNTT+ marca LP/pro-B; MS4A1 marca la transición a pre-B-II/B inmadura;
  #EBF1 marca el subtipo LyP comprometido a linaje B (LyP-B); 
  #SPINK2 marca el subtipo LyP más multipotente/"stemness" (LyP-S); 
  #Aunque inicialmente habíamos añadido CD24 que marca LMPP, 
  #este no está incluido en la matriz de ARN, si está el marcador de proteína (ADT))
  "DNTT", "RAG1", "RAG2", "FLT3", "EBF1", "SPINK2", 
  
  # Mieloide / monocitos
  "LYZ", "CD14", "LST1", "S100A8", "S100A9", "FCGR3A", "FCN1",
  
  # Granulocítico (precursores inmaduros, característicos de médula ósea)
  "MPO", "ELANE", "AZU1", "PRTN3", "CSF3R",
  
  # Neutrófilos maduros (típico de sangre periférica; pierden MPO/ELANE/AZU1/PRTN3 al madurar,

  "FCGR3B", "CXCR2",
  
  # Linfocitos T (#TRAC no está en la matriz de ARN por eso no está incluido)
  "CD3D", "CD3E", "CD2", "IL7R",
  
  # Linfocitos NK / citotóxicas
  "NKG7", "GNLY", "GZMA", "GZMB", "KLRD1", "PRF1",
  
  # Eritroide
  "HBB", "HBA1", "HBA2", "ALAS2", "GYPA",
  
  # Megacariocitos/ plaquetas
  "PPBP", "PF4", "ITGA2B", "GP9",
  
  # Genes asociados a proliferación
  "MKI67", "TOP2A", "STMN1", "PCNA"
)

dotplot_marker_genes<- DotPlot(
  dataset_integrado,
  features = marker_genes
) +
  RotatedAxis() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("Genes marcadores por clúster")

dotplot_marker_genes

ggsave(
  filename = here("plots1", "DotPlot_marker_genes_clusters.png"),
  plot = dotplot_marker_genes,
  width = 16,
  height = 8,
  dpi = 300
)



#Para mejorar la anotación de tipos celulares, utilizamos la informacoin de proteína (ADT)
#Juntamos el ADT de las librerías en una matriz (adt_list)
adt_L5 <- adt_list$L5
adt_L6 <- adt_list$L6
adt_L7 <- adt_list$L7
adt_L8 <- adt_list$L8
adt_L9 <- adt_list$L9
adt_L10 <- adt_list$L10

adt_L5_L10 <- cbind(adt_L5, adt_L6, adt_L7, adt_L8, adt_L9, adt_L10)

#La utilizamos con el dataset integrado
adt_L5_L10 <- adt_L5_L10[, colnames(dataset_integrado)]

dataset_integrado[["ADT"]] <- CreateAssayObject(counts = adt_L5_L10)

#Utilizamos el método CLR (centered log-ratio) como forma de normalización estándar para proteína, distinta de la de ARN
dataset_integrado <- NormalizeData(dataset_integrado, assay = "ADT", normalization.method = "CLR")

DefaultAssay(dataset_integrado) <- "ADT"


dotplot_adt <- DotPlot(
  dataset_integrado,
  features = rownames(dataset_integrado[["ADT"]])
) +
  RotatedAxis() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("Marcadores de proteína (ADT) por clúster")

dotplot_adt

ggsave(
  filename = here("plots1", "DotPlot_ADT_clusters.png"),
  plot = dotplot_adt,
  width = 14,
  height = 8,
  dpi = 300
)


# Sacamos tabla con los datos del DotPlot de ARN
dotplot_data_rna <- dotplot_marker_genes$data
write.csv(dotplot_data_rna, file = here("dotplot_data_rna.csv"), row.names = FALSE)

# Lo mismo para el de proteína (ADT)
dotplot_data_adt <- dotplot_adt$data
write.csv(dotplot_data_adt, file = here("dotplot_data_adt.csv"), row.names = FALSE)

#Volvemos a dejar RNA como assay por defecto, por si se sigue trabajando con el script después de esto
DefaultAssay(dataset_integrado) <- "RNA"

#Generamos una tabla que combine los genes y/o proteínas más expresados en cada cluster
# Cargamos las dos tablas de datos del DotPlot
dotplot_rna <- read.csv(here("dotplot_data_rna.csv"))
dotplot_adt <- read.csv(here("dotplot_data_adt.csv"))


#Marcadores ARN por cluster: top 6 genes con pct.exp > 20% y mayor avg.exp.scaled 
# "id" es Factor (no texto), lo convertimos a character 
dotplot_rna$id <- as.character(dotplot_rna$id)

resultado_rna <- dotplot_rna %>%
  filter(pct.exp > 20) %>%
  group_by(id) %>%
  slice_max(order_by = avg.exp.scaled, n = 6) %>%
  summarise(
    texto = paste0(features.plot, " (", round(pct.exp), "%/", round(avg.exp.scaled, 2), ")", collapse = "; "),
    .groups = "drop"
  )

top_rna <- data.frame(
  cluster = resultado_rna$id,
  marcadores_ARN = resultado_rna$texto,
  stringsAsFactors = FALSE
)


# Marcadores ADT (proteína) por cluster: top 6 con pct.exp > 50% y mayor avg.exp.scaled

dotplot_adt$id <- as.character(dotplot_adt$id)

resultado_adt <- dotplot_adt %>%
  filter(pct.exp > 50) %>%
  group_by(id) %>%
  slice_max(order_by = avg.exp.scaled, n = 6) %>%
  summarise(
    texto = paste0(features.plot, " (", round(pct.exp), "%/", round(avg.exp.scaled, 2), ")", collapse = "; "),
    .groups = "drop"
  )

top_adt <- data.frame(
  cluster = resultado_adt$id,
  marcadores_ADT = resultado_adt$texto,
  stringsAsFactors = FALSE
)

colnames(top_rna)
colnames(top_adt)

# Nº de células por cluster
n_cells <- as.data.frame(table(dataset_integrado$seurat_clusters))
colnames(n_cells) <- c("cluster", "n_celulas")
n_cells$cluster <- as.character(n_cells$cluster)

# Tabla de resultado que muestra los marcadores que definen cada cluster 
tabla_markers_by_cluster <- n_cells %>%
  left_join(top_rna, by = "cluster") %>%
  left_join(top_adt, by = "cluster")

tabla_markers_by_cluster <- tabla_markers_by_cluster[order(as.numeric(tabla_markers_by_cluster$cluster)), ]

tabla_markers_by_cluster

write.csv(tabla_markers_by_cluster, file = here("tabla_markers_by_cluster.csv"), row.names = FALSE)


#Asignación de identidades celulares a cada cluster
cluster_labels <- c(
  "0"  = "Linfocitos T (CD4/CD8 mixto)",
  "1"  = "Sin identidad clara",
  "2"  = "Linfocitos B maduros",
  "3"  = "Sin identidad clara",
  "4"  = "Linfocitos NK citotóxicas",
  "5"  = "Linfocitos NK (CD16+)",
  "6"  = "Progenitor linfoide multipotente (LyP-S)",
  "7"  = "Progenitor B comprometido (LyP-B/pro-B)",
  "8"  = "LMPP",
  "9"  = "Progenitor B comprometido, proliferando",
  "10" = "Sin identidad clara",
  "11" = "Progenitor B muy proliferativo (G2/M)",
  "12" = "Linfocitos T CD8+",
  "13" = "Monocitos",
  "14" = "Monocitos/mieloide maduro",
  "15" = "Progenitor pro-B/pre-B en recombinación V(D)J",
  "16" = "Linaje eritroide",
  "17" = "Precursor granulocítico",
  "18" = "HSC",
  "19" = "Linfocitos T CD4+"
)


#Añadimos una columna en el objeto seurat con las etiquetas
#unname() es necesario: sin él, el vector resultante queda con los números de cluster como "nombres"
#(en vez de vacío/posicional) y aparece un error

dataset_integrado$cell_type <- unname(cluster_labels[as.character(dataset_integrado$seurat_clusters)])

#UMAP final con nombres de tipo celular 

celltype_colors <- c(
  "Linaje eritroide" = "#E41A1C",
  "HSC" = "#FFD700",
  "Linfocitos B maduros" = "#377EB8",
  "Linfocitos T (CD4/CD8 mixto)" = "#4DAF4A",
  "Linfocitos T CD4+" = "#A6D854",
  "Linfocitos T CD8+" = "#1B9E77",
  "LMPP" = "#FF7F00",
  "Monocitos" = "#984EA3",
  "Monocitos/mieloide maduro" = "#B15928",
  "Linfocitos NK (CD16+)" = "#66C2A5",
  "Linfocitos NK citotóxicas" = "#8DA0CB",
  "Precursor granulocítico" = "#FDBF6F",
  "Progenitor pro-B/pre-B en recombinación V(D)J" = "#E7298A",
  "Progenitor B comprometido (LyP-B/pro-B)" = "#7570B3",
  "Progenitor B comprometido, proliferando" = "#D95F02",
  "Progenitor B muy proliferativo (G2/M)" = "#E6AB02",
  "Progenitor linfoide multipotente (LyP-S)" = "#A65628",
  "Sin identidad clara" = "#999999"
)

umap_celltype <- DimPlot(
  dataset_integrado,
  reduction = "umap",
  group.by = "cell_type"
) +
  scale_color_manual(values = celltype_colors) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("UMAP anotado por tipo celular") +
  theme(legend.position = "right") #legenda a la derecha

umap_celltype

ggsave(
  filename = here("plots1", "UMAP_cell_type_final.png"),
  plot = umap_celltype,
  width = 12,
  height = 8,
  dpi = 300
)

#Añadimos las etiquetas en la tabla que habíamos creado
tabla_markers <- read.csv(here("tabla_markers_by_cluster.csv"))
tabla_markers$cluster <- as.character(tabla_markers$cluster)

tabla_markers$id_asignada <- cluster_labels[tabla_markers$cluster]

# La reordenamos para que quede: cluster, identidad, n_celulas, marcadores...
tabla_markers <- tabla_markers[, c("cluster", "id_asignada", "n_celulas", "marcadores_ARN", "marcadores_ADT")]

write.csv(tabla_markers, file = here("tabla_markers_by_cluster.csv"), row.names = FALSE)


#Guardamos datos: dataset_integrado ya con cell_type y el assay ADT añadidos
saveRDS(dataset_integrado, file = here("4_dataset_integrado_anotado.rds"))

#Información de la sesión
sessionInfo()
