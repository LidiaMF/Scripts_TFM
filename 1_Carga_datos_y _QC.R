#===============================================================================
# # TFM - Análisis integrado de la expresión génica y de proteínas de superficie 
# a nivel unicelular en muestras de médula ósea y sangre periférica de pacientes 
# con LLA-B (GSE230295)
#
# Se cargan 6 librerías (L5-L10) y se separa la información correspondiente a AR,
# proteína(ADT) y los hashtags(HTO), se lleva a cabo el cálculo de distintas métricas 
# de control de calidad, se filtran y eliminan células de baja calidad, y se detecta
# y eliminan dobletes..
#
#===============================================================================

#Cargamos los paquetes
library(Seurat)
library(ggplot2)
library(dplyr)
library(here)
library(scDblFinder)
library(SingleCellExperiment)
library(devtools)


here() #para que no haya rutas absolutas

# Creamos una carpeta para guardar las figuras
if (!dir.exists(here("plots1"))) {
  dir.create(here("plots1"))
}

#===============================================================================

###CARGAMOS LOS DATOS Y CREAMOS EL OBJETO SEURAT###


#Cargamos los datasets.Son archivos archivos .mtx.gz/.tsv.gz originales de Cell Ranger.
#Los archivos matrix.mtx son matrices de conteo.
#En los archivos features.tsv hay 3 tipos de filas:  "Gene Expression", "Antibody capture" en los que se indican los 5 hashtag utilizados 
#en el multiplexado y "Antibody Capture" que hacen referencia a los marcadores de superficie utilizado en el CITE-seq
#Dividimos estos tres tipos de filas de cada librería (L5-L10)


#Creamos un vectos con los nombres de las etiquetas usadas en el multiplexado
hashtags <- c("Hash1","Hash2","Hash3","Hash4","Hash5")

#L5
features_L5 <- read.delim(here("L5", "GSE230295_L5features.tsv.gz"), header = FALSE)
gene_rows_L5 <- features_L5$V3 == "Gene Expression"
adt_rows_L5 <- features_L5$V3 == "Antibody Capture" & !(features_L5$V2 %in% hashtags)
hto_rows_L5 <- features_L5$V2 %in% hashtags

data_10x_L5_raw <- ReadMtx(
  mtx = here("L5/GSE230295_L5matrix.mtx.gz"),
  features = here("L5/GSE230295_L5features.tsv.gz"),
  cells = here("L5/GSE230295_L5barcodes.tsv.gz"),
  feature.column = 2
)
adt_L5 <- data_10x_L5_raw[adt_rows_L5, ]
colnames(adt_L5) <- paste0("L5_", colnames(adt_L5)) #mismo prefijo que pondrá RenameCells más abajo
hto_L5 <- data_10x_L5_raw[hto_rows_L5, ]
colnames(hto_L5) <- paste0("L5_", colnames(hto_L5))
data_10x_L5 <- data_10x_L5_raw[gene_rows_L5, ]
data_10x_L5

#L6
features_L6 <- read.delim(here("L6", "GSE230295_L6features.tsv.gz"), header = FALSE)
gene_rows_L6 <- features_L6$V3 == "Gene Expression"
adt_rows_L6 <- features_L6$V3 == "Antibody Capture" & !(features_L6$V2 %in% hashtags)
hto_rows_L6 <- features_L6$V2 %in% hashtags

data_10x_L6_raw <- ReadMtx(
  mtx = here("L6/GSE230295_L6matrix.mtx.gz"),
  features = here("L6/GSE230295_L6features.tsv.gz"),
  cells = here("L6/GSE230295_L6barcodes.tsv.gz"),
  feature.column = 2
)
adt_L6 <- data_10x_L6_raw[adt_rows_L6, ]
colnames(adt_L6) <- paste0("L6_", colnames(adt_L6))
hto_L6 <- data_10x_L6_raw[hto_rows_L6, ]
colnames(hto_L6) <- paste0("L6_", colnames(hto_L6))
data_10x_L6 <- data_10x_L6_raw[gene_rows_L6, ]

#L7
features_L7 <- read.delim(here("L7", "GSE230295_L7features.tsv.gz"), header = FALSE)
gene_rows_L7 <- features_L7$V3 == "Gene Expression"
adt_rows_L7 <- features_L7$V3 == "Antibody Capture" & !(features_L7$V2 %in% hashtags)
hto_rows_L7 <- features_L7$V2 %in% hashtags

data_10x_L7_raw <- ReadMtx(
  mtx = here("L7/GSE230295_L7matrix.mtx.gz"),
  features = here("L7/GSE230295_L7features.tsv.gz"),
  cells = here("L7/GSE230295_L7barcodes.tsv.gz"),
  feature.column = 2
)
adt_L7 <- data_10x_L7_raw[adt_rows_L7, ]
colnames(adt_L7) <- paste0("L7_", colnames(adt_L7))
hto_L7 <- data_10x_L7_raw[hto_rows_L7, ]
colnames(hto_L7) <- paste0("L7_", colnames(hto_L7))
data_10x_L7 <- data_10x_L7_raw[gene_rows_L7, ]


#L8
features_L8 <- read.delim(here("L8", "GSE230295_L8features.tsv.gz"), header = FALSE)
gene_rows_L8 <- features_L8$V3 == "Gene Expression"
adt_rows_L8 <- features_L8$V3 == "Antibody Capture" & !(features_L8$V2 %in% hashtags)
hto_rows_L8 <- features_L8$V2 %in% hashtags

data_10x_L8_raw <- ReadMtx(
  mtx = here("L8/GSE230295_L8matrix.mtx.gz"),
  features = here("L8/GSE230295_L8features.tsv.gz"),
  cells = here("L8/GSE230295_L8barcodes.tsv.gz"),
  feature.column = 2
)
adt_L8 <- data_10x_L8_raw[adt_rows_L8, ]
colnames(adt_L8) <- paste0("L8_", colnames(adt_L8))
hto_L8 <- data_10x_L8_raw[hto_rows_L8, ]
colnames(hto_L8) <- paste0("L8_", colnames(hto_L8))
data_10x_L8 <- data_10x_L8_raw[gene_rows_L8, ]


#L9
features_L9 <- read.delim(here("L9", "GSE230295_L9features.tsv.gz"), header = FALSE)
gene_rows_L9 <- features_L9$V3 == "Gene Expression"
adt_rows_L9 <- features_L9$V3 == "Antibody Capture" & !(features_L9$V2 %in% hashtags)
hto_rows_L9 <- features_L9$V2 %in% hashtags

data_10x_L9_raw <- ReadMtx(
  mtx = here("L9/GSE230295_L9matrix.mtx.gz"),
  features = here("L9/GSE230295_L9features.tsv.gz"),
  cells = here("L9/GSE230295_L9barcodes.tsv.gz"),
  feature.column = 2
)
adt_L9 <- data_10x_L9_raw[adt_rows_L9, ]
colnames(adt_L9) <- paste0("L9_", colnames(adt_L9))
hto_L9 <- data_10x_L9_raw[hto_rows_L9, ]
colnames(hto_L9) <- paste0("L9_", colnames(hto_L9))
data_10x_L9 <- data_10x_L9_raw[gene_rows_L9, ]


#L10
features_L10 <- read.delim(here("L10", "GSE230295_L10features.tsv.gz"), header = FALSE)
gene_rows_L10 <- features_L10$V3 == "Gene Expression"
adt_rows_L10 <- features_L10$V3 == "Antibody Capture" & !(features_L10$V2 %in% hashtags)
hto_rows_L10 <- features_L10$V2 %in% hashtags

data_10x_L10_raw <- ReadMtx(
  mtx = here("L10/GSE230295_L10matrix.mtx.gz"),
  features = here("L10/GSE230295_L10features.tsv.gz"),
  cells = here("L10/GSE230295_L10barcodes.tsv.gz"),
  feature.column = 2
)
adt_L10 <- data_10x_L10_raw[adt_rows_L10, ]
colnames(adt_L10) <- paste0("L10_", colnames(adt_L10))
hto_L10 <- data_10x_L10_raw[hto_rows_L10, ]
colnames(hto_L10) <- paste0("L10_", colnames(hto_L10))
data_10x_L10 <- data_10x_L10_raw[gene_rows_L10, ]


#Comprobamos el tipo de datos que tenemos
dim(data_10x_L5)
head(rownames(data_10x_L5))
head(colnames(data_10x_L5))

dim(data_10x_L6)
head(rownames(data_10x_L6))
head(colnames(data_10x_L6))

dim(data_10x_L7)
head(rownames(data_10x_L7))
head(colnames(data_10x_L7))

dim(data_10x_L8)
head(rownames(data_10x_L8))
head(colnames(data_10x_L8))

dim(data_10x_L9)
head(rownames(data_10x_L9))
head(colnames(data_10x_L9))

dim(data_10x_L10)
head(rownames(data_10x_L10))
head(colnames(data_10x_L10))

#Creamos los objetos Seurat con los datos de cada librería (L5-L10) y añadimos el prefijo de la librería al barcode con RenameCells() para mantenerlas identificadas
seu_L5 <- CreateSeuratObject(counts = data_10x_L5,
                             project = 'L5')
seu_L5 <- RenameCells(seu_L5, add.cell.id = "L5")

seu_L6 <- CreateSeuratObject(counts = data_10x_L6,
                             project = 'L6')
seu_L6 <- RenameCells(seu_L6, add.cell.id = "L6")

seu_L7 <- CreateSeuratObject(counts = data_10x_L7,
                             project = 'L7')
seu_L7 <- RenameCells(seu_L7, add.cell.id = "L7")
seu_L8 <- CreateSeuratObject(counts = data_10x_L8,
                             project = 'L8')
seu_L8 <- RenameCells(seu_L8, add.cell.id = "L8")
seu_L9 <- CreateSeuratObject(counts = data_10x_L9,
                             project = 'L9')
seu_L9 <- RenameCells(seu_L9, add.cell.id = "L9")

seu_L10 <- CreateSeuratObject(counts = data_10x_L10,
                              project = 'L10')
seu_L10 <- RenameCells(seu_L10, add.cell.id = "L10")


# Asignamos metadatos a cada librería. Estos incluyen el origen de las muestras (médula ósea o sangre) que se son los mismos para todas las células de la librería

seu_L5$library <- "L5";  seu_L5$tissue <- "Bone_marrow"
seu_L6$library <- "L6";  seu_L6$tissue <- "Bone_marrow"
seu_L7$library <- "L7";  seu_L7$tissue <- "Bone_marrow"
seu_L8$library <- "L8";  seu_L8$tissue <- "Bone_marrow"
seu_L9$library <- "L9";  seu_L9$tissue <- "Blood"
seu_L10$library <- "L10";seu_L10$tissue <- "Blood"

#===============================================================================

###CONTROL DE CALIDAD###

#===============================================================================
#Eliminamos las gotas que están vacías y que tienen 0 conteos de ARN

seu_L5 <- subset(seu_L5, subset = nCount_RNA > 0)
seu_L6 <- subset(seu_L6, subset = nCount_RNA > 0)
seu_L7 <- subset(seu_L7, subset = nCount_RNA > 0)
seu_L8 <- subset(seu_L8, subset = nCount_RNA > 0)
seu_L9 <- subset(seu_L9, subset = nCount_RNA > 0)
seu_L10 <- subset(seu_L10, subset = nCount_RNA > 0)


#Comprobamos que se han incluido los datos y que tenemos las métricas de QC correspondientes a el nº de genes/célula, nº de UMIs/célula en el objeto de seurat

head(seu_L5@meta.data)
head(seu_L6@meta.data)
head(seu_L7@meta.data)
head(seu_L8@meta.data)
head(seu_L9@meta.data)
head(seu_L10@meta.data)

#Calculamos el % de transcritos(UMIs) mitocondriales ( (UMIs mitocondriales/UMIs totales)*100; indicativo de la viabilidad celular)

seu_L5[["percent.mt"]] <- PercentageFeatureSet(seu_L5, pattern = "^MT-")
seu_L6[["percent.mt"]] <- PercentageFeatureSet(seu_L6, pattern = "^MT-")
seu_L7[["percent.mt"]] <- PercentageFeatureSet(seu_L7, pattern = "^MT-")
seu_L8[["percent.mt"]] <- PercentageFeatureSet(seu_L8, pattern = "^MT-")
seu_L9[["percent.mt"]] <- PercentageFeatureSet(seu_L9, pattern = "^MT-")
seu_L10[["percent.mt"]] <- PercentageFeatureSet(seu_L10, pattern = "^MT-")



#Visualizamos las métricas de QC, es decir, respresentamos con histogramas de freceuncia los nFeatures= nºgenes/célula, 
#y nCounts = nº de transcritos detectados, así como el %genes mitocondriales
#También representamos la relación enter ncounts(transcritos) y nFeatures (genes) incorporando una recta de regresión lineal.
#Los esperable es una correlación positiva, es decir, a más moléculas de RNA más genes detectados


#L5
Histogram_nFeature_L5 <- ggplot(seu_L5@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nFeature_RNA - L5")

Histogram_nCount_L5 <- ggplot(seu_L5@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA - L5")

Histogram_mt_L5 <- ggplot(seu_L5@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("percent.mt - L5")

HistPlot5 <- Histogram_nFeature_L5 + Histogram_nCount_L5 + Histogram_mt_L5
HistPlot5

Scatterplot5 <- FeatureScatter(seu_L5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" ) +
  geom_smooth(method = "lm") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA vs nFeature_RNA -L5") #comprobamos si hay una relacion lineal entre el nº de tránscritos y de genes detectados, como QC
Scatterplot5 

#L6
Histogram_nFeature_L6 <- ggplot(seu_L6@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nFeature_RNA - L6")

Histogram_nCount_L6 <- ggplot(seu_L6@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA - L6")

Histogram_mt_L6 <- ggplot(seu_L6@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("percent.mt - L6")

HistPlot6 <- Histogram_nFeature_L6 + Histogram_nCount_L6 + Histogram_mt_L6
HistPlot6
Scatterplot6 <- FeatureScatter(seu_L6, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" ) +
  geom_smooth(method = "lm") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA vs nFeature_RNA -L6") 
Scatterplot6 

#L7
Histogram_nFeature_L7 <- ggplot(seu_L7@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nFeature_RNA - L7")

Histogram_nCount_L7 <- ggplot(seu_L7@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA - L7")

Histogram_mt_L7 <- ggplot(seu_L7@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("percent.mt - L7")

HistPlot7 <- Histogram_nFeature_L7 + Histogram_nCount_L7 + Histogram_mt_L7
HistPlot7


Scatterplot7 <- FeatureScatter(seu_L7, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" ) +
  geom_smooth(method = "lm") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA vs nFeature_RNA -L7") 
Scatterplot7 

#L8
Histogram_nFeature_L8 <- ggplot(seu_L8@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nFeature_RNA - L8")

Histogram_nCount_L8 <- ggplot(seu_L8@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA - L8")

Histogram_mt_L8 <- ggplot(seu_L8@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("percent.mt - L8")

HistPlot8 <- Histogram_nFeature_L8 + Histogram_nCount_L8 + Histogram_mt_L8
HistPlot8

Scatterplot8 <- FeatureScatter(seu_L8, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" ) +
  geom_smooth(method = "lm") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA vs nFeature_RNA -L8") 
Scatterplot8 

#L9

Histogram_nFeature_L9 <- ggplot(seu_L9@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nFeature_RNA - L9")

Histogram_nCount_L9 <- ggplot(seu_L9@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA - L9")

Histogram_mt_L9 <- ggplot(seu_L9@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("percent.mt - L9")

HistPlot9 <- Histogram_nFeature_L9 + Histogram_nCount_L9 + Histogram_mt_L9
HistPlot9

Scatterplot9 <- FeatureScatter(seu_L9, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" ) +
  geom_smooth(method = "lm")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA vs nFeature_RNA -L9") 
Scatterplot9 

#L10
Histogram_nFeature_L10 <- ggplot(seu_L10@meta.data, aes(x = nFeature_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nFeature_RNA - L10")

Histogram_nCount_L10 <- ggplot(seu_L10@meta.data, aes(x = nCount_RNA)) +
  geom_histogram(bins = 100) +
  scale_x_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA - L10")

Histogram_mt_L10 <- ggplot(seu_L10@meta.data, aes(x = percent.mt)) +
  geom_histogram(bins = 100) +
  scale_y_log10() +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("percent.mt - L10")

HistPlot10 <- Histogram_nFeature_L10 + Histogram_nCount_L10 + Histogram_mt_L10
HistPlot10

Scatterplot10 <- FeatureScatter(seu_L10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA" ) +
  geom_smooth(method = "lm")+
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11)
  ) +
  ggtitle("nCount_RNA vs nFeature_RNA -L10") 
Scatterplot10 



#Guardamos los plots
ggsave(here("plots1/HistPlot_L5.png"), plot = HistPlot5, width = 12, height = 4, dpi = 300)
ggsave(here("plots1/HistPlot_L6.png"), plot = HistPlot6, width = 12, height = 4, dpi = 300)
ggsave(here("plots1/HistPlot_L7.png"), plot = HistPlot7, width = 12, height = 4, dpi = 300)
ggsave(here("plots1/HistPlot_L8.png"), plot = HistPlot8, width = 12, height = 4, dpi = 300)
ggsave(here("plots1/HistPlot_L9.png"), plot = HistPlot9, width = 12, height = 4, dpi = 300)
ggsave(here("plots1/HistPlot_L10.png"), plot = HistPlot10, width = 12, height = 4, dpi = 300)
ggsave(here("plots1/Scatterplot_L5.png"), plot = Scatterplot5, width = 6, height = 5, dpi = 300)
ggsave(here("plots1/Scatterplot_L6.png"), plot = Scatterplot6, width = 6, height = 5, dpi = 300)
ggsave(here("plots1/Scatterplot_L7.png"), plot = Scatterplot7, width = 6, height = 5, dpi = 300)
ggsave(here("plots1/Scatterplot_L8.png"), plot = Scatterplot8, width = 6, height = 5, dpi = 300)
ggsave(here("plots1/Scatterplot_L9.png"), plot = Scatterplot9, width = 6, height = 5, dpi = 300)
ggsave(here("plots1/Scatterplot_L10.png"), plot = Scatterplot10, width = 6, height = 5, dpi = 300)


#Creamos una lista con todas las librerias para poder aplicar los siguientes pasos a todas a la vez

seu_list <- list(L5 = seu_L5, L6 = seu_L6, L7 = seu_L7, L8 = seu_L8, L9 = seu_L9, L10 = seu_L10) 

# Se establecen los umbrales de QC, considerandose células de "buena calidad aquellas en las que se detectan:
#>200 y < 6000 genes y un porcentaje de transcritos mitocondriales < 10%)


#Calculamos el % de células excluidas en cada librería con estos criterios
qc_exclusion_by_library <- lapply(names(seu_list), function(nm) {
  obj <- seu_list[[nm]]
  
  candidates <- subset(obj, subset = nFeature_RNA > 200)
  
  data.frame(
    library = nm,
    n_candidates = ncol(candidates),
    pct_excluidos_por_max_feat = round(sum(candidates$nFeature_RNA >= 6000) / ncol(candidates) * 100, 2),
    pct_excluidos_por_mt = round(sum(candidates$percent.mt >= 10) / ncol(candidates) * 100, 2)
  )
})

qc_exclusion_by_library <- do.call(rbind, qc_exclusion_by_library)
qc_exclusion_by_library

write.csv(qc_exclusion_by_library, file = here("qc_exclusión_por_libreria.csv"), row.names = FALSE)


#Filtramos para conservar sólo las celulas con entre 200 y 6000 genes y menos de un 10% de genes mitocondriales

qc_filter <- function(obj, min_feat = 200, max_feat = 6000, max_mt = 10) {
  subset(obj, subset = nFeature_RNA > min_feat & nFeature_RNA < max_feat & percent.mt < max_mt)
} #revisar estos parámetros y ajustar

seu_list_filtered <- lapply(seu_list, qc_filter) #aplicamos filtro a todas las librerías


# comprobamos cuantas células había antes y después del filtrado en cada librería
qc_summary <- data.frame(
  library = names(seu_list),
  cells_before = sapply(seu_list, ncol),
  cells_after  = sapply(seu_list_filtered, ncol)
)

# Añadimos el porcentaje de células retenidas
qc_summary$percent_retained <- round(
  qc_summary$cells_after / qc_summary$cells_before * 100,
  2
)

qc_summary
write.csv(qc_summary, file = here("qc_summary.csv"), row.names = FALSE)

#Obtenemos datos de las células,genes y UMIs de cada librería tras el QC

summary_postQC <- rbind(
  data.frame(
    Libreria = "L5",
    Celulas = ncol(seu_list_filtered$L5),
    Media_genes = round(mean(seu_list_filtered$L5$nFeature_RNA), 1),
    Mediana_genes = median(seu_list_filtered$L5$nFeature_RNA),
    Media_UMIs = round(mean(seu_list_filtered$L5$nCount_RNA), 1)
  ),

  data.frame(
    Libreria = "L6",
    Celulas = ncol(seu_list_filtered$L6),
    Media_genes = round(mean(seu_list_filtered$L6$nFeature_RNA), 1),
    Mediana_genes = median(seu_list_filtered$L6$nFeature_RNA),
    Media_UMIs = round(mean(seu_list_filtered$L6$nCount_RNA), 1)
  ),

  data.frame(
    Libreria = "L7",
    Celulas = ncol(seu_list_filtered$L7),
    Media_genes = round(mean(seu_list_filtered$L7$nFeature_RNA), 1),
    Mediana_genes = median(seu_list_filtered$L7$nFeature_RNA),
    Media_UMIs = round(mean(seu_list_filtered$L7$nCount_RNA), 1)
  ),

  data.frame(
    Libreria = "L8",
    Celulas = ncol(seu_list_filtered$L8),
    Media_genes = round(mean(seu_list_filtered$L8$nFeature_RNA), 1),
    Mediana_genes = median(seu_list_filtered$L8$nFeature_RNA),
    Media_UMIs = round(mean(seu_list_filtered$L8$nCount_RNA), 1)
  ), 
  
  data.frame(
    Libreria = "L9",
    Celulas = ncol(seu_list_filtered$L9),
    Media_genes = round(mean(seu_list_filtered$L9$nFeature_RNA), 1),
    Mediana_genes = median(seu_list_filtered$L9$nFeature_RNA),
    Media_UMIs = round(mean(seu_list_filtered$L9$nCount_RNA), 1)
  ),
  
  data.frame(
    Libreria = "L10",
    Celulas = ncol(seu_list_filtered$L10),
    Media_genes = round(mean(seu_list_filtered$L10$nFeature_RNA), 1),
    Mediana_genes = median(seu_list_filtered$L10$nFeature_RNA),
    Media_UMIs = round(mean(seu_list_filtered$L10$nCount_RNA), 1)
  )
)

summary_postQC
write.csv( summary_postQC,"summary_postQC.csv",row.names = FALSE)

#Detección de dobletes con scDblFinder()

set.seed(1234) #para que la parte aleatoria de scDblFinder se reproducible

for (nm in names(seu_list_filtered)) {
  
  obj <- seu_list_filtered[[nm]]
  
  sce <- as.SingleCellExperiment(obj) #scDblFinder necesita un SingleCellExperiment, no un objeto Seurat
  sce <- scDblFinder(sce)
  
  obj$scDblFinder.class <- sce$scDblFinder.class #"singlet" o "doublet"
  
  seu_list_filtered[[nm]] <- obj
}

#Verificamos los dobletes
doublet_summary <- data.frame(
  library = names(seu_list_filtered),
  n_cells = sapply(seu_list_filtered, ncol),
  n_doublets = sapply(seu_list_filtered, function(x) sum(x$scDblFinder.class == "doublet"))
)
doublet_summary$percent_doublets <- round(doublet_summary$n_doublets / doublet_summary$n_cells * 100, 2)
doublet_summary
write.csv(doublet_summary, file = here("doublet_summary.csv"), row.names = FALSE)


#Visualizamos los dobletes, para ver como se reparten en el UMAP

for (nm in names(seu_list_filtered)) {
  
  obj <- seu_list_filtered[[nm]]
  
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 20, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:20, verbose = FALSE)
  
  Plot_dobletes <- DimPlot(obj, group.by = "scDblFinder.class", reduction = "umap") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11)
    ) +
    ggtitle(paste0("Dobletes - ", nm))
  print(Plot_dobletes)
  
  ggsave(here("plots1", paste0("QC_doublets_", nm, ".png")), plot = Plot_dobletes, width = 7, height = 5, dpi = 300)
}

#Eliminación de dobletes
for (nm in names(seu_list_filtered)) {
  
  obj <- seu_list_filtered[[nm]]
  obj <- subset(obj, subset = scDblFinder.class == "singlet")
  seu_list_filtered[[nm]] <- obj
  
}

#Comprobamos tras la eliminación
doublet_summary$n_after_doublet_removal <- sapply(seu_list_filtered, ncol)
doublet_summary

write.csv(doublet_summary, file = here("doublet_summary.csv"), row.names = FALSE)


#Guardamos los datos: seu_list_filtered (ARN, ya con QC y eliminación de dobletes aplicados)
#listas de ADT (proteína) y HTO (hashtags) de las 6 librerías

adt_list <- list(L5 = adt_L5, L6 = adt_L6, L7 = adt_L7, L8 = adt_L8, L9 = adt_L9, L10 = adt_L10)
hto_list <- list(L5 = hto_L5, L6 = hto_L6, L7 = hto_L7, L8 = hto_L8, L9 = hto_L9, L10 = hto_L10)

saveRDS(seu_list_filtered, file = here("1_seu_list_filtered.rds"))
saveRDS(adt_list, file = here("1_adt_list.rds"))
saveRDS(hto_list, file = here("1_hto_list.rds"))

#Información de la sesión
sessionInfo()

