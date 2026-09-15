#===============================================================================
# # TFM - Análisis integrado de la expresión génica y de proteínas de superficie 
# a nivel unicelular en muestras de médula ósea y sangre periférica de pacientes 
# con LLA-B (GSE230295)
#
# DEMULTIPLEXADO:
# Asignación de identificacion del paciente donante a cada muestra y del momento 
# de obtención de esta (diagnóstico / 15 días / 2 días) 
# A partir de la señal de hashtag (HTO) y de los metadatos de muestra individual (GSM) de GEO.
# Se añade "donor" y "timepoint" como metadato a cada célula de cada librería,
#
#===============================================================================

library(Seurat)
library(here)

#===============================================================================

#Cargamos datos

seu_list_filtered <- readRDS(here("1_seu_list_filtered.rds"))
hto_list <- readRDS(here("1_hto_list.rds"))

hto_L5 <- hto_list$L5
hto_L6 <- hto_list$L6
hto_L7 <- hto_list$L7
hto_L8 <- hto_list$L8
hto_L9 <- hto_list$L9
hto_L10 <- hto_list$L10


#Tabla manual de correspondencia hashtag -> paciente -> momento clínico, obtenida a partir de
#los registros GSM de cada muestra de GEO (2 por librería, para las 6 librerías)
tabla_hashtags <- data.frame(
  library = c("L5","L5","L6","L6","L7","L7","L8","L8","L9","L9","L10","L10"),
  hashtag = c("Hash1","Hash2","Hash3","Hash4","Hash1","Hash2","Hash3","Hash4","Hash1","Hash2","Hash3","Hash4"),
  donor = c("D5","D6","D7","D8","D6","D5","D8","D7","D8","D7","D5","D7"),
  timepoint = c("Dx","15d","Dx","15d","Dx","15d","Dx","15d","Dx","2d","2d","Dx"),
  stringsAsFactors = FALSE
)
tabla_hashtags


#HTODemux por librería, usando solo los 2 hashtags reales de esa librería (según tabla_hashtags) -

#L5
hashtags_L5 <- tabla_hashtags$hashtag[tabla_hashtags$library == "L5"]
hto_L5_qc <- hto_L5[hashtags_L5, colnames(seu_list_filtered$L5)]
celulas_con_senal_L5 <- colnames(hto_L5_qc)[colSums(hto_L5_qc) > 0]

seu_L5_temp <- subset(seu_list_filtered$L5, cells = celulas_con_senal_L5) #seleccionamos sólo las células con señal y creamos un objeto temporal para el HTODemux()
hto_L5_qc <- hto_L5_qc[, celulas_con_senal_L5]
seu_L5_temp[["HTO"]] <- CreateAssayObject(counts = hto_L5_qc) #añadimos el hashtag al objeto seurat
seu_L5_temp <- NormalizeData(seu_L5_temp, assay = "HTO", normalization.method = "CLR") #Normalización de los datos
seu_L5_temp <- HTODemux(seu_L5_temp, assay = "HTO", positive.quantile = 0.99)
table(seu_L5_temp$HTO_classification.global) #Se cuenta las células incluidas en cada categoría

#Representación de los patrones de expresión en Ridgeplot
L5_RidgePlot <- RidgePlot(seu_L5_temp, assay = "HTO", features = rownames(seu_L5_temp[["HTO"]]),ncol = 2)
ggsave(here("plots1/L5_HTO_ridgeplot.png"),plot = L5_RidgePlot, width = 8, height = 4, dpi = 300)


#Traducimos hash.ID a paciente/momento clínico y lo añadimos a la librería completa (no solo
#a las células con señal de HTO) - las que no tienen dato quedan con NA, pero siguen en el
#análisis principal de tipos celulares, solo se excluyen de la comparación por momento clínico

mapa_L5 <- data.frame(cell = colnames(seu_L5_temp), hashtag = as.character(seu_L5_temp$hash.ID), stringsAsFactors = FALSE)
mapa_L5 <- merge(mapa_L5, tabla_hashtags[tabla_hashtags$library == "L5", ], by = "hashtag", all.x = TRUE)

posiciones_L5 <- match(colnames(seu_list_filtered$L5), mapa_L5$cell) #incorporamos las equivalencias de los hashtags al objeto seurat original
seu_list_filtered$L5$donor <- mapa_L5$donor[posiciones_L5]
seu_list_filtered$L5$timepoint <- mapa_L5$timepoint[posiciones_L5]
table(seu_list_filtered$L5$timepoint, useNA = "ifany")


#L6
hashtags_L6 <- tabla_hashtags$hashtag[tabla_hashtags$library == "L6"]
hto_L6_qc <- hto_L6[hashtags_L6, colnames(seu_list_filtered$L6)]
celulas_con_senal_L6 <- colnames(hto_L6_qc)[colSums(hto_L6_qc) > 0]

seu_L6_temp <- subset(seu_list_filtered$L6, cells = celulas_con_senal_L6)
hto_L6_qc <- hto_L6_qc[, celulas_con_senal_L6]
seu_L6_temp[["HTO"]] <- CreateAssayObject(counts = hto_L6_qc)
seu_L6_temp <- NormalizeData(seu_L6_temp, assay = "HTO", normalization.method = "CLR")
seu_L6_temp <- HTODemux(seu_L6_temp, assay = "HTO", positive.quantile = 0.99)
table(seu_L6_temp$HTO_classification.global)


L6_RidgePlot <- RidgePlot(seu_L6_temp, assay = "HTO", features = rownames(seu_L6_temp[["HTO"]]),ncol = 2)
L6_RidgePlot 
ggsave(here("plots1/L6_HTO_ridgeplot.png"),plot = L6_RidgePlot, width = 8, height = 4, dpi = 300)


mapa_L6 <- data.frame(cell = colnames(seu_L6_temp), hashtag = as.character(seu_L6_temp$hash.ID), stringsAsFactors = FALSE)
mapa_L6 <- merge(mapa_L6, tabla_hashtags[tabla_hashtags$library == "L6", ], by = "hashtag", all.x = TRUE)

posiciones_L6 <- match(colnames(seu_list_filtered$L6), mapa_L6$cell)
seu_list_filtered$L6$donor <- mapa_L6$donor[posiciones_L6]
seu_list_filtered$L6$timepoint <- mapa_L6$timepoint[posiciones_L6]
table(seu_list_filtered$L6$timepoint, useNA = "ifany")


#L7
hashtags_L7 <- tabla_hashtags$hashtag[tabla_hashtags$library == "L7"]
hto_L7_qc <- hto_L7[hashtags_L7, colnames(seu_list_filtered$L7)]
celulas_con_senal_L7 <- colnames(hto_L7_qc)[colSums(hto_L7_qc) > 0]

seu_L7_temp <- subset(seu_list_filtered$L7, cells = celulas_con_senal_L7)
hto_L7_qc <- hto_L7_qc[, celulas_con_senal_L7]
seu_L7_temp[["HTO"]] <- CreateAssayObject(counts = hto_L7_qc)
seu_L7_temp <- NormalizeData(seu_L7_temp, assay = "HTO", normalization.method = "CLR")
seu_L7_temp <- HTODemux(seu_L7_temp, assay = "HTO", positive.quantile = 0.99)
table(seu_L7_temp$HTO_classification.global)


L7_RidgePlot <- RidgePlot(seu_L7_temp, assay = "HTO", features = rownames(seu_L7_temp[["HTO"]]),ncol = 2)
L7_RidgePlot
ggsave(here("plots1/L7_HTO_ridgeplot.png"),plot = L7_RidgePlot, width = 8, height = 4, dpi = 300)


mapa_L7 <- data.frame(cell = colnames(seu_L7_temp), hashtag = as.character(seu_L7_temp$hash.ID), stringsAsFactors = FALSE)
mapa_L7 <- merge(mapa_L7, tabla_hashtags[tabla_hashtags$library == "L7", ], by = "hashtag", all.x = TRUE)

posiciones_L7 <- match(colnames(seu_list_filtered$L7), mapa_L7$cell)
seu_list_filtered$L7$donor <- mapa_L7$donor[posiciones_L7]
seu_list_filtered$L7$timepoint <- mapa_L7$timepoint[posiciones_L7]
table(seu_list_filtered$L7$timepoint, useNA = "ifany")


#L8: EXCLUIDA del demultiplexado, ya que con HTODemux() da el error "Cells with zero counts
#exist as a cluster". La señal de hashtag en esta librería es, en conjunto, más débil que en el resto
# que podría ser una explicación técnica plausible para este fallo del código.

hashtags_L8 <- tabla_hashtags$hashtag[tabla_hashtags$library == "L8"]
hto_L8_qc <- hto_L8[hashtags_L8, colnames(seu_list_filtered$L8)]
celulas_con_senal_L8 <- colnames(hto_L8_qc)[colSums(hto_L8_qc) > 0]

seu_L8_temp <- subset(seu_list_filtered$L8, cells = celulas_con_senal_L8)
hto_L8_qc <- hto_L8_qc[, celulas_con_senal_L8]
seu_L8_temp[["HTO"]] <- CreateAssayObject(counts = hto_L8_qc)
seu_L8_temp <- NormalizeData(seu_L8_temp, assay = "HTO", normalization.method = "CLR")
seu_L8_temp <- HTODemux(seu_L8_temp, assay = "HTO", positive.quantile = 0.99) #no se puede hacer el demultiplexado. Error:Cells with zero counts exist as a cluster

L8_RidgePlot <- RidgePlot(seu_L8_temp, assay = "HTO", features = rownames(seu_L8_temp[["HTO"]]),ncol = 2)
L8_RidgePlot 
table(
  Hash3 = hto_L8_qc["Hash3", ] > 0,
  Hash4 = hto_L8_qc["Hash4", ] > 0
) #La mayoría de las células expresan ambos hashtags a la vez (55%)


# Las células de L8 se quedan sin donor/real_timepoint (NA)
seu_list_filtered$L8$donor <- NA
seu_list_filtered$L8$timepoint <- NA

#L9
hashtags_L9 <- tabla_hashtags$hashtag[tabla_hashtags$library == "L9"]
hto_L9_qc <- hto_L9[hashtags_L9, colnames(seu_list_filtered$L9)]
celulas_con_senal_L9 <- colnames(hto_L9_qc)[colSums(hto_L9_qc) > 0]

seu_L9_temp <- subset(seu_list_filtered$L9, cells = celulas_con_senal_L9)
hto_L9_qc <- hto_L9_qc[, celulas_con_senal_L9]
seu_L9_temp[["HTO"]] <- CreateAssayObject(counts = hto_L9_qc)
seu_L9_temp <- NormalizeData(seu_L9_temp, assay = "HTO", normalization.method = "CLR")
seu_L9_temp <- HTODemux(seu_L9_temp, assay = "HTO", positive.quantile = 0.99)
table(seu_L9_temp$HTO_classification.global)


L9_RidgePlot <- RidgePlot(seu_L9_temp, assay = "HTO", features = rownames(seu_L9_temp[["HTO"]]),ncol = 2)
L9_RidgePlot
ggsave(here("plots1/L9_HTO_ridgeplot.png"),plot = L9_RidgePlot, width = 8, height = 4, dpi = 300)


mapa_L9 <- data.frame(cell = colnames(seu_L9_temp), hashtag = as.character(seu_L9_temp$hash.ID), stringsAsFactors = FALSE)
mapa_L9 <- merge(mapa_L9, tabla_hashtags[tabla_hashtags$library == "L9", ], by = "hashtag", all.x = TRUE)

posiciones_L9 <- match(colnames(seu_list_filtered$L9), mapa_L9$cell)
seu_list_filtered$L9$donor <- mapa_L9$donor[posiciones_L9]
seu_list_filtered$L9$timepoint <- mapa_L9$timepoint[posiciones_L9]
table(seu_list_filtered$L9$timepoint, useNA = "ifany")


#L10
hashtags_L10 <- tabla_hashtags$hashtag[tabla_hashtags$library == "L10"]
hto_L10_qc <- hto_L10[hashtags_L10, colnames(seu_list_filtered$L10)]
celulas_con_senal_L10 <- colnames(hto_L10_qc)[colSums(hto_L10_qc) > 0]

seu_L10_temp <- subset(seu_list_filtered$L10, cells = celulas_con_senal_L10)
hto_L10_qc <- hto_L10_qc[, celulas_con_senal_L10]
seu_L10_temp[["HTO"]] <- CreateAssayObject(counts = hto_L10_qc)
seu_L10_temp <- NormalizeData(seu_L10_temp, assay = "HTO", normalization.method = "CLR")
seu_L10_temp <- HTODemux(seu_L10_temp, assay = "HTO", positive.quantile = 0.99)
table(seu_L10_temp$HTO_classification.global)


L10_RidgePlot <- RidgePlot(seu_L10_temp, assay = "HTO", features = rownames(seu_L10_temp[["HTO"]]),ncol = 2)
L10_RidgePlot
ggsave(here("plots1/L10_HTO_ridgeplot.png"),plot = L10_RidgePlot, width = 8, height = 4, dpi = 300)

mapa_L10 <- data.frame(cell = colnames(seu_L10_temp), hashtag = as.character(seu_L10_temp$hash.ID), stringsAsFactors = FALSE)
mapa_L10 <- merge(mapa_L10, tabla_hashtags[tabla_hashtags$library == "L10", ], by = "hashtag", all.x = TRUE)

posiciones_L10 <- match(colnames(seu_list_filtered$L10), mapa_L10$cell)
seu_list_filtered$L10$donor <- mapa_L10$donor[posiciones_L10]
seu_list_filtered$L10$timepoint <- mapa_L10$timepoint[posiciones_L10]
table(seu_list_filtered$L10$timepoint, useNA = "ifany")


# Resumen de clasificación HTO (singlet/negative/doublet) en porcentaje de las 5 librerías demultiplexadas
# (L8 no entra porque no se pudo demultiplexar)

resumen_L5 <- as.data.frame(round(100 * table(seu_L5_temp$HTO_classification.global) / length(seu_L5_temp$HTO_classification.global), 1))
resumen_L6 <- as.data.frame(round(100 * table(seu_L6_temp$HTO_classification.global) / length(seu_L6_temp$HTO_classification.global), 1))
resumen_L7 <- as.data.frame(round(100 * table(seu_L7_temp$HTO_classification.global) / length(seu_L7_temp$HTO_classification.global), 1))
resumen_L9 <- as.data.frame(round(100 * table(seu_L9_temp$HTO_classification.global) / length(seu_L9_temp$HTO_classification.global), 1))
resumen_L10 <- as.data.frame(round(100 * table(seu_L10_temp$HTO_classification.global) / length(seu_L10_temp$HTO_classification.global), 1))

colnames(resumen_L5) <- c("Clasificación", "L5")
colnames(resumen_L6) <- c("Clasificación", "L6")
colnames(resumen_L7) <- c("Clasificación", "L7")
colnames(resumen_L9) <- c("Clasificación", "L9")
colnames(resumen_L10) <- c("Clasificación", "L10")

resumen_demultiplexado <- merge(resumen_L5, resumen_L6, by = "Clasificación")
resumen_demultiplexado <- merge(resumen_demultiplexado, resumen_L7, by = "Clasificación")
resumen_demultiplexado <- merge(resumen_demultiplexado, resumen_L9, by = "Clasificación")
resumen_demultiplexado <- merge(resumen_demultiplexado, resumen_L10, by = "Clasificación")
resumen_demultiplexado

write.csv(resumen_demultiplexado, file = here("resumen_demultiplexado.csv"), row.names = FALSE)


#Guardamos seu_list_filtered, ya con donor y timepoint añadidos
saveRDS(seu_list_filtered, file = here("2_seu_list_filtered_demux.rds"))


#Información de la sesión
sessionInfo()

