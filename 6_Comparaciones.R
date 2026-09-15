#===============================================================================
# # TFM - Análisis integrado de la expresión génica y de proteínas de superficie 
# a nivel unicelular en muestras de médula ósea y sangre periférica de pacientes 
# con LLA-B (GSE230295)
#
# Comparación de la composición de tipos celulares entre médula ósea y sangre
# periférica, y entre diagnóstico y post-tratamiento usando "donor"/timepoint" tanto de forma
# agregada como por paciente individual, contrastando con el blast_cell_percentage
# clínico de los metadatos GSM de GEO.
#

#===============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(here)

#===============================================================================

#Cargamos datos

dataset_integrado <- readRDS(here("5_dataset_integrado_validado.rds"))

#Paleta de colores por tipo celular 

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


#COMPARACIÓN DE COMPOSICIÓN CELULAR: MÉDULA ÓSEA VS SANGRE

#Tabla de conteos: tipo celular x tejido
tabla_celltype_tissue <- table(dataset_integrado$cell_type, dataset_integrado$tissue)
tabla_celltype_tissue

#Proporción de cada tipo celular en de cada tejido (para poder comparar composición aunque
#médula y sangre tengan un nº total de células distinto)
prop_celltype_tissue <- prop.table(tabla_celltype_tissue, margin = 2) * 100
prop_celltype_tissue

prop_celltype_tissue_df <- as.data.frame(prop_celltype_tissue)
colnames(prop_celltype_tissue_df) <- c("cell_type", "tissue", "percentage")

write.csv(prop_celltype_tissue_df, file = here("composicion_celltype_por_tejido.csv"), row.names = FALSE)


#Gráfico de barras apiladas: composición celular por tejido
plot_composicion_tejido <- ggplot(prop_celltype_tissue_df, aes(x = factor(tissue, levels = c("Bone_marrow", "Blood")), y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.3) +
  scale_fill_manual(values = celltype_colors) +
  scale_x_discrete(labels = c("Bone_marrow" = "Médula ósea", "Blood" = "Sangre")) +
  theme_minimal(base_size = 12) +
    theme(
    plot.title = element_text(face = "bold", margin = margin(b = 15), hjust = 0.5),
    plot.margin = margin(t = 15, r = 10, b = 10, l = 10),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.key.size = unit(0.4, "cm"),
    axis.text.x = element_text(size = 14)
  ) +
  labs(x = "Tejido", y = "% de células", fill = "Tipo celular",
       title = "Composición celular por tejido")

plot_composicion_tejido

ggsave(
  filename = here("plots", "composicion_celltype_por_tejido.png"),
  plot = plot_composicion_tejido,
  width = 11,
  height = 8,
  dpi = 300
)

ggsave(filename = here("plots1", "composicion_celltype_por_tejido.png"), plot = plot_composicion_tejido, width = 10, height = 8, dpi = 300)


#Test estadístico: ¿la composición celular difiere significativamente entre médula y sangre?

chisq_test_tissue <- chisq.test(tabla_celltype_tissue)
chisq_test_tissue


#Tabla de diferencias: para cada tipo celular, su % en médula vs en sangre, ordenada por
#la diferencia absoluta 
bm_pct <- prop_celltype_tissue_df[prop_celltype_tissue_df$tissue == "Bone_marrow", c("cell_type", "percentage")]
blood_pct <- prop_celltype_tissue_df[prop_celltype_tissue_df$tissue == "Blood", c("cell_type", "percentage")]
colnames(bm_pct)[2] <- "pct_bone_marrow"
colnames(blood_pct)[2] <- "pct_blood"

comparacion_tejidos <- merge(bm_pct, blood_pct, by = "cell_type")
comparacion_tejidos$diferencia_abs <- abs(comparacion_tejidos$pct_bone_marrow - comparacion_tejidos$pct_blood)
comparacion_tejidos <- comparacion_tejidos[order(-comparacion_tejidos$diferencia_abs), ]

comparacion_tejidos

write.csv(comparacion_tejidos, file = here("comparacion_celltype_tejido_diferencias.csv"), row.names = FALSE)


#COMPARACIÓN DE COMPOSICIÓN CELULAR: DIAGNÓSTICO VS. POST-TRATAMIENTO (timepoint)


#Vemos primero, paciente a paciente, qué combinaciones de tejido/momento clínico tenemos 

table(dataset_integrado$donor, dataset_integrado$timepoint, useNA = "ifany")


# MÉDULA ÓSEA: Diagnóstico vs. Día 15
datos_medula_timepoint <- subset(dataset_integrado, subset = tissue == "Bone_marrow" & !is.na(timepoint))

tabla_celltype_medula_tp <- table(datos_medula_timepoint$cell_type, datos_medula_timepoint$timepoint)
tabla_celltype_medula_tp

prop_celltype_medula_tp <- prop.table(tabla_celltype_medula_tp, margin = 2) * 100
prop_celltype_medula_tp_df <- as.data.frame(prop_celltype_medula_tp)
colnames(prop_celltype_medula_tp_df) <- c("cell_type", "timepoint", "percentage")

write.csv(prop_celltype_medula_tp_df, file = here("composicion_celltype_medula_dx_vs_15d.csv"), row.names = FALSE)

#Contamos cuántos pacientes distintos hay detrás de cada grupo 
donantes_medula <- unique(data.frame(donor = datos_medula_timepoint$donor, timepoint = datos_medula_timepoint$timepoint))
donantes_medula
n_pacientes_medula <- table(donantes_medula$timepoint)
n_pacientes_medula

composicion_medula_tp <- ggplot(prop_celltype_medula_tp_df, aes(x = factor(timepoint, levels = c("Dx", "15d")), y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.3) +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", margin = margin(b = 15), hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.key.size = unit(0.4, "cm")
  ) +
  labs(
    x = "Momento clínico", y = "% de células", fill = "Tipo celular",
    title = "Composición celular: médula ósea, diagnóstico vs. día 15",
    subtitle = paste0("n = ", n_pacientes_medula["Dx"], " pacientes (Dx), n = ", n_pacientes_medula["15d"], " pacientes (15d)")
  )

composicion_medula_tp

ggsave(filename = here("plots1", "composicion_medula_dx_vs_15d.png"), plot = composicion_medula_tp, width = 10, height = 8, dpi = 300)


# SANGRE: Diagnóstico vs. Día 2
datos_sangre_timepoint <- subset(dataset_integrado, subset = tissue == "Blood" & !is.na(timepoint))

tabla_celltype_sangre_tp <- table(datos_sangre_timepoint$cell_type, datos_sangre_timepoint$timepoint)
tabla_celltype_sangre_tp

prop_celltype_sangre_tp <- prop.table(tabla_celltype_sangre_tp, margin = 2) * 100
prop_celltype_sangre_tp_df <- as.data.frame(prop_celltype_sangre_tp)
colnames(prop_celltype_sangre_tp_df) <- c("cell_type", "timepoint", "percentage")

write.csv(prop_celltype_sangre_tp_df, file = here("composicion_celltype_sangre_dx_vs_2d.csv"), row.names = FALSE)

donantes_sangre <- unique(data.frame(donor = datos_sangre_timepoint$donor, timepoint = datos_sangre_timepoint$timepoint))
donantes_sangre
n_pacientes_sangre <- table(donantes_sangre$timepoint)
n_pacientes_sangre

composicion_sangre_tp <- ggplot(prop_celltype_sangre_tp_df, aes(x = factor(timepoint, levels = c("Dx", "2d")), y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.3) +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", margin = margin(b = 15), hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 18),
    legend.key.size = unit(0.4, "cm")
  ) +
  labs(
    x = "Momento clínico", y = "% de células", fill = "Tipo celular",
    title = "Composición celular: sangre, diagnóstico vs. día 2",
    subtitle = paste0("n = ", n_pacientes_sangre["Dx"], " pacientes (Dx), n = ", n_pacientes_sangre["2d"], " pacientes (2d)")
  )

composicion_sangre_tp

ggsave(filename = here("plots1", "composicion_sangre_dx_vs_2d.png"), plot = composicion_sangre_tp, width = 10, height = 8, dpi = 300)


#Tablas de diferencias (similares a las que hemos hecho para médula vs sangre)

#Médula: Dx vs 15d
dx_pct_medula <- prop_celltype_medula_tp_df[prop_celltype_medula_tp_df$timepoint == "Dx", c("cell_type","percentage")]
d15_pct_medula <- prop_celltype_medula_tp_df[prop_celltype_medula_tp_df$timepoint == "15d", c("cell_type","percentage")]
colnames(dx_pct_medula)[2] <- "pct_Dx"
colnames(d15_pct_medula)[2] <- "pct_15d"
comparacion_medula_tp <- merge(dx_pct_medula, d15_pct_medula, by = "cell_type")
comparacion_medula_tp$diferencia_abs <- abs(comparacion_medula_tp$pct_Dx - comparacion_medula_tp$pct_15d)
comparacion_medula_tp <- comparacion_medula_tp[order(-comparacion_medula_tp$diferencia_abs), ]
comparacion_medula_tp

write.csv(comparacion_medula_tp, file = here("comparacion_medula_dx_15d_diferencias.csv"), row.names = FALSE)

#Sangre: Dx vs 2d
dx_pct_sangre <- prop_celltype_sangre_tp_df[prop_celltype_sangre_tp_df$timepoint == "Dx", c("cell_type","percentage")]
d2_pct_sangre <- prop_celltype_sangre_tp_df[prop_celltype_sangre_tp_df$timepoint == "2d", c("cell_type","percentage")]
colnames(dx_pct_sangre)[2] <- "pct_Dx"
colnames(d2_pct_sangre)[2] <- "pct_2d"
comparacion_sangre_tp <- merge(dx_pct_sangre, d2_pct_sangre, by = "cell_type")
comparacion_sangre_tp$diferencia_abs <- abs(comparacion_sangre_tp$pct_Dx - comparacion_sangre_tp$pct_2d)
comparacion_sangre_tp <- comparacion_sangre_tp[order(-comparacion_sangre_tp$diferencia_abs), ]
comparacion_sangre_tp

write.csv(comparacion_sangre_tp, file = here("comparacion_sangre_dx_2d_diferencias.csv"), row.names = FALSE)


#COMPOSICIÓN POR PACIENTE 
#Miramos cada paciente por separado y comparamos también con blast_cell_percentage,
#el dato clínico real reportado en los metadatos GSM de GEO.

categorias_blasto <- c(
  "Progenitor B comprometido (LyP-B/pro-B)",
  "Progenitor B comprometido, proliferando",
  "Progenitor B muy proliferativo (G2/M)",
  "Progenitor pro-B/pre-B en recombinación V(D)J"
)


#SANGRE, por paciente 

datos_D7_dx_sangre <- subset(dataset_integrado, subset = donor == "D7" & timepoint == "Dx" & tissue == "Blood")
datos_D7_2d <- subset(dataset_integrado, subset = donor == "D7" & timepoint == "2d")
datos_D8_dx_sangre <- subset(dataset_integrado, subset = donor == "D8" & timepoint == "Dx" & tissue == "Blood")
datos_D5_2d <- subset(dataset_integrado, subset = donor == "D5" & timepoint == "2d")

prop_D7_dx_sangre <- round(prop.table(table(datos_D7_dx_sangre$cell_type)) * 100, 2)
prop_D7_2d <- round(prop.table(table(datos_D7_2d$cell_type)) * 100, 2)
prop_D8_dx_sangre <- round(prop.table(table(datos_D8_dx_sangre$cell_type)) * 100, 2)
prop_D5_2d <- round(prop.table(table(datos_D5_2d$cell_type)) * 100, 2)

blasto_D7_dx_sangre <- sum(prop_D7_dx_sangre[names(prop_D7_dx_sangre) %in% categorias_blasto])
blasto_D7_2d <- sum(prop_D7_2d[names(prop_D7_2d) %in% categorias_blasto])
blasto_D8_dx_sangre <- sum(prop_D8_dx_sangre[names(prop_D8_dx_sangre) %in% categorias_blasto])
blasto_D5_2d <- sum(prop_D5_2d[names(prop_D5_2d) %in% categorias_blasto])

#blast_percentage_clinico_GEO viene de los registros GSM (EG17-EG20)
resumen_blasto_sangre <- data.frame(
  paciente = c("D7", "D7", "D8", "D5"),
  momento = c("Dx", "2d", "Dx", "2d"),
  n_celulas = c(ncol(datos_D7_dx_sangre), ncol(datos_D7_2d), ncol(datos_D8_dx_sangre), ncol(datos_D5_2d)),
  pct_blasto_scRNAseq = c(blasto_D7_dx_sangre, blasto_D7_2d, blasto_D8_dx_sangre, blasto_D5_2d),
  blast_percentage_clinico_GEO = c(55, 10, 25, 0)
)
resumen_blasto_sangre

write.csv(resumen_blasto_sangre, file = here("resumen_blasto_por_paciente_sangre.csv"), row.names = FALSE)


#MÉDULA ÓSEA, por paciente 

datos_D5_dx_medula <- subset(dataset_integrado, subset = donor == "D5" & timepoint == "Dx" & tissue == "Bone_marrow")
datos_D5_15d <- subset(dataset_integrado, subset = donor == "D5" & timepoint == "15d")
datos_D6_dx_medula <- subset(dataset_integrado, subset = donor == "D6" & timepoint == "Dx" & tissue == "Bone_marrow")
datos_D6_15d <- subset(dataset_integrado, subset = donor == "D6" & timepoint == "15d")
datos_D7_dx_medula <- subset(dataset_integrado, subset = donor == "D7" & timepoint == "Dx" & tissue == "Bone_marrow")
datos_D8_15d_medula <- subset(dataset_integrado, subset = donor == "D8" & timepoint == "15d")

prop_D5_dx_medula <- round(prop.table(table(datos_D5_dx_medula$cell_type)) * 100, 2)
prop_D5_15d <- round(prop.table(table(datos_D5_15d$cell_type)) * 100, 2)
prop_D6_dx_medula <- round(prop.table(table(datos_D6_dx_medula$cell_type)) * 100, 2)
prop_D6_15d <- round(prop.table(table(datos_D6_15d$cell_type)) * 100, 2)
prop_D7_dx_medula <- round(prop.table(table(datos_D7_dx_medula$cell_type)) * 100, 2)
prop_D8_15d_medula <- round(prop.table(table(datos_D8_15d_medula$cell_type)) * 100, 2)

blasto_D5_dx_medula <- sum(prop_D5_dx_medula[names(prop_D5_dx_medula) %in% categorias_blasto])
blasto_D5_15d <- sum(prop_D5_15d[names(prop_D5_15d) %in% categorias_blasto])
blasto_D6_dx_medula <- sum(prop_D6_dx_medula[names(prop_D6_dx_medula) %in% categorias_blasto])
blasto_D6_15d <- sum(prop_D6_15d[names(prop_D6_15d) %in% categorias_blasto])
blasto_D7_dx_medula <- sum(prop_D7_dx_medula[names(prop_D7_dx_medula) %in% categorias_blasto])
blasto_D8_15d_medula <- sum(prop_D8_15d_medula[names(prop_D8_15d_medula) %in% categorias_blasto])

#blast_percentage_clinico_GEO (de los registros GSM (EG9, EG10, EG11, EG12, EG13, EG14))
resumen_blasto_medula <- data.frame(
  paciente = c("D5", "D5", "D6", "D6", "D7", "D8"),
  momento = c("Dx", "15d", "Dx", "15d", "Dx", "15d"),
  n_celulas = c(ncol(datos_D5_dx_medula), ncol(datos_D5_15d), ncol(datos_D6_dx_medula), ncol(datos_D6_15d), ncol(datos_D7_dx_medula), ncol(datos_D8_15d_medula)),
  pct_blasto_scRNAseq = c(blasto_D5_dx_medula, blasto_D5_15d, blasto_D6_dx_medula, blasto_D6_15d, blasto_D7_dx_medula, blasto_D8_15d_medula),
  blast_percentage_clinico_GEO = c(70, 1, 90, 1, 90, 10)
)
resumen_blasto_medula

write.csv(resumen_blasto_medula, file = here("resumen_blasto_por_paciente_medula.csv"), row.names = FALSE)



#Guardamos datos: dataset_integrado con todo el análisis completo
saveRDS(dataset_integrado, file = here("6_dataset_integrado_final.rds"))

#Información de la sesión
sessionInfo()

