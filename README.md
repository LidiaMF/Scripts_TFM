# TFM:Análisis integrado de la expresión génica y de proteínas de superficie a nivel unicelular en muestras de médula ósea y sangre periférica de pacientes con LLA-B
Autora: Lidia Martínez Fernández de Sevilla
## Descripción
Se incluyen los scripts utilizados para el desarrollo y análisis de datos del TFM titulado: Análisis integrado de la expresión génica y de proteínas de superficie a nivel unicelular en muestras de médula ósea y sangre periférica de pacientes con LLA-B. 
Se ha utilizado para el paquete Seurat en R.

## Contenido
###1_Carga_datos_QC
  - Se cargan matrices de conteo de la serie GSE230295) descargadas de GEO
  - Carga y preprocesamiento de los datos. Control de calidad. Eliminación de dobletes
    
###2_Demultiplexado
  - Asignación muestra y paciente
  - Incorporación del tejido y momento clínico
    
###3_Integración_clustering
  - Normalización, selección de genes variables, escalado y PCA por librería
  - Integración de las 6 librerías (RPCA) para corregir el efecto de lote
  - Clustering
  - Proyección UMAP
    
###4_Anotación
  -  Identificación de genes marcadores por clúster (FindAllMarkers)
  -  Construcción de un panel de marcadores de ARN combinado con la expresión de proteína de superficie
  -  Anotación manual la identidad de cada clúster
    
###5_Anotación_singleR
  -Anotacion celular automática con SingleR, usando dos referencias de celldex: NovershternHematopoieticData y HumanPrimaryCellAtlasData  

###6_Comparaciones
  - Comparación de la composición de tipos celulares entre médula ósea y sangre periférica
  - Comparación entre diagnóstico y post-tratamiento
  - Comparación intraindividuo
  - Comparación % blastos con metadatos de GEO
