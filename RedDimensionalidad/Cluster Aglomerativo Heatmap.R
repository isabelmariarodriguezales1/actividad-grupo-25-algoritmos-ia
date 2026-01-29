set.seed(123)

library(tidyverse)
library(caret)
library(factoextra)

# CARGA DE LOS ARCHIVOS

setwd("C:/Users/julia/Desktop/BIOINF MASTER/Algoritmos e Inteligencia Artifical/actividad 3/")

gene_expression <- read.table("gene_expression.csv",
                              sep = ";",
                              header = FALSE)

##carga de nombres de genes
column_names <- read.table("column_names.txt",
                           header = FALSE,
                           stringsAsFactors = FALSE)
##cargar las clases
classes <- read.table("classes.csv",
                      sep = ";",
                      header = FALSE,
                      stringsAsFactors = FALSE)
colnames(classes) <- c("SampleID", "Class")


#CONSTRUIR EL DATAFRAME FINAL

##Asignar nombres de genes como nombres de columnas
colnames(gene_expression) <- column_names$V1
##Asignar IDs de muestra como nombres de filas
rownames(gene_expression) <- classes$SampleID
##Añadir la clase como última columna
data <- gene_expression %>%
  mutate(Class = as.factor(classes[,2]))


#DEPURADO DE DATOS
##calcular la proporcion de NA por gen
na_prop <- colMeans(is.na(data))
##eliminar genes con mas del 20% de NA
data <- data[, na_prop < 0.2] #(no se elimina ningun gen)

#IMPUTACION DE VALORES NA
##No haria falta ya que no hay ningun NA, pero en un caso real cabria hacerlo
###Aqui se estima un valor razonable para los NA usando la informacion disponible
###Se realiza con KNN (k-nearest neighbors)
##separar datos y clase ya que nunca imputamos la clase, solo los datos
x <- data[, colnames(data) != "Class"]
y <- data$Class

#NORMALIZACION
##centrar y escalar
preproc_scale <- preProcess(x, method = c("center", "scale"))
x_scaled <- predict(preproc_scale, x)


##dataset final
data_final <- x_scaled
data_final$Class <- y


# -------------

# Variables con desviación estandar 0
# Columnas con SD = 0 que indican que son constantes

desviaciones <- apply(x_scaled, 2, sd) # La función apply con margen 2 calcula la SD por columna.


columnas_a_mantener <- desviaciones > 0 # se mantienen las columnas que su desviación estándar es > 0
x_scaled_sd <- x_scaled # se crea una copia para no modificar el original
data.filtrada2 <- x_scaled_sd[, columnas_a_mantener]

genes_eliminados <- sum(!columnas_a_mantener)
nombres_genes_eliminados <- colnames(x_scaled_sd)[!columnas_a_mantener]



cat("Genes (columnas) eliminados (SD = 0):", genes_eliminados, "\n")
cat("Genes:",paste(nombres_genes_eliminados, collapse = ", "), "\n")


# filtro sobre el dataset original para quedarte solo con los genes útiles
x_filtrado_sd <- x[, columnas_a_mantener]

#  escalado y centrado sobre el dataset filtrado
preproc_scale_sd <- preProcess(x_filtrado_sd, method = c("center", "scale"))
X_scaled_sd <- predict(preproc_scale_sd, x_filtrado_sd)



#===============================
# CLUSTERIZACIÓN JERÁRQUICA AGLOMERATIVA: HEATMAP
#===============================

library(pheatmap)
library(RColorBrewer)

# Calcular varianza de cada gen
varianzas <- apply(x_filtrado_sd, 2, var)

# Seleccionar top genes más representativos

# Calculamos un ANOVA para cada gen respecto a la Clase
p_values_anova <- apply(x_filtrado_sd, 2, function(gen) {
  res_anova <- aov(gen ~ y)
  summary(res_anova)[[1]][["Pr(>F)"]][1]
})

# Elegimos los 20 genes con el p-valor más bajo (los más significativos)
top_genes_anova <- names(sort(p_values_anova))[1:20]
matriz_anova <- x_filtrado_sd[, top_genes_anova]


# Preparar anotaciones para mostrar las clases reales
annotation_col <- data.frame(
  Clase = y
)

rownames(annotation_col) <- rownames(matriz_anova)

# Colores para las clases
ann_colors <- list(
  Clase = setNames(
    brewer.pal(length(unique(y)), "Set1"),
    unique(y)
  )
)


# Heatmap con clustering jerárquico (método Ward)
pheatmap(t(matriz_anova),  # Transponer: genes en filas, muestras en columnas
         scale = "none",
         clustering_method = "ward.D",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Clustering Jerárquico Aglomerativo - Método Ward",
         fontsize = 10,
         border_color = NA)



#  Motivo de selección:
#  En bioinformática, el Heatmap no es solo una técnica de clusterización, sino una herramienta de visualización bidimensional. Lo hemos seleccionado porque permite realizar un análisis simultáneo: agrupa las muestras (columnas) para ver similitudes entre pacientes y, al mismo tiempo, agrupa los genes (filas) para identificar firmas moleculares.

#El Método de Ward se ha elegido porque, a diferencia de otros métodos como el single linkage (que tiende a crear grupos alargados tipo "cadena"), Ward minimiza la varianza total dentro de los clústeres. Esto genera grupos mucho más compactos y esféricos, lo cual es ideal para datos genéticos donde buscamos perfiles de expresión muy homogéneos dentro de una misma patología (BRCA, KIRC, etc.).

#Ventajas:
  
#  Jerarquía visual: Gracias al dendrograma, podemos ver no solo que dos muestras son parecidas, sino qué tan parecidas son en relación con el resto del dataset sin prefijar el número de grupos (a diferencia de K-means).

#Robustez: Como bien dicen tus apuntes, Ward es más resistente a valores atípicos (outliers) porque su métrica se basa en la suma de cuadrados, lo que suaviza el impacto de ruidos puntuales en la medición de un gen.

#Interpretación Dual: Es la única técnica que permite ver el "porqué" de un clúster: podemos señalar exactamente qué genes están "en rojo" (sobreexpresados) para justificar la agrupación de una muestra.

#Limitaciones:
  
 # Coste Computacional: Al calcular una matriz de distancias completa (N×N), si el dataset fuera de millones de muestras en lugar de 801, el ordenador se quedaría sin memoria.

#Inestabilidad (Sensibilidad): Si eliminamos algunas muestras de una clase (ej. 10 muestras de BRCA), el árbol jerárquico podría reestructurarse de forma distinta, lo que complica la replicabilidad exacta en otros estudios.

#Subjetividad: La elección de la distancia (Euclídea) y el linkage (Ward) es una decisión del investigador; otro analista podría usar distancia Manhattan y obtener grupos diferentes.

#3. ¿Son estos clústeres los mejores posibles?
  
#  No se puede afirmar con certeza absoluta que sean los "mejores posibles". En el aprendizaje no supervisado no existe una "verdad absoluta" (ground truth) de la misma forma que en el supervisado, ya que el algoritmo no intenta "predecir", sino "encontrar estructura".

  
 # Dependencia de la Métrica: La clusterización es esclava de la métrica de disimilitud. Hemos usado la distancia Euclídea, que asume relaciones lineales. Si los datos tuvieran relaciones no lineales complejas, una métrica de correlación de Pearson podría haber sido "mejor".

# Óptimo Local vs. Global: Los métodos jerárquicos toman decisiones en cada paso (unir las dos muestras más cercanas) que no pueden deshacerse después. Esto significa que el algoritmo es "codicioso" (greedy) y puede quedar atrapado en una solución buena, pero no necesariamente la óptima global para todo el conjunto.

# Validación: Para decir que son los mejores, necesitaríamos calcular métricas adicionales como el Coeficiente de Correlación Cofenética (que mide qué tan bien el dendrograma preserva las distancias originales) o el Análisis de Silueta.


