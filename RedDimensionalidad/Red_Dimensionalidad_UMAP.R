set.seed(123)

library(tidyverse)
library(caret)
library(factoextra)

# CARGA DE LOS ARCHIVOS

setwd("C:/Users/socla/MASTER/algoritmos e IA/ACT3 grupal/")
##carga de datos de expresion genetica
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

##imputacion KNN simple
preproc <- preProcess(x, method = "knnImpute") #esto busca muestras parecidas y usa sus valores para rellenar los NA
x_imputed <- predict(preproc, x) #aqui se aplica la imputacion

##creamos un nuveo dataset con los datos imputados
data_imputed <- x_imputed
data_imputed$Class <- y

#NORMALIZACION
##centrar y escalar
preproc_scale <- preProcess(x_imputed, method = c("center", "scale"))
x_scaled <- predict(preproc_scale, x_imputed)

##dataset final
data_final <- x_scaled
data_final$Class <- y

#===============================
# UMAP - Reducción dimensional
#===============================

# Cargar librerías adicionales necesarias
library(uwot)
library(ggplot2)

# Separar datos y etiquetas para UMAP
x_umap <- data_final[, colnames(data_final) != "Class"]
labels_umap <- data_final$Class

# Aplicar UMAP con parámetros ajustados
# n_neighbors: número de vecinos cercanos (20% del total de muestras)
# n_components: dimensiones de salida (2 para visualización)
# min_dist: distancia mínima entre puntos
umap_results <- umap(x_umap, 
                     n_neighbors = 0.2 * nrow(x_umap),
                     n_components = 2, 
                     min_dist = 0.1, 
                     local_connectivity = 1,
                     ret_model = TRUE, 
                     verbose = TRUE)

# Crear dataframe con resultados de UMAP
umap_df <- data.frame(umap_results$embedding)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# Añadir las etiquetas para el gráfico
umap_df$Class <- labels_umap

# Visualizar resultados UMAP
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Class)) +
  geom_point(size = 3) +
  labs(title = "UMAP - Expresión Genética", 
       x = "UMAP1", 
       y = "UMAP2", 
       color = "Tipo") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray95"), 
        plot.title = element_text(hjust = 0.5))

## Justificación: ¿Por qué UMAP?

#La principal razón para elegir UMAP (Uniform Manifold Approximation and Projection) en estudios de expresión génica es su capacidad para preservar tanto la estructura local (qué células son similares entre sí) como la estructura global (cómo se relacionan los grupos grandes entre sí).

#Naturaleza de los datos: Los datos de expresión génica suelen tener miles de dimensiones (genes). UMAP logra "comprimir" esa información en 2D sin perder las relaciones no lineales que métodos como el PCA (lineal) suelen ignorar.

#Separabilidad para SVM: Al agrupar los datos de forma tan compacta (como se ve en el gráfico), UMAP facilita enormemente la tarea posterior de los clasificadores SVM, ya que los límites de decisión se vuelven más evidentes.

## Análisis de Ventajas y Limitaciones
  #Aspectos Positivos (Ventajas)

    #Conservación de la topología: A diferencia de t-SNE, UMAP suele mantener mejor las distancias entre clústeres lejanos. Por ejemplo, en tu gráfico, la distancia de AGH al resto tiene un significado biológico real de diferenciación.

    #Escalabilidad y Velocidad: Es computacionalmente más eficiente que otros algoritmos no lineales, lo que permite procesar grandes matrices de expresión génica rápidamente.

    #Claridad Visual: Genera agrupamientos muy definidos, lo que ayuda a identificar subtipos celulares o estados patológicos de forma intuitiva.

  #Aspectos Negativos (Limitaciones)

    #Sensibilidad a Hiperparámetros: la elección de n_neighbors (vecinos) y min_dist (distancia mínima) cambia drásticamente el resultado.

    #Un "n_neighbors" bajo se enfoca en detalles muy pequeños.

    #Un "n_neighbors" alto prioriza la visión general.

    #Interpretación de los ejes: Los ejes UMAP1 y UMAP2 no tienen una unidad física ni una magnitud biológica directa (como sí la tienen los Componentes Principales en un PCA). Son unidades abstractas.

    #Naturaleza Estocástica: Si no se fija una "semilla" (seed), el gráfico puede variar ligeramente cada vez que se ejecuta, lo que puede afectar la reproducibilidad si no se es cuidadoso.
