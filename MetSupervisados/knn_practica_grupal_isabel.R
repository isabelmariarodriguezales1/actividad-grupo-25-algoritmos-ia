set.seed(123)

library(tidyverse)
library(caret)
library(factoextra)


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

#IMPLEMENTACION DE KNN

# Se construye el dataset final para KNN:
# - Variables predictoras: genes escalados y filtrados (X_scaled_sd)
# - Variable respuesta: Class

data_knn <- as.data.frame(X_scaled_sd)
data_knn$Class <- y


# -------------------------------
# División en entrenamiento y test
# -------------------------------

set.seed(123)

# Se mantiene la proporción de clases
train_index <- createDataPartition(data_knn$Class,
                                   p = 0.8,
                                   list = FALSE)

train_data <- data_knn[train_index, ]
test_data  <- data_knn[-train_index, ]


# -------------------------------
# Entrenamiento del modelo KNN
# -------------------------------

# Se utiliza validación cruzada (10-fold CV)
# para seleccionar el valor óptimo de k

control <- trainControl(
  method = "cv",
  number = 10
)

# Se prueban distintos valores de k (solo impares)
grid_k <- expand.grid(
  k = seq(3, 21, by = 2)
)

set.seed(123)

knn_model <- train(
  Class ~ .,
  data = train_data,
  method = "knn",
  trControl = control,
  tuneGrid = grid_k
)

# Mostrar el mejor valor de k
knn_model

#El mejor valor de K que se ha obtenido ha sido o bien k = 9, k = 11, k = 13, pero caret elige 13 (primer máximo estable)

# -------------------------------
# Evaluación del modelo
# -------------------------------

# Predicciones sobre el conjunto de test
pred_knn <- predict(knn_model, newdata = test_data)

# Matriz de confusión y métricas de evaluación
cm <- confusionMatrix(pred_knn, test_data$Class)
cm

cm_df <- as.data.frame(cm$table)

library(ggplot2)

ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Matriz de confusión del clasificador KNN",
       x = "Clase real",
       y = "Clase predicha") +
  theme_minimal()

#INTERPRETACION DE LOS RESULTADOS - 
# El clasificador KNN obtuvo un rendimiento muy elevado tanto en validación cruzada como en el conjunto de test. 
#El valor óptimo de k fue 13, seleccionado mediante validación cruzada 10-fold. En el conjunto de test, el modelo alcanzó una 
#precisión del 100%, con valores de sensibilidad, especificidad y score F1 iguales a 1 para todas las clases.

#Este resultado sugiere una alta separabilidad entre las clases, probablemente debida a diferencias biológicas claras en los 
#perfiles de expresión génica. No obstante, aunque el rendimiento es excelente, es importante destacar que KNN puede verse 
#favorecido en datasets donde las clases están bien definidas y correctamente escaladas.

#A pesar de los resultados obtenidos, cabe destacar que en escenarios con mayor solapamiento entre clases o con un mayor número 
#de muestras, el rendimiento del modelo podría verse reducido, por lo que sería recomendable contrastar estos resultados con otros 
#clasificadores.