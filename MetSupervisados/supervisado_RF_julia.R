set.seed(123)

library(tidyverse)
library(caret)
library(factoextra)

# CARGA DE LOS ARCHIVOS
#setwd("C:/Users/socla/MASTER/algoritmos e IA/ACT3 grupal/")
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

#APRENDIZAJE NO SUPERVISADO RANDOM FOREST
library(randomForest)

##division de datos (0.7 train y 0.3 test)
trainIndex <- createDataPartition(data_final$Class, p = 0.7, lis = FALSE)
train_data <- data_final[trainIndex, ]
test_data <- data_final[-trainIndex, ]

cat("Dimension de train_data:", dim(train_data), "\n")
cat("Dimension de test_data:", dim(test_data), "\n")
cat("Dimension de clases en train:\n")
print(table(train_data$Class))
cat("\nDistribucion de clases en test:\n")
print(table(test_data$Class))

#ENTRENAMIENTO DEL MODELO
##entrenar Random Forest con validacion cruzada
rf_model <- train(
  Class ~ .,
  data = train_data,
  method = "rf",
  trControl = trainControl(
    method = "cv", #validacion cruzada
    number = 5, #5 folds
    savePredictions = TRUE,
    classProbs = TRUE #para calcular curvas ROC
  ),
  ntree = 75, #numero de arboles
  importance = TRUE #calcular importancia de variables
)

print(rf_model)
cat("\nMejor valor de mtry:", rf_model$bestTune$mtry, "\n")

#PREDICCIONES
##predicciones en test_data
predictions <- predict(rf_model, test_data)
predictions
predictions_prob <- predict(rf_model, test_data, type = "prob")
predictions_prob

#METRICAS DE EVALUACION
##matriz de confusion
conf_matrix <- confusionMatrix(predictions, test_data$Class)
print(conf_matrix)

##extraer metricas clave
cat("Accuracy:", round(conf_matrix$overall['Accuracy'], 4), "\n")
cat("\nMétricas por clase:\n")
print(conf_matrix$byClass[, c('Sensitivity', 'Specificity', 'Precision', 'F1')])

#metricas promedio
cat("\nMétricas promedio (macro):\n")
cat("Sensitivity promedio:", round(mean(conf_matrix$byClass[,'Sensitivity'], na.rm=TRUE), 4), "\n")
cat("Specificity promedio:", round(mean(conf_matrix$byClass[,'Specificity'], na.rm=TRUE), 4), "\n")
cat("F1-Score promedio:", round(mean(conf_matrix$byClass[,'F1'], na.rm=TRUE), 4), "\n")

#IMPORTANCIA DE VARIABLES (GENES)
importance <- varImp(rf_model, scale = FALSE)
cat("\nTop 20 genes mas importances:\n")
print(head(importance$importance, 20))

##visualizacion de importancia
pdf("random_forest_importance.pdf", width = 10, height = 8)
plot(importance, top = 20, main = "Top 20 genes mas importances")
dev.off()
cat("\nGrafico guardado en: random_forest_importance.pdf\n")


#ANALISIS DE ERRORES
##Identificar muestras mal clasificadas
errors <- which(predictions != test_data$Class)
cat("Numero de muestras mal clasificadas:", length(errors), "\n")
cat("Porcentaje de error:", round(length(errors)/nrow(test_data)*100, 2), "%\n")

if(length(errors) > 0 && length(errors) <= 10) {
  cat("\nMuestras mal clasificadas:\n")
  error_df <- data.frame(
    Muestra = rownames(test_data)[errors],
    Clase_Real = test_data$Class[errors],
    Prediccion = predictions[errors]
  )
  print(error_df)
}
