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

