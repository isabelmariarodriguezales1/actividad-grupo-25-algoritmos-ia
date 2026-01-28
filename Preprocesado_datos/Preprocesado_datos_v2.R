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


