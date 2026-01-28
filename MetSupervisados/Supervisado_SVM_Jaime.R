# Previamente ejecutado el preprocesado de los datos

### ------------------------------------------------------------
### ACTIVIDAD GRUPAL AIA – Método supervisado: SVM (kernel gaussiano)
### Jaime Hernán Gil
### ------------------------------------------------------------

# Comprobación de que el entorno está preparado
exists("data_final")
dim(data_final)
names(data_final)[1:5]
table(data_final$Class)

# Preparar X (predictores) e y (clases)
X <- data_final[, colnames(data_final) != "Class"]
y <- as.factor(data_final$Class)

dim(X)
length(y)

# Cargar librería (partición estratificada + entrenamiento + evaluación)
library(caret)

# Eliminar predictores con varianza cero (evita warnings al centrar/escalar)
nzv <- nearZeroVar(X, saveMetrics = TRUE)
X <- X[, !nzv$zeroVar]

# Partición train/test estratificada
set.seed(123)
train_index <- createDataPartition(y, p = 0.80, list = FALSE)

X_train <- X[train_index, ]
X_test  <- X[-train_index, ]

y_train <- y[train_index]
y_test  <- y[-train_index]

prop.table(table(y_train))
prop.table(table(y_test))

# Control del entrenamiento: validación cruzada 10-fold
ctrl <- trainControl(
  method = "cv",
  number = 10,
  classProbs = TRUE
)

# Entrenar el modelo SVM con kernel gaussiano (RBF) -> en caret: svmRadial
set.seed(123)
svm_model <- train(
  x = X_train,
  y = y_train,
  method = "svmRadial",
  trControl = ctrl,
  preProcess = c("center", "scale"),
  metric = "Accuracy"
)

svm_model
plot(svm_model)  # rendimiento vs hiperparámetro (C)

# Predicción en test y matriz de confusión
pred_test <- predict(svm_model, newdata = X_test)
cm <- confusionMatrix(pred_test, y_test)
cm

# Extracción de métricas generales
acc   <- cm$overall["Accuracy"]
kappa <- cm$overall["Kappa"]

acc
kappa

# Métricas por clase (sensibilidad y especificidad)
sens <- cm$byClass[, "Sensitivity"]
spec <- cm$byClass[, "Specificity"]

sens
spec

# F1 por clase (cálculo manual: 2 * (precision * recall) / (precision + recall))
precision <- cm$byClass[, "Pos Pred Value"]
recall    <- cm$byClass[, "Sensitivity"]

f1 <- 2 * (precision * recall) / (precision + recall)
f1[is.nan(f1)] <- NA

f1

# F1 macro-promedio (útil en multiclase)
f1_macro <- mean(f1, na.rm = TRUE)
f1_macro

### CHECKLIST FINAL 
# Accuracy (test): 0.855
# Kappa: 0.808
# F1 macro: 0.79
# Clase peor clasificada: CHC
# Confusión principal: CHC ↔ CGC
# Modelo: SVM kernel gaussiano (RBF)
# Validación: 10-fold CV + test independiente

### INTERPRETACIÓN PLOT 
# El análisis del parámetro de coste (C) mostró una mejora del rendimiento 
# del modelo hasta C = 0.5, a partir del cual la accuracy se estabiliza, 
# indicando un ajuste óptimo sin evidencia de sobreajuste. Este comportamiento 
# sugiere un modelo robusto y estable para los datos de expresión
# génica analizados.

### INTERPRETACIÓN
# El modelo SVM con kernel gaussiano alcanzó un accuracy del 85.5 % y un 
# coeficiente Kappa de 0.81, indicando un alto grado de concordancia entre 
# predicciones y clases reales. El análisis por clase mostró un rendimiento 
# excelente para la mayoría de los subtipos tumorales, aunque se observó una 
# menor sensibilidad para la clase CHC, probablemente debido al solapamiento 
# de perfiles de expresión génica con la clase CGC.
