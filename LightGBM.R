library(data.table)
library(lightgbm)
library(caret)
library(dplyr)
library(hydroGOF)
library(Matrix)
library(stringr)

# Configuración
set.seed(1234)
directory = "C:/2. ARTICULOS/Relleno_datos/BD_Procesada/BD_HORARIA/DATA_COMPLETE"
output_dir = file.path(directory, "IMPUTACION_RESULTADOS_LightGBM_ValidacionCorrecta")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

VARIABLES_PRIORITARIAS = c("VelocidadViento", "RS_W_Avg", "BP_mbar_Avg", "HumAire_Avg", "TempAire_Avg")

# Función para calcular métricas GOF
calculate_gof_metrics <- function(observed, predicted) {
  valid_indices = complete.cases(observed, predicted)
  if (sum(valid_indices) < 10) return(NULL)
  
  obs_clean = observed[valid_indices]
  pred_clean = predicted[valid_indices]
  
  tryCatch({
    gof_results = gof(sim = pred_clean, obs = obs_clean, na.rm = TRUE)
    if (is.null(names(gof_results))) {
      names(gof_results) = c("ME", "MAE", "MSE", "RMSE", "NRMSE", "PBIAS", "RSR", 
                             "rSD", "NSE", "mNSE", "rNSE", "d", "md", "rd", "cp", 
                             "r", "R2", "bR2", "KGE", "VE")
    }
    return(gof_results)
  }, error = function(e) {
    return(c(
      ME = mean(pred_clean - obs_clean, na.rm = TRUE),
      MAE = mean(abs(pred_clean - obs_clean), na.rm = TRUE),
      RMSE = sqrt(mean((pred_clean - obs_clean)^2, na.rm = TRUE)),
      R2 = cor(obs_clean, pred_clean, use = "complete.obs")^2,
      NSE = 1 - sum((obs_clean - pred_clean)^2) / sum((obs_clean - mean(obs_clean))^2),
      PBIAS = 100 * sum(pred_clean - obs_clean) / sum(obs_clean)
    ))
  })
}

# Función para cargar datos
load_data <- function(directory) {
  archivos = list.files(directory, pattern = "*.csv", full.names = TRUE)
  datos_list = list()
  for (archivo in archivos) {
    nombre_estacion = str_extract(basename(archivo), "^[^.]+")
    datos_list[[nombre_estacion]] = fread(archivo)
  }
  return(datos_list)
}

# Función para encontrar covariables correlacionadas
find_covariables <- function(datos_list, target_var, target_station, threshold = 0.5) {
  covariables = c()
  
  for (station_name in names(datos_list)) {
    if (station_name == target_station) next
    
    df = datos_list[[station_name]]
    numeric_vars = names(df)[sapply(df, is.numeric)]
    
    for (var in numeric_vars) {
      covar_name = paste0(var, "_", station_name)
      df_target = datos_list[[target_station]]
      
      if (target_var %in% names(df_target) && var %in% names(df)) {
        df_temp1 = df_target[, .(TIMESTAMP, target_value = get(target_var))]
        df_temp2 = df[, .(TIMESTAMP, predictor_value = get(var))]
        temp_merge = merge(df_temp1, df_temp2, by = "TIMESTAMP", all = FALSE)
        
        if (nrow(temp_merge) > 30) {
          temp_merge = temp_merge[complete.cases(temp_merge)]
          if (nrow(temp_merge) > 30) {
            cor_value = cor(temp_merge$target_value, temp_merge$predictor_value, use = "complete.obs")
            if (!is.na(cor_value) && abs(cor_value) >= threshold) {
              covariables = c(covariables, covar_name)
            }
          }
        }
      }
    }
  }
  
  # Agregar variables de la misma estación si es necesario
  if (length(covariables) < 3) {
    df_target = datos_list[[target_station]]
    same_station_vars = names(df_target)[sapply(df_target, is.numeric)]
    same_station_vars = same_station_vars[same_station_vars != target_var]
    
    for (var in same_station_vars) {
      if (var %in% names(df_target)) {
        temp_cor = cor(df_target[[target_var]], df_target[[var]], use = "complete.obs")
        if (!is.na(temp_cor) && abs(temp_cor) >= 0.3) {
          covar_name = paste0(var, "_", target_station)
          covariables = c(covariables, covar_name)
        }
      }
    }
  }
  
  return(covariables)
}

# Función para crear dataset combinado
create_combined_dataset <- function(datos_list, target_station, target_var, selected_covariables) {
  df_combined = copy(datos_list[[target_station]])
  
  for (covar in selected_covariables) {
    covar_split = strsplit(covar, "_")[[1]]
    if (length(covar_split) >= 2) {
      station_name = covar_split[length(covar_split)]
      var_name = paste(covar_split[1:(length(covar_split)-1)], collapse = "_")
      
      if (station_name %in% names(datos_list) && var_name %in% names(datos_list[[station_name]])) {
        df_temp = datos_list[[station_name]][, .(TIMESTAMP, temp_var = get(var_name))]
        setnames(df_temp, "temp_var", covar)
        df_combined = merge(df_combined, df_temp, by = "TIMESTAMP", all.x = TRUE)
      }
    }
  }
  
  return(df_combined)
}

# Función principal de imputación con validación correcta
impute_with_proper_validation <- function(df_combined, target_var, station_name, selected_covariables) {
  # Crear características temporales
  df_combined[, `:=`(
    Year = as.numeric(format(TIMESTAMP, "%Y")),
    Month = as.numeric(format(TIMESTAMP, "%m")),
    Day = as.numeric(format(TIMESTAMP, "%d")),
    Hour = as.numeric(format(TIMESTAMP, "%H")),
    DayOfYear = as.numeric(format(TIMESTAMP, "%j"))
  )]
  
  # Identificar índices con NA
  na_indices = which(is.na(df_combined[[target_var]]))
  df_train = df_combined[!is.na(df_combined[[target_var]]), ]
  
  # Seleccionar predictores
  temporal_vars = c("Year", "Month", "Day", "Hour", "DayOfYear")
  predictor_vars = c(selected_covariables, temporal_vars)
  predictor_vars = predictor_vars[predictor_vars %in% names(df_train)]
  
  # Asegurar predictores mínimos
  if (length(predictor_vars) < 2) {
    all_numeric = names(df_train)[sapply(df_train, is.numeric)]
    all_numeric = all_numeric[!all_numeric %in% c(target_var, "TIMESTAMP")]
    predictor_vars = unique(c(predictor_vars, all_numeric, temporal_vars))
    predictor_vars = predictor_vars[predictor_vars %in% names(df_train)]
  }
  
  # Preparar datos
  X = df_train[, ..predictor_vars]
  y = df_train[[target_var]]
  
  # Imputar NA en predictores
  for (col in names(X)) {
    if (is.numeric(X[[col]])) {
      median_val = median(X[[col]], na.rm = TRUE)
      if (is.na(median_val)) median_val = 0
      X[is.na(get(col)), (col) := median_val]
    }
  }
  
  # Remover filas con NA
  complete_rows = complete.cases(X)
  X = X[complete_rows, ]
  y = y[complete_rows]
  
  if (length(y) < 50) {
    cat("Pocos datos para", station_name, target_var, ". Usando imputación por mediana.\n")
    median_val = median(df_combined[[target_var]], na.rm = TRUE)
    if (is.na(median_val)) median_val = 0
    df_combined[na_indices, (target_var) := median_val]
    return(list(df_imputed = df_combined, cv_metrics = NULL, test_metrics = NULL, covariables_used = predictor_vars))
  }
  
  # VALIDACIÓN CORREGIDA: División 80-20 + Validación cruzada robusta
  set.seed(1234)
  
  # División inicial 80-20 (semilla fija para reproducibilidad)
  train_size = floor(0.8 * length(y))
  train_indices = sample(1:length(y), train_size)
  test_indices = setdiff(1:length(y), train_indices)
  
  X_train = X[train_indices, ]
  y_train = y[train_indices]
  X_test = X[test_indices, ]
  y_test = y[test_indices]
  
  # Validación cruzada 5-fold SOLO en entrenamiento (SIN semilla fija)
  folds = createFolds(y_train, k = 5, list = TRUE, returnTrain = TRUE)  # returnTrain = TRUE
  cv_predictions = c()
  cv_actuals = c()
  
  # Parámetros LightGBM más conservadores
  lgb_params = list(
    objective = "regression",
    metric = "rmse",
    num_leaves = 15,           # Reducido para evitar overfitting
    learning_rate = 0.01,      # Más conservador
    feature_fraction = 0.8,    # Más regularización
    bagging_fraction = 0.7,    # Más regularización
    bagging_freq = 5,
    min_data_in_leaf = 50,     # Incrementado para mayor robustez
    lambda_l1 = 0.1,           # Regularización L1
    lambda_l2 = 0.1,           # Regularización L2
    verbose = -1
  )
  
  # Ejecutar validación cruzada
  for (i in 1:5) {
    fold_train_indices = folds[[i]]  # Ahora son índices de entrenamiento
    fold_test_indices = setdiff(1:nrow(X_train), fold_train_indices)
    
    X_fold_train = as.matrix(X_train[fold_train_indices, ])
    y_fold_train = y_train[fold_train_indices]
    X_fold_test = as.matrix(X_train[fold_test_indices, ])
    y_fold_test = y_train[fold_test_indices]
    
    # Crear datasets con validación
    lgb_train = lgb.Dataset(data = X_fold_train, label = y_fold_train)
    lgb_valid = lgb.Dataset(data = X_fold_test, label = y_fold_test, reference = lgb_train)
    
    # Entrenar con early stopping
    lgb_model = lgb.train(
      params = lgb_params, 
      data = lgb_train, 
      nrounds = 1000,
      valids = list(valid = lgb_valid),
      early_stopping_rounds = 50,
      verbose = -1
    )
    
    fold_predictions = predict(lgb_model, X_fold_test)
    cv_predictions = c(cv_predictions, fold_predictions)
    cv_actuals = c(cv_actuals, y_fold_test)
  }
  
  # Entrenar modelo final con early stopping
  X_train_matrix = as.matrix(X_train)
  lgb_train_final = lgb.Dataset(data = X_train_matrix, label = y_train)
  
  # Usar una porción del entrenamiento como validación para early stopping
  set.seed(5678)  # Semilla diferente para validación interna
  val_indices = sample(1:nrow(X_train), floor(0.15 * nrow(X_train)))
  train_final_indices = setdiff(1:nrow(X_train), val_indices)
  
  X_train_final = as.matrix(X_train[train_final_indices, ])
  y_train_final = y_train[train_final_indices]
  X_val_final = as.matrix(X_train[val_indices, ])
  y_val_final = y_train[val_indices]
  
  lgb_train_final = lgb.Dataset(data = X_train_final, label = y_train_final)
  lgb_val_final = lgb.Dataset(data = X_val_final, label = y_val_final, reference = lgb_train_final)
  
  lgb_model_final = lgb.train(
    params = lgb_params, 
    data = lgb_train_final, 
    nrounds = 1000,
    valids = list(valid = lgb_val_final),
    early_stopping_rounds = 50,
    verbose = -1
  )
  
  # Evaluación final en conjunto de prueba
  X_test_matrix = as.matrix(X_test)
  test_predictions = predict(lgb_model_final, X_test_matrix)
  
  # Calcular métricas
  cv_metrics = calculate_gof_metrics(cv_actuals, cv_predictions)
  test_metrics = calculate_gof_metrics(y_test, test_predictions)
  
  # Mostrar resultados
  cat("\n=== VALIDACIÓN:", station_name, "-", target_var, "===\n")
  cat("Datos - Total:", length(y), "| Entrenamiento:", length(y_train), "| Prueba:", length(y_test), "\n")
  
  if (!is.null(cv_metrics)) {
    cat("CV - R²:", round(cv_metrics["R2"], 3), "| NSE:", round(cv_metrics["NSE"], 3), 
        "| RMSE:", round(cv_metrics["RMSE"], 3), "\n")
  }
  
  if (!is.null(test_metrics)) {
    cat("TEST - R²:", round(test_metrics["R2"], 3), "| NSE:", round(test_metrics["NSE"], 3), 
        "| RMSE:", round(test_metrics["RMSE"], 3), "\n")
  }
  
  # Imputar valores faltantes (sin cambios)
  if (length(na_indices) > 0) {
    X_impute = df_combined[na_indices, ..predictor_vars]
    
    for (col in names(X_impute)) {
      if (is.numeric(X_impute[[col]])) {
        median_val = median(df_combined[[col]], na.rm = TRUE)
        if (is.na(median_val)) median_val = 0
        X_impute[is.na(get(col)), (col) := median_val]
      }
    }
    
    X_impute_matrix = as.matrix(X_impute)
    imputados = predict(lgb_model_final, X_impute_matrix)
    df_combined[na_indices, (target_var) := imputados]
    
    cat("Imputados:", length(imputados), "(", round(length(imputados) / nrow(df_combined) * 100, 2), "%)\n")
  }
  
  # Verificar NA restantes
  remaining_na = sum(is.na(df_combined[[target_var]]))
  if (remaining_na > 0) {
    median_val = median(df_combined[[target_var]], na.rm = TRUE)
    if (is.na(median_val)) median_val = 0
    df_combined[is.na(get(target_var)), (target_var) := median_val]
  }
  
  return(list(
    df_imputed = df_combined,
    cv_metrics = cv_metrics,
    test_metrics = test_metrics,
    covariables_used = predictor_vars
  ))
}

# Función principal de procesamiento
process_all_stations <- function(directory, correlation_threshold = 0.5) {
  datos_list = load_data(directory)
  estaciones = names(datos_list)
  
  all_results = list()
  all_metrics = data.table()
  
  for (estacion in estaciones) {
    cat("\n=== PROCESANDO ESTACIÓN:", estacion, "===\n")
    
    df_station_imputed = copy(datos_list[[estacion]])
    station_metrics = data.table()
    
    # Procesar SOLO variables prioritarias
    for (variable in VARIABLES_PRIORITARIAS) {
      if (!(variable %in% names(df_station_imputed))) next
      
      na_count = sum(is.na(df_station_imputed[[variable]]))
      if (na_count == 0) {
        cat("Variable", variable, "sin valores faltantes. Saltando...\n")
        next
      }
      
      # Encontrar covariables
      covariables = find_covariables(datos_list, variable, estacion, correlation_threshold)
      
      if (length(covariables) == 0) {
        cat("No se encontraron covariables para", variable, ". Usando mediana.\n")
        median_val = median(df_station_imputed[[variable]], na.rm = TRUE)
        if (is.na(median_val)) median_val = 0
        df_station_imputed[is.na(get(variable)), (variable) := median_val]
        next
      }
      
      # Crear dataset y realizar imputación
      df_combined = create_combined_dataset(datos_list, estacion, variable, covariables)
      
      tryCatch({
        result = impute_with_proper_validation(df_combined, variable, estacion, covariables)
        
        # Actualizar dataframe
        df_station_imputed[[variable]] = result$df_imputed[[variable]]
        all_results[[paste(estacion, variable, sep = "_")]] = result
        
        # Procesar métricas
        metrics_types = c("cv_metrics", "test_metrics")
        for (metric_type in metrics_types) {
          if (!is.null(result[[metric_type]])) {
            gof_vector = if (is.matrix(result[[metric_type]])) as.vector(result[[metric_type]]) else result[[metric_type]]
            if (is.matrix(result[[metric_type]])) names(gof_vector) = rownames(result[[metric_type]])
            
            metric_names = names(gof_vector)
            metric_values = as.numeric(gof_vector)
            
            if (length(metric_names) > 0) {
              metrics_df = data.table(
                Estacion = rep(estacion, length(metric_names)),
                Variable = rep(variable, length(metric_names)),
                Tipo_Metrica = rep(ifelse(metric_type == "cv_metrics", "Validacion_Cruzada", "Prueba_Final"), length(metric_names)),
                Metrica = metric_names,
                Valor = metric_values,
                Covariables_Utilizadas = rep(paste(result$covariables_used, collapse = "; "), length(metric_names))
              )
              station_metrics = rbind(station_metrics, metrics_df)
            }
          }
        }
        
      }, error = function(e) {
        cat("Error procesando", estacion, variable, ":", e$message, "\n")
        median_val = median(df_station_imputed[[variable]], na.rm = TRUE)
        if (is.na(median_val)) median_val = 0
        df_station_imputed[is.na(get(variable)), (variable) := median_val]
      })
    }
    
    # Guardar resultados
    fwrite(df_station_imputed, file.path(output_dir, paste0(estacion, "_imputed.csv")))
    if (nrow(station_metrics) > 0) {
      fwrite(station_metrics, file.path(output_dir, paste0(estacion, "_metrics.csv")))
      all_metrics = rbind(all_metrics, station_metrics)
    }
  }
  
  # Guardar métricas consolidadas
  if (nrow(all_metrics) > 0) {
    fwrite(all_metrics, file.path(output_dir, "all_metrics_consolidated.csv"))
  }
  
  # Generar informe
  generate_report(all_results, all_metrics)
  
  return(list(results = all_results, metrics = all_metrics))
}

# Función para generar informe
generate_report <- function(all_results, all_metrics) {
  report_file = file.path(output_dir, "informe_imputacion_validacion_correcta.txt")
  
  sink(report_file)
  cat("=== INFORME DE IMPUTACIÓN CON VALIDACIÓN CORRECTA ===\n")
  cat("Fecha:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Método: División 80-20 + Validación Cruzada 5-Fold\n")
  cat("Variables procesadas:", paste(VARIABLES_PRIORITARIAS, collapse = ", "), "\n\n")
  
  if (nrow(all_metrics) > 0) {
    cat("RESUMEN GENERAL\n")
    cat("===============\n")
    cat("Estaciones procesadas:", length(unique(all_metrics$Estacion)), "\n")
    cat("Variables procesadas:", length(unique(all_metrics$Variable)), "\n\n")
    
    for (estacion in unique(all_metrics$Estacion)) {
      cat("ESTACIÓN:", estacion, "\n")
      cat("==================\n")
      
      metrics_estacion = all_metrics[Estacion == estacion]
      
      for (variable in unique(metrics_estacion$Variable)) {
        cat("Variable:", variable, "\n")
        metrics_var = metrics_estacion[Variable == variable]
        
        # Mostrar métricas clave por tipo
        for (tipo in unique(metrics_var$Tipo_Metrica)) {
          cat("  ", tipo, ":\n")
          tipo_metrics = metrics_var[Tipo_Metrica == tipo]
          
          key_metrics = c("R2", "NSE", "RMSE", "MAE", "PBIAS")
          for (key_metric in key_metrics) {
            metric_row = tipo_metrics[Metrica == key_metric]
            if (nrow(metric_row) > 0) {
              cat("    ", key_metric, ":", round(metric_row$Valor[1], 4), "\n")
            }
          }
        }
        
        if (nrow(metrics_var) > 0) {
          cat("  Covariables:", unique(metrics_var$Covariables_Utilizadas)[1], "\n")
        }
        cat("\n")
      }
    }
  }
  
  sink()
  cat("Informe generado en:", report_file, "\n")
}

# EJECUCIÓN PRINCIPAL
cat("=== INICIANDO IMPUTACIÓN CON VALIDACIÓN CORRECTA ===\n")
cat("Solo procesando variables prioritarias:", paste(VARIABLES_PRIORITARIAS, collapse = ", "), "\n")

resultados_finales = process_all_stations(directory)

cat("\n=== PROCESO COMPLETADO ===\n")
cat("Archivos generados:\n")
cat("- Datos imputados: *_imputed.csv\n")
cat("- Métricas por estación: *_metrics.csv\n")
cat("- Métricas consolidadas: all_metrics_consolidated.csv\n")
cat("- Informe: informe_imputacion_validacion_correcta.txt\n")
