
library(terra)
library(sf)
library(maxnet)
library(dplyr)
library(caret)
library(Metrics)
library(ggplot2)
library(gridExtra)
library(blockCV)  
library(raster)
library(dismo)
library(rJava)
library(pROC)
library(ENMeval)
library(maxnet)
library(spatialsample)


###-------------------- Functions -----------------------###

# function to extract the env data based on the selected nlm layers found in the species
get_env_data <- function(species, nlms, used_vars){
  # map variables to nlms
  mapped_vars <- paste0("NLM", regmatches(used_vars, regexpr("[0-9]+$", used_vars)))
  
  # extract layers/data
  env_stack <- nlms[[mapped_vars]]
  pa_df <- species[[3]]$sample.points
  env_data <- terra::extract(env_stack, pa_df[, c("x", "y")])
  data <- cbind(env_data[, -1], response = pa_df$Observed, x = pa_df$x, y = pa_df$y)  
  data <- na.omit(data)
  return(data)
}


# run model for different CV methods and species
run_species_cv <- function(species_data, predictors_cols, response_col, cv_folds, species_name) {
  
  # prepare data
  predictors <- st_drop_geometry(species_data[, predictors_cols])
  response <- species_data[[response_col]]
  
  # prepare for the results
  all_predictions <- numeric(nrow(species_data))
  all_observed <- numeric(nrow(species_data))
  fold_metrics <- data.frame(
    fold = 1:length(cv_folds$folds_list),
    RMSE = NA,
    MAE = NA,
    Pearson = NA,
    AUC = NA
  )
  
  # Loop over all folds
  for(i in 1:length(cv_folds$folds_list)) {
    
    # train and test indices
    train_idx <- cv_folds$folds_list[[i]][[1]]
    test_idx <- cv_folds$folds_list[[i]][[2]]
    
    # train and test data 
    train_predictors <- predictors[train_idx, ]
    train_response <- response[train_idx]
    
    test_predictors <- predictors[test_idx, ]
    test_response <- response[test_idx]
    
    # fit train data to maxnet model
    model <- maxnet(train_response, train_predictors)
    
    # predict data
    test_predictions <- predict(model, test_predictors, type = "logistic")
    
    # write results
    all_predictions[test_idx] <- test_predictions
    all_observed[test_idx] <- test_response
    
    # calculate metrics (rmse, mae, pearson and auc) for each fold
    fold_metrics$RMSE[i] <- sqrt(mean((test_response - test_predictions)^2))
    fold_metrics$MAE[i] <- mean(abs(test_response - test_predictions))
    fold_metrics$Pearson[i] <- cor(test_response, test_predictions)
    
    if(length(unique(test_response)) > 1) {
      fold_metrics$AUC[i] <- as.numeric(auc(test_response, test_predictions))
    }
  }
  
  # calculate metrics (rmse, mae, pearson and auc) for all folds combined
  cv_metrics <- data.frame(
    species = species_name,
    RMSE = sqrt(mean((all_observed - all_predictions)^2)),
    MAE = mean(abs(all_observed - all_predictions)),
    Pearson = cor(all_observed, all_predictions),
    AUC = as.numeric(auc(all_observed, all_predictions))
  )
  
  # calculate mean of all metrics (rmse, mae, pearson and auc) for each fold
  mean_fold_metrics <- data.frame(
    species = species_name,
    RMSE_mean = mean(fold_metrics$RMSE, na.rm = TRUE),
    MAE_mean = mean(fold_metrics$MAE, na.rm = TRUE),
    Pearson_mean = mean(fold_metrics$Pearson, na.rm = TRUE),
    AUC_mean = mean(fold_metrics$AUC, na.rm = TRUE)
  )
  
  cat("Species", species_name, "completed\n")
  
  return(list(
    cv_metrics = cv_metrics,
    mean_fold_metrics = mean_fold_metrics,
    fold_metrics = fold_metrics,
    predictions = all_predictions,
    observed = all_observed
  ))
}


# define all 5 cross validation methods 
get_cv_methods <- function(suitability) {
  list(
    # random k-cross validation (k=10 folds)
    krandom = function(data) {
      n <- nrow(data)
      k = 10
      fold_assignment <- sample(rep(1:k, length.out = n))
      folds <- lapply(1:k, function(i) list(
        which(fold_assignment != i),
        which(fold_assignment == i)
      ))
      structure(list(folds = folds, folds_list = folds, k = k, records = list(method = "random", k = k, cells = n)), class = c("cv_random", "list"))
    },
    
    # Block CV (separates samples in hexagonal "blocks", k=10 folds, at size 50, assignment of blocks into fold at random)
    block = function(data) cv_spatial(x = data,
                                      column = "response",
                                      r = suitability,
                                      k = 10,
                                      size = 50,
                                      selection = "random",
                                      iteration = 50,
                                      biomod2 = FALSE),
    
    # Spacial Cross Validation (clustering by euclidean distance, 10 folds)
    scv = function(data) cv_cluster(x = data,
                                    column = "response",
                                    k = 10),
    
    # Ecological Cross Validation (similar to scv, but clustering is done based on the suitablity raster, "object of covariates to identify environmental groups")
    ecv = function(data) cv_cluster(x = data,
                                    column = "response",
                                    r = suitability,
                                    k = 10,
                                    scale = TRUE),
    
    # Latitudinal Cross Validation (separated into k (10) "stripes", and always select one as a fold)
    latitudinal = function(data) {
      coords <- st_coordinates(data)
      latitude <- coords[, 2]
      
      k_folds <- 10
      lat_range <- range(latitude)
      lat_breaks <- seq(lat_range[1], lat_range[2], length.out = k_folds + 1)
      fold_assignment <- cut(latitude, breaks = lat_breaks, labels = 1:k_folds, include.lowest = TRUE)
      
      folds <- lapply(1:k_folds, function(i) list(
        which(as.numeric(fold_assignment) != i),
        which(as.numeric(fold_assignment) == i)
      ))
      
      result <- list()
      result$folds <- folds
      result$k <- k_folds
      result$records <- list(
        method = "latitudinal",
        k = k_folds,
        cells = nrow(data)
      )
      result$folds_list <- folds
      
      class(result) <- c("cv_latitudinal", "list")
      
      return(result)
    }
  )
}


# function to calculate map accuracy against ground truth 
evaluate_map_accuracy <- function(species_obj, nlm_stack, used_vars, species_name) {
  
  cat("Evaluiere Species:", species_name, "\n")
  cat("Verwendete Variablen:", paste(used_vars, collapse=", "), "\n")
  
  # extract env data
  env_data <- tryCatch({
    get_env_data(species_obj, nlm_stack, used_vars)
  }, error = function(e) {
    cat("Fehler bei get_env_data:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(env_data)) {
    return(NULL)
  }
  
  cat("Environmental data extrahiert:", nrow(env_data), "Punkte", "\n")
  
  # BEIDE Ground Truth Raster laden
  ground_truth_pa <- unwrap(species_obj[[2]]$pa.raster)        # Für AUC
  ground_truth_suit <- unwrap(species_obj[[1]]$suitab.raster)  # Für RMSE, MAE, Pearson
  
  # Training data as sf object
  training_data <- sf::st_as_sf(env_data, coords = c("x", "y"), crs = 32633)
  
  # prepare maxent model
  predictors <- st_drop_geometry(training_data[, used_vars])
  response <- training_data$response
  
  cat("Predictor-Spalten:", paste(colnames(predictors), collapse=", "), "\n")
  cat("Erwartete Variablen:", paste(used_vars, collapse=", "), "\n")
  
  missing_vars <- used_vars[!used_vars %in% colnames(predictors)]
  if (length(missing_vars) > 0) {
    cat("Fehlende Predictor-Spalten:", paste(missing_vars, collapse=", "), "\n")
    return(NULL)
  }
  
  # fit model
  maxent_model <- tryCatch({
    maxnet(response, predictors)
  }, error = function(e) {
    cat("Fehler bei MaxEnt Modell:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(maxent_model)) {
    return(NULL)
  }
  
  cat("MaxEnt Modell erstellt", "\n")
  
  # predict the full raster
  prediction_raster <- tryCatch({
    predict(nlm_stack[[used_vars]], maxent_model, type = "logistic", na.rm=TRUE)
  }, error = function(e) {
    cat("Fehler bei Raster-Prediction:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(prediction_raster)) {
    return(NULL)
  }
  
  cat("Prediction Raster erstellt", "\n")
  
  # resample prediction zu beiden ground truth rastern
  pred_for_suit <- resample(prediction_raster, ground_truth_suit, method = "bilinear")
  pred_for_pa <- resample(prediction_raster, ground_truth_pa, method = "bilinear")
  
  # === BERECHNUNG FÜR SUITABILITY-METRIKEN (RMSE, MAE, Pearson) ===
  gt_suit_values <- values(ground_truth_suit, na.rm = TRUE)
  pred_suit_values <- values(pred_for_suit, na.rm = TRUE)
  
  valid_suit_idx <- !is.na(gt_suit_values) & !is.na(pred_suit_values)
  gt_suit_clean <- gt_suit_values[valid_suit_idx] 
  pred_suit_clean <- pred_suit_values[valid_suit_idx]
  
  # === BERECHNUNG FÜR PA-METRIKEN (AUC) ===
  gt_pa_values <- values(ground_truth_pa, na.rm = TRUE)
  pred_pa_values <- values(pred_for_pa, na.rm = TRUE)
  
  valid_pa_idx <- !is.na(gt_pa_values) & !is.na(pred_pa_values)
  gt_pa_clean <- gt_pa_values[valid_pa_idx] 
  pred_pa_clean <- pred_pa_values[valid_pa_idx]
  
  if (length(gt_suit_clean) == 0 || length(gt_pa_clean) == 0) {
    cat("Keine gültigen Pixel für Evaluation", "\n")
    return(NULL)
  }
  
  # === METRIKEN BERECHNEN ===
  # Suitability-basierte Metriken
  rmse <- sqrt(mean((gt_suit_clean - pred_suit_clean)^2))
  mae <- mean(abs(gt_suit_clean - pred_suit_clean))
  pearson <- cor(gt_suit_clean, pred_suit_clean)
  
  # PA-basierte AUC
  auc <- tryCatch({
    round(as.numeric(auc(gt_pa_clean, pred_pa_clean)), 4)
  }, error = function(e) {
    cat("Warnung bei AUC-Berechnung:", e$message, "\n")
    return(NA)
  })
  
  # write results
  results <- data.frame(
    species = species_name,
    RMSE = rmse,
    MAE = mae,
    Pearson = pearson,
    AUC = auc,
    n_pixels_suit = length(gt_suit_clean),
    n_pixels_pa = length(gt_pa_clean)
  )
  
  cat("Species", species_name, "completed - AUC:", auc, "\n")
  
  return(list(
    metrics = results,
    model = maxent_model,
    prediction_raster = prediction_raster,
    ground_truth_pa = ground_truth_pa,
    ground_truth_suit = ground_truth_suit
  ))
}

