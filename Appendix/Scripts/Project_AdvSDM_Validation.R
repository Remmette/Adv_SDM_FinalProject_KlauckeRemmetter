#load packages
library(terra)
library(sf)
library(maxnet)
library(dplyr)
library(blockCV)
library(pROC)


#initialize functions-Script 
source("D:/Master/Adv_SDM/Project/Script/Project_AdvSDM_functions.R")

#set output directory
output_dir_maxent_cv <- "D:/Master/Adv_SDM/Project/Data/MaxEntCV/"


set.seed(21)

# load and name nlm-stacks 
nlms_remmetter <- rast("D:/Master/Adv_SDM/Data/Artificial_Landscape/Assignment03_NLMs_Roshar.tif")
names(nlms_remmetter) <- paste0("NLM", 1:12)
nlms_klaucke <- rast("D:/Master/Adv_SDM/Project/Data/Virtual Species/Klaucke/fantasyworld.tif")
names(nlms_klaucke) <- paste0("NLM", 1:12)
nlms_schneider <- rast("D:/Master/Adv_SDM/Project/Data/Virtual Species/Schneider/helmarcia.tif")
names(nlms_schneider) <- paste0("NLM", 1:12)


# load all species
ordner <- "D:/Master/Adv_SDM/Project/Data/Virtual Species/All_Species"

# list all RDS-files
files <- list.files(ordner, pattern = "\\.RDS$", full.names = TRUE)

species_list <- list()
envstacks_list <- list()
nlms_list <- list()

# loop over all files
for (f in files) {
  species <- readRDS(f)
  used_vars <- species[[1]]$details$variables
  
  # assign nlms depending on species origin
  if (grepl("remmetter", f)) {
    nlms_current <- nlms_remmetter
  } else if (grepl("schneider", f)) {
    nlms_current <- nlms_schneider
  } else if (grepl("klaucke", f)) {
    nlms_current <- nlms_klaucke
  } else {
    stop("Kein passender NLM-Stack für Datei: ", f)
  }
  
  env_stack <- get_env_data(species, nlms_current, used_vars)
  
  # extract species name from file
  species_name <- tools::file_path_sans_ext(basename(f))
  
  # save results in lists
  species_list[[species_name]] <- species
  envstacks_list[[species_name]] <- env_stack
  nlms_list[[species_name]] <- nlms_current
}



list_pa_data <- list() 

# create presence/absence datasets from envstacks (as sf)
for (i in seq_along(envstacks_list)) {
  pa_data <- sf::st_as_sf(envstacks_list[[i]], coords = c("x", "y"), crs = 32633)
  list_pa_data[[i]] <- pa_data
}

suitability_list <- list()

# extract suitability raster information for cv methods
for (i in seq_along(species_list)) {
  suitability <- unwrap(species_list[[i]][[1]]$suitab.raster)
  suitability_list[[i]] <- suitability
}


#Test
species_list[[i]][[1]]
suitability
plot(suitability_list[[5]])


###################################################### 
# extract species names and add to lists 
species_names <- names(species_list)

names(list_pa_data) <- species_names
names(suitability_list) <- species_names

# create predictor sets by extracting variables dynamicly (based on species)
predictor_sets <- lapply(species_list, function(species) {
  vars <- paste0("NLM", regmatches(species[[1]]$details$variables, regexpr("[0-9]+$", species[[1]]$details$variables)))
  if (is.null(vars) || length(vars) == 0) {
    warning(paste("Keine Variablen für Species", species, "gefunden"))
    return(character(0))
  }
  return(vars)
})
names(predictor_sets) <- species_names


# create species configurations by adding pa_data, predictor sets and response for each species into a list
species_config <- list()
for (species_name in species_names) {
  species_config[[species_name]] <- list(
    data = list_pa_data[[species_name]],
    predictors = predictor_sets[[species_name]],
    response = "response"
  )
}

# new lists to save the results in 
all_results <- list()
all_cv_metrics <- data.frame()
all_mean_fold_metrics <- data.frame()

# outer loop: run through all species and create needed variables
for (species_name in names(species_config)) {
  config <- species_config[[species_name]]
  suitability_raster <- suitability_list[[species_name]] 
  
  # get the cross validation methods based on the suitability raster for the species
  cv_methods <- get_cv_methods(suitability_raster)
  
  cat("Verarbeite Species:", species_name, "\n")
  
  # inner loop: run through all cross validation methods
  for (cv_method_name in names(cv_methods)) {
    
    cv_folds <- cv_methods[[cv_method_name]](config$data)
    combination_name <- paste0(species_name, "_", cv_method_name)
    
    cat("  - CV-Methode:", cv_method_name, "\n")
    
    result <- tryCatch({
      run_species_cv(
        species_data = config$data,
        predictors_cols = config$predictors,
        response_col = config$response,
        cv_folds = cv_folds,
        species_name = combination_name
      )
    }, error = function(e) {
      warning(paste("Fehler in", combination_name, ":", e$message))
      return(NULL)
    })
    
    # Nur hinzufügen wenn erfolgreich
    if (!is.null(result)) {
      all_results[[combination_name]] <- result
      all_cv_metrics <- rbind(all_cv_metrics, result$cv_metrics)
      all_mean_fold_metrics <- rbind(all_mean_fold_metrics, result$mean_fold_metrics)
    }
  }
}

# checks, tests and result overview
cat("\nVerarbeitung abgeschlossen für", length(species_names), "Species\n")
cat("Anzahl erfolgreiche Kombinationen:", nrow(all_cv_metrics), "\n\n")

print("Cross-Validation Metriken für alle Species:")
print(all_cv_metrics)
print("\nDurchschnittliche Fold-Metriken für alle Species:")
print(all_mean_fold_metrics)

# save as csv-files
write.csv(all_cv_metrics, "D:/Master/Adv_SDM/Project/Data/cv_metrics_all_species.csv", row.names = FALSE)
write.csv(all_mean_fold_metrics, "D:/Master/Adv_SDM/Project/Data/mean_fold_metrics_all_species.csv", row.names = FALSE)


