
#load packages
library(terra)
library(sf)
library(maxnet)
library(dplyr)
library(pROC)

#initialize functions-Script 
source("D:/Master/Adv_SDM/Project/Script/Project_AdvSDM_functions.R")

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
nlms_list <- list()

# loop over all files
for (f in files) {
  species <- readRDS(f)
  
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
  
  # extract species name from file
  species_name <- tools::file_path_sans_ext(basename(f))
  
  # save results in lists
  species_list[[species_name]] <- species
  nlms_list[[species_name]] <- nlms_current
}

species_names <- names(species_list)

# create predictor sets based on species variables 
predictor_sets <- lapply(species_list, function(species) {
  original_vars <- species[[1]]$details$variables
  # map to nlm format
  numbers <- regmatches(original_vars, regexpr("[0-9]+$", original_vars))
  mapped_vars <- paste0("NLM", numbers)
  
  if (is.null(mapped_vars) || length(mapped_vars) == 0) {
    warning("Keine Variablen für Species gefunden")
    return(character(0))
  }
  return(mapped_vars)
})

names(predictor_sets) <- species_names


map_accuracy_results <- list()
all_metrics <- data.frame()
failed_species <- character()

cat("=== MAP ACCURACY EVALUATION ===")
cat("Gefundene Species:", length(species_names))
cat("Species Namen:", paste(species_names, collapse=", "))

# loop throug species and evaluate map accuracy
for (species_name in species_names) {
  
  result <- tryCatch({
    evaluate_map_accuracy(
      species_obj = species_list[[species_name]],
      nlm_stack = nlms_list[[species_name]], 
      used_vars = predictor_sets[[species_name]],
      species_name = species_name
    )
  }, error = function(e) {
    cat("Fehler bei", species_name, ":", e$message)
    return(NULL)
  })
  
  if (!is.null(result)) {
    map_accuracy_results[[species_name]] <- result
    all_metrics <- rbind(all_metrics, result$metrics)
  } else {
    failed_species <- c(failed_species, species_name)
  }
}



# overview over processing
cat("=== ZUSAMMENFASSUNG MAP ACCURACY ===")
cat("Anzahl berechneter Species:", nrow(all_metrics))
cat("Anzahl fehlende Species:", length(failed_species))
if (length(failed_species) > 0) {
  cat("Fehlgeschlagene Species:", paste(failed_species, collapse=", "))
}


# show results and write to csv-file
if (nrow(all_metrics) > 0) {
  print("Map Accuracy Metriken für alle Species:")
  print(all_metrics)
  
  write.csv(all_metrics, "D:/Master/Adv_SDM/Project/Data/map_accuracy_metrics.csv", row.names = FALSE)
  
  # save prediction rasters as csv
  output_dir <- "D:/Master/Adv_SDM/Project/Data/MaxEnt_Predictions/"
  dir.create(output_dir, showWarnings = FALSE)
  
  for (species_name in names(map_accuracy_results)) {
    if (!is.null(map_accuracy_results[[species_name]]$prediction_raster)) {
      writeRaster(
        map_accuracy_results[[species_name]]$prediction_raster,
        file.path(output_dir, paste0("prediction_", species_name, ".tif")),
        overwrite = TRUE
      )
      cat("Prediction Raster gespeichert für:", species_name)
    }
  }
} 


#Tests and visualizations
plot(map_accuracy_results[["remmetter_species2"]]$prediction_raster)

map_accuracy_results[["remmetter_species2"]]

species2 <- readRDS("D:/Master/Adv_SDM/Data/Virtual_Species/species2.RDS")

plot(rast(species2[[1]]$suitab.raster), main = "Suitability")
plot(rast(species2[[2]]$pa.raster), main = "Presence/Absence")


test_species <- species_list[["remmetter_species2"]]
raster <- unwrap(test_species[[2]]$pa.raster)
plot(raster)
raster
