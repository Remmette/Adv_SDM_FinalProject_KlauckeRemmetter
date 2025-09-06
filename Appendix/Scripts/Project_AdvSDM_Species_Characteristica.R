
#load packages
library(terra)
library(sf)
library(dplyr)

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


species_test <- species_list[["remmetter_species2"]]


df <- species_test[[3]]$sample.points
species_test[[3]]
nrow(species_test[[3]]$sample.points)


# Prävalenz berechnen (Anteil Presence == 1)
prevalence <- mean(df$Observed == 1, na.rm = TRUE)
prevalence

##############################

library(terra)
library(spdep)
library(sp)
library(raster)

analyze_species <- function(loaded_data) {
  results <- list()
  
  for (name in names(loaded_data)) {
    obj <- loaded_data[[name]]
    
    # Sample points extrahieren
    if (length(obj) >= 3 && is.data.frame(obj[[3]]$sample.points)) {
      df <- obj[[3]]$sample.points
      
      # Prävalenz
      prevalence <- mean(df$Observed == 1, na.rm = TRUE)
      
      # Moran's I
      coords <- df[, c("x", "y")]
      sp_obj <- SpatialPoints(coords)
      nb <- dnearneigh(sp_obj, 0, 50)  # Distanz anpassen
      lw <- nb2listw(nb, zero.policy = TRUE)
      moran_result <- moran.test(df$Observed, lw, zero.policy = TRUE)
      moran_i <- moran_result$estimate[1]
      
      # Fragmentierung
      rast <- terra::rast(
        terra::ext(range(df$x), range(df$y)),
        resolution = 1
      )
      values(rast) <- NA
      
      spdf <- df[, c("x", "y", "Observed")]
      rast_pts <- terra::vect(spdf, geom = c("x", "y"))
      pres_rast <- terra::rasterize(rast_pts, rast, field = "Observed", fun = "last")
      
      pres_r <- raster::raster(pres_rast)
      pres_r[is.na(pres_r)] <- 0
      
      clumped <- clump(pres_r, directions = 8)
      n_fragments <- length(unique(na.omit(values(clumped))))
      
      results[[name]] <- data.frame(
        species_name = name,
        prevalence = prevalence,
        moran_i = moran_i,
        fragments = n_fragments
      )
    }
  }
  
  return(do.call(rbind, results))
}

results <- analyze_species(species_list)
results

# Visualisierung
plot_results <- function(results) {
  library(ggplot2)
  library(gridExtra)
  
  # Prevalence Plot
  p1 <- ggplot(results, aes(x = reorder(species_name, prevalence), y = prevalence)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Prävalenz", x = "Species", y = "Prevalence") +
    theme_minimal()
  
  # Moran's I Plot
  p2 <- ggplot(results, aes(x = reorder(species_name, moran_i), y = moran_i)) +
    geom_col(fill = "darkgreen", alpha = 0.7) +
    coord_flip() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = "Moran's I", x = "Species", y = "Moran's I") +
    theme_minimal()
  
  # Fragments Plot
  p3 <- ggplot(results, aes(x = reorder(species_name, fragments), y = fragments)) +
    geom_col(fill = "orange", alpha = 0.7) +
    coord_flip() +
    labs(title = "Fragmente", x = "Species", y = "Anzahl Fragmente") +
    theme_minimal()
  
  # Kombinierter Plot
  grid.arrange(p1, p2, p3, ncol = 1)
}

# Scatterplot Matrix
plot_correlations <- function(results) {
  library(GGally)
  
  ggpairs(results[, c("prevalence", "moran_i", "fragments")],
          title = "Korrelationen zwischen Metriken") +
    theme_minimal()
}

# Simple Korrelationsmatrix
plot_correlations <- function(results) {
  # Korrelationsmatrix berechnen
  cor_matrix <- cor(results[, c("prevalence", "moran_i", "fragments")])
  
  # Heatmap
  corrplot::corrplot(cor_matrix, 
                     method = "color",
                     type = "upper",
                     addCoef.col = "black",
                     tl.col = "black")
}

show_correlations <- function(results) {
  cor(results[, c("prevalence", "moran_i", "fragments")])
}

plot_results(results)
plot_correlations(results)
show_correlations(results)

# write results as csv 
write.csv(results, "D:/Master/Adv_SDM/Project/Data/Species_characteristica.csv", row.names = FALSE)

