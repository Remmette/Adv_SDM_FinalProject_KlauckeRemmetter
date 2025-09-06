# ============================================================================
# SDM Cross-Validation vs Map Accuracy Analysis
# ============================================================================

# Pakete laden
library(dplyr)
library(ggplot2)
library(tidyr)
library(corrplot)

# Daten einlesen

map_accuracy <- read.csv("D:/Master/Adv_SDM/Project/Data/map_accuracy_metrics.csv")
cv_metrics <- read.csv("D:/Master/Adv_SDM/Project/Data/mean_fold_metrics_all_species.csv")
species_types <- read.csv("D:/Master/Adv_SDM/Project/Data/Species_characteristica.csv")

## !! Metriken je nach gewünschter Auswertung mit "Find and Replace" (Strg+F) ersetzen (ist leider ein bisschen viel :) ) -> schlauer wäre wohl die Metrik im Code als Variable aufzubauen.... ##

# === 1. DATENAUFBEREITUNG ===
# Spalten umbenennen
cv_metrics <- rename(cv_metrics, RMSE=RMSE_mean, MAE=MAE_mean, Pearson=Pearson_mean, AUC=AUC_mean)
cv_metrics

# Map Accuracy aufbereiten
map_clean <- map_accuracy %>%
  separate(species, into = c("author", "species_num"), sep = "_") %>%
  mutate(species_id = paste(author, species_num, sep = "_"))

# CV Metrics aufbereitung  
cv_clean <- cv_metrics %>%
  separate(species, into = c("author", "species_num", "cv_method"), sep = "_") %>%
  mutate(species_id = paste(author, species_num, sep = "_"))

# === 2. BEST METHOD FREQUENCY ===
# Für jede Art beste CV-Methode finden (basierend auf Metrik)
best_methods <- cv_clean %>%
  group_by(species_id) %>%
  slice_max(AUC, n = 1) %>%
  ungroup() %>%
  count(cv_method) %>%
  mutate(percentage = n/sum(n)*100)

print("Best Method Frequency:")
print(best_methods)

# === 3. BIAS BERECHNUNG ===
# Map Accuracy mit CV-Metriken joinen
combined <- cv_clean %>%
  left_join(map_clean %>% dplyr::select(species_id, map_auc = AUC), 
            by = "species_id") %>%
  mutate(bias = AUC - map_auc)

# Durchschnittlicher Bias pro Methode
avg_bias <- combined %>%
  group_by(cv_method) %>%
  summarise(
    mean_bias = mean(bias, na.rm = TRUE),
    sd_bias = sd(bias, na.rm = TRUE),
    .groups = "drop"
  )

print("Average Bias (CV_AUC - Map_AUC):")
print(avg_bias)

# === 4. RANKING CORRELATION ===
# Ranking-Korrelation zwischen CV und Map Accuracy
ranking_data <- combined %>%
  dplyr::select(species_id, cv_method, cv_auc = AUC, map_auc) %>%
  pivot_wider(names_from = cv_method, values_from = cv_auc) %>%
  dplyr::select(-species_id)

ranking_cors <- cor(ranking_data %>% dplyr::select(-map_auc), 
                    ranking_data$map_auc, 
                    method = "spearman", 
                    use = "complete.obs")

print("Ranking Correlations (Spearman):")
print(ranking_cors)

# === 5. VISUALISIERUNGEN ===

# Plot 1: Best Method Frequency
p1 <- ggplot(best_methods, aes(x = reorder(cv_method, n), y = n)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = paste0(n, " (", round(percentage, 1), "%)")), 
            hjust = -0.1) +
  coord_flip() +
  labs(title = "Best Method Frequency", 
       x = "CV Method", y = "Count") +
  theme_minimal()

# Plot 2: Bias Distribution
p2 <- ggplot(combined, aes(x = cv_method, y = bias)) +
  geom_boxplot(fill = "lightcoral", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Bias Distribution (CV - Map Accuracy)", 
       x = "CV Method", y = "Bias (AUC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 3: Scatter Plot - CV vs Map Accuracy
p3 <- ggplot(combined, aes(x = map_auc, y = AUC, color = cv_method)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  facet_wrap(~cv_method) +
  labs(title = "CV vs Map Accuracy", 
       x = "Map Accuracy (AUC)", y = "CV AUC") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 4: Ranking Correlation Heatmap
ranking_cors_df <- data.frame(
  method = names(ranking_cors[,1]),
  correlation = ranking_cors[,1]
)

p4 <- ggplot(ranking_cors_df, aes(x = method, y = 1, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = round(correlation, 3)), color = "black", size = 4) +
  scale_fill_gradient2(low = "red", mid = "yellow", high = "green", 
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Ranking Correlation (Spearman) with Map Accuracy", 
       x = "CV Method", y = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

# Alle Plots anzeigen
print(p1)
print(p2)
print(p3)
print(p4)

# === 6. ZUSAMMENFASSUNG ===
summary_stats <- combined %>%
  group_by(cv_method) %>%
  summarise(
    n_species = n(),
    mean_bias = round(mean(bias, na.rm = TRUE), 4),
    rmse = round(sqrt(mean((AUC - map_auc)^2, na.rm = TRUE)), 4),
    correlation = round(cor(AUC, map_auc, use = "complete.obs"), 3),
    .groups = "drop"
  ) %>%
  arrange(abs(mean_bias))

print("SUMMARY STATISTICS:")
print(summary_stats)


##################################################################################
# Spezies Eigenschaften 

# Species characteristics mit den CV-Ergebnissen verbinden
combined_traits <- combined %>%
  left_join(species_types, by = c("species_id" = "species_name"))

# --- Scatterplots: Bias vs Species Eigenschaften ---
p5 <- ggplot(combined_traits, aes(x = prevalence, y = bias, color = cv_method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  labs(title = "Bias vs Prävalenz", x = "Prävalenz", y = "Bias (CV - Map AUC)") +
  theme_minimal()

p6 <- ggplot(combined_traits, aes(x = moran_i, y = bias, color = cv_method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  labs(title = "Bias vs Moran's I", x = "Moran's I", y = "Bias (CV - Map AUC)") +
  theme_minimal()

p7 <- ggplot(combined_traits, aes(x = fragments, y = bias, color = cv_method)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() + # oft log besser für Fragmente
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  labs(title = "Bias vs Fragmentierung", x = "Anzahl Fragmente (log-Skala)", y = "Bias (CV - Map AUC)") +
  theme_minimal()

print(p5)
print(p6)
print(p7)

# --- Korrelationen numerisch berechnen ---
correlations <- combined_traits %>%
  group_by(cv_method) %>%
  summarise(
    cor_prev = cor(prevalence, bias, use = "complete.obs"),
    cor_moran = cor(moran_i, bias, use = "complete.obs"),
    cor_frag = cor(fragments, bias, use = "complete.obs"),
    .groups = "drop"
  )

print("Korrelationen zwischen Bias und Spezies-Eigenschaften:")
print(correlations)

# --- Heatmap der Korrelationen ---
cor_long <- correlations %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "trait", values_to = "correlation")

p8 <- ggplot(cor_long, aes(x = trait, y = cv_method, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = round(correlation, 2)), color = "black") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "Korrelationen: Bias vs Species-Eigenschaften",
       x = "Eigenschaft", y = "CV-Methode") +
  theme_minimal()

print(p8)


# --- CV-Ergebnisse direkt mit Species-Eigenschaften verbinden ---
cv_traits <- cv_clean %>%
  left_join(species_types, by = c("species_id" = "species_name"))

# --- Scatterplots: AUC vs Species Eigenschaften ---
p9 <- ggplot(cv_traits, aes(x = prevalence, y = AUC, color = cv_method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  labs(title = "AUC vs Prävalenz", x = "Prävalenz", y = "CV-AUC") +
  theme_minimal()

p10 <- ggplot(cv_traits, aes(x = moran_i, y = AUC, color = cv_method)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  labs(title = "AUC vs Moran's I", x = "Moran's I", y = "CV-AUC") +
  theme_minimal()

p11 <- ggplot(cv_traits, aes(x = fragments, y = AUC, color = cv_method)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "black") +
  labs(title = "AUC vs Fragmentierung", 
       x = "Anzahl Fragmente (log-Skala)", 
       y = "CV-AUC") +
  theme_minimal()

print(p9)
print(p10)
print(p11)

# --- Korrelationen numerisch berechnen (CV-Ergebnisse) ---
correlations_raw <- cv_traits %>%
  group_by(cv_method) %>%
  summarise(
    cor_prev = cor(prevalence, AUC, use = "complete.obs"),
    cor_moran = cor(moran_i, AUC, use = "complete.obs"),
    cor_frag = cor(fragments, AUC, use = "complete.obs"),
    .groups = "drop"
  )

print("Korrelationen zwischen AUC und Species-Eigenschaften:")
print(correlations_raw)


# --- Heatmap der Korrelationen ---
cor_long_raw <- correlations_raw %>%
  pivot_longer(cols = starts_with("cor_"), 
               names_to = "trait", 
               values_to = "correlation")

p12 <- ggplot(cor_long_raw, aes(x = trait, y = cv_method, fill = correlation)) +
  geom_tile() +
  geom_text(aes(label = round(correlation, 2)), color = "black") +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue", midpoint = 0) +
  labs(title = "Korrelationen: CV-AUC vs Species-Eigenschaften",
       x = "Eigenschaft", y = "CV-Methode") +
  theme_minimal()

print(p12)
