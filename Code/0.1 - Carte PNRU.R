library(sf)
library(dplyr)
library(tmap)
library(haven)
library(readr)
library(readxl)

# 1. Tes chemins exacts
out_path <- "C:/Users/camil/Documents/PNRU-et-Resultats-scolaires/Output groupes delim"
shp_path <- "C:/Users/camil/Dropbox/2026_Camille_Favre/Analyse/input_ANRU/data_Quartiers/data_Shapefiles_quartiers"

# 2. Charger les Stable Units déjà générées (1 451 unités en Lambert-93 / EPSG:2154)
su_new15 <- readRDS(file.path(out_path, "versions_new15.rds")) %>%
  st_make_valid() %>%
  st_cast("MULTIPOLYGON") # homogénéise POLYGON et MULTIPOLYGON pour tmap

# Vérification : on retrouve bien les 735 unités PNRU et 716 non-PNRU
table(su_new15$has_pnru)

# 3. Charger les ZUS brutes pour comparer le contour historique PNRU au contour harmonisé
zus_pnru <- st_read(file.path(shp_path, "zus_art6/zus-art6.shp"), quiet = TRUE) %>%
  filter(Categ_quar %in% c("PNRU conventionnés", "PNRU non-conventionnés")) %>%
  st_make_valid() %>%
  st_transform(2154) %>%
  st_cast("MULTIPOLYGON") %>%
  select(CODQUART, Categ_quar)

# Séparer les unités traitées PNRU (735) et les unités non-PNRU (716)
su_pnru     <- su_new15 %>% filter(has_pnru == TRUE)
su_non_pnru <- su_new15 %>% filter(has_pnru == FALSE)

# 4. Générer la carte interactive HTML des quartiers PNRU
tmap_mode("view")

carte_pnru <-
  # Couche 1 : Les 735 Stable Units PNRU
  tm_shape(st_transform(su_pnru, 4326), name = "Stable Units PNRU (735)") +
  tm_polygons(
    id = "stable_unit_id",
    col = "#e41a1c",
    alpha = 0.45,
    border.col = "#99000d",
    lwd = 1.2,
    popup.vars = c(
      "ID Unité Stable" = "stable_unit_id",
      "Codes quartiers" = "members",
      "Sources"         = "sources",
      "Surface (m²)"    = "area"
    )
  ) +
  # Couche 2 : Périmètre ZUS PNRU brut (sans la case Missing)
  tm_shape(st_transform(zus_pnru, 4326), name = "Contours ZUS PNRU bruts") +
  tm_polygons(
    id = "CODQUART",
    col = "Categ_quar",
    title = "Conventionnement ZUS",
    palette = c("PNRU conventionnés" = "#ff7f00", "PNRU non-conventionnés" = "#984ea3"),
    showNA = FALSE,
    alpha = 0.30,
    border.col = "#ff7f00",
    lwd = 1.8,
    popup.vars = c("Code ZUS" = "CODQUART", "Type PNRU" = "Categ_quar")
  ) +
  # Couche 3 : Unités non-PNRU (716) en bleu distinct
  tm_shape(st_transform(su_non_pnru, 4326), name = "Quartiers non-PNRU (716)") +
  tm_polygons(
    id = "stable_unit_id",
    col = "#377eb8",
    alpha = 0.40,
    border.col = "#08519c",
    lwd = 1.2,
    popup.vars = c(
      "ID Unité Stable" = "stable_unit_id",
      "Codes quartiers" = "members",
      "Sources"         = "sources"
    )
  ) +
  tm_layout(
    title = "Délimitation des 735 quartiers PNRU (Stable Units new15)",
    legend.outside = TRUE
  )

tmap_save(carte_pnru, file.path(out_path, "carte_quartiers_pnru_new15.html"))