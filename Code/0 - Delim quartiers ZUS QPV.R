# ==============================================================================
# 0. Setup (code exact de 03_d_Stable_units_updated.qmd)
# ==============================================================================
library(sf)
library(dplyr)
library(purrr)
library(readr)
library(igraph)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rmapshaper)
library(tmap)
library(haven)
library(tidyr)
library(fuzzyjoin)
library(stringr)
library(tidyverse)
# library(RStata)   # unused in this document

# SEULE MODIFICATION : tes chemins locaux
data_path   <- "C:/Users/camil/Dropbox/2026_Camille_Favre/Analyse/input_ANRU/data_Quartiers/data_Shapefiles_quartiers/"
output_path <- "C:/Users/camil/Documents/PNRU-et-Resultats-scolaires/Output groupes delim" # dossier contenant quartiers_npnru.dta


# ==============================================================================
# A. Stable units (code exact lettre pour lettre)
# ==============================================================================

### Helper functions

THR    <- 0.15      # overlap threshold for merging polygons
METRIC <- "small"  

collapse_unique_tokens <- function(x, sep = "|") {
  toks <- unlist(strsplit(x[!is.na(x) & x != ""], "\\|"), use.names = FALSE)
  toks <- trimws(toks); toks <- unique(toks[toks != "" & toks != "NA"])
  paste(toks, collapse = sep)
}

# Standardise a layer 
std_layer <- function(x, code_col, src, keep = character()) {
  x <- x %>% st_make_valid() %>% st_transform(2154)
  out <- x %>%
    mutate(key = paste0(src, "_", row_number()),
           CODQUART = as.character(.data[[code_col]]),
           source = src) %>%
    select(key, CODQUART, source, any_of(keep))
  out$area <- as.numeric(st_area(out))
  out
}

# All pairwise overlaps between A and B (both need `key`, `area`). POLYGON part only.
cross_pairs <- function(A, B) {
  A2 <- A %>% transmute(key_a = key, area_a = area)
  B2 <- B %>% transmute(key_b = key, area_b = area)
  ix <- suppressWarnings(st_intersection(A2, B2))
  empty <- tibble(key_a = character(), key_b = character(), area_a = numeric(),
                  area_b = numeric(), overlap_m2 = numeric(), share_a = numeric(),
                  share_b = numeric(), share_small = numeric(), jaccard = numeric())
  if (nrow(ix) == 0) return(empty)
  ix <- suppressWarnings(st_collection_extract(ix, "POLYGON"))
  ix$overlap_m2 <- as.numeric(st_area(ix))
  ix %>% st_drop_geometry() %>%
    group_by(key_a, key_b, area_a, area_b) %>%
    summarise(overlap_m2 = sum(overlap_m2), .groups = "drop") %>%
    filter(overlap_m2 > 0) %>%
    mutate(share_a     = overlap_m2 / area_a,
           share_b     = overlap_m2 / area_b,
           share_small = overlap_m2 / pmin(area_a, area_b),
           jaccard     = overlap_m2 / (area_a + area_b - overlap_m2))
}

passes <- function(p, thr = THR, metric = METRIC) {
  switch(metric,
         small   = p$share_small >= thr,
         jaccard = p$jaccard >= thr,
         both    = pmin(p$share_a, p$share_b) >= thr)
}

# Edges within one layer set. min_overlap_m2 = NULL -> touching counts  
poly_edges <- function(P, min_overlap_m2 = NULL) {
  if (is.null(min_overlap_m2)) {
    nb <- st_intersects(P)
    e  <- tibble(i = rep(seq_along(nb), lengths(nb)), j = unlist(nb)) %>% filter(i < j)
    return(tibble(key_a = P$key[e$i], key_b = P$key[e$j]))
  }
  
  P2 <- P %>% select(key)
  ix <- suppressWarnings(st_intersection(P2, P2)) %>%     # columns: key, key.1
    filter(key < key.1)                                    # each unordered pair once, no self pairs
  if (nrow(ix) == 0) return(tibble(key_a = character(), key_b = character()))
  
  ix <- suppressWarnings(st_collection_extract(ix, "POLYGON"))   # polygon part only
  ix$ov <- as.numeric(st_area(ix))
  
  ix %>% st_drop_geometry() %>%
    group_by(key_a = key, key_b = key.1) %>%
    summarise(ov = sum(ov), .groups = "drop") %>%
    filter(ov > min_overlap_m2) %>%
    select(key_a, key_b)
}

components_from_edges <- function(keys, edges) {
  stopifnot(!anyDuplicated(keys))
  g <- graph_from_data_frame(edges[, c("key_a", "key_b")], directed = FALSE,
                             vertices = data.frame(name = keys))
  m <- components(g)$membership
  tibble(key = names(m), comp = as.integer(m))
}

as_nodes <- function(x) x %>% transmute(key, members = CODQUART, sources = source, n_poly = 1L, area)

dissolve_nodes <- function(P, comp_tbl, prefix) {
  P %>% left_join(comp_tbl, by = "key") %>%
    group_by(comp) %>%
    summarise(members  = collapse_unique_tokens(members),
              sources  = collapse_unique_tokens(sources),
              n_poly   = sum(n_poly),
              geometry = st_union(geometry), .groups = "drop") %>%
    st_as_sf() %>%
    mutate(key = sprintf("%s%05d", prefix, comp)) %>%
    { .$area <- as.numeric(st_area(.)); . } %>%
    select(key, members, sources, n_poly, area)
}

# stable ID: smallest member code (not the igraph component number)
finalize_ids <- function(u) {
  u %>%
    mutate(first_code = vapply(strsplit(members, "\\|"), function(x) sort(x)[1], ""),
           n_codes    = lengths(strsplit(members, "\\|"))) %>%
    arrange(first_code) %>%
    mutate(stable_unit_id = make.unique(paste0("SU_", first_code))) %>%
    select(stable_unit_id, members, sources, n_codes, n_poly, area)
}


### 1. Load boundaries and treatment tables

zus_raw <- st_read(file.path(data_path, "zus_art6/zus-art6.shp"), quiet = TRUE)
q15_raw <- st_read(file.path(data_path, "qpv-2015/qp-politiquedelaville-shp/QP_METROPOLE_LB93.shp"), quiet = TRUE)
q24_raw <- st_read(file.path(data_path, "qpv-2024/QPV_2024/SHP/QP2024_France_hexagonale_LB93.shp"), quiet = TRUE)

zus <- std_layer(zus_raw, "CODQUART", "ZUS")
q15 <- std_layer(q15_raw, "CODE_QP",  "QPV2015")
q24 <- std_layer(q24_raw, "code_qp",  "QPV2024", keep = "insee_com")


qrr_raw <- st_read(file.path(data_path, "qrr/contours_QRR.shp"), quiet = TRUE)
qrr <- std_layer(qrr_raw, "code_qrr", "QRR")


cucs <- std_layer(st_read(file.path(data_path, "Cucs/HZUS_Scan25_LB93.shp"), quiet = TRUE),
                  "CODQUART", "CUCS")

polys <- bind_rows(zus, q15, q24)
stopifnot(!anyDuplicated(polys$key))

# duplicate codes ? should be empty
polys %>% st_drop_geometry() %>% count(source, CODQUART) %>% filter(n > 1)


### 2. Treatment lookup

pnru_codes <- zus_raw %>% st_drop_geometry() %>%
  filter(Categ_quar %in% c("PNRU conventionnés", "PNRU non-conventionnés")) %>%
  transmute(CODQUART = as.character(CODQUART), pnru_subtype = Categ_quar) %>% distinct()

### 2. Treatment lookup (sans quartiers_npnru.dta pour l'instant)
pnru_codes <- zus_raw %>% st_drop_geometry() %>%
  filter(Categ_quar %in% c("PNRU conventionnés", "PNRU non-conventionnés")) %>%
  transmute(CODQUART = as.character(CODQUART), pnru_subtype = Categ_quar) %>% distinct()

npnru_codes <- tibble(CODQUART = character())

treat_lookup <- full_join(
  pnru_codes  %>% mutate(has_pnru  = TRUE),
  npnru_codes %>% mutate(has_npnru = TRUE), by = "CODQUART") %>%
  mutate(has_pnru  = coalesce(has_pnru,  FALSE),
         has_npnru = coalesce(has_npnru, FALSE))
stopifnot(!anyDuplicated(treat_lookup$CODQUART))

# Treated codes with no polygon (overseas, QP999 ...) -> not part of the units
missing_tbl <- tibble(CODQUART = setdiff(treat_lookup$CODQUART, polys$CODQUART)) %>%
  mutate(dep = substr(CODQUART, 3, 4))  
table(missing_tbl$dep)
if (any(!missing_tbl$dep %in% c("97", "99")))
  warning("Unexpected missing treated codes: ", paste(head(missing_tbl$CODQUART[!missing_tbl$dep %in% c("97","99")]), collapse = ", "))


### 3. Building stable units

classify_units <- function(units, lookup) {
  flags <- units %>% st_drop_geometry() %>%
    select(stable_unit_id, members) %>%
    separate_rows(members, sep = "\\|") %>%
    rename(CODQUART = members) %>% distinct() %>%
    left_join(lookup %>% select(CODQUART, has_pnru, has_npnru), by = "CODQUART") %>%
    group_by(stable_unit_id) %>%
    summarise(has_pnru = any(has_pnru %in% TRUE),
              has_npnru = any(has_npnru %in% TRUE),
              n_treated_codes = sum(!is.na(has_pnru)), .groups = "drop")
  units %>% left_join(flags, by = "stable_unit_id") %>%
    mutate(treatment_category = factor(case_when(
      has_pnru & has_npnru  ~ "PNRU + NPNRU",
      has_pnru              ~ "PNRU only",
      has_npnru             ~ "NPNRU only",
      TRUE                  ~ "Control"),
      levels = c("Control", "PNRU only", "NPNRU only", "PNRU + NPNRU")))
}

# Baseline: all three layers, linked if they touch (min_overlap_m2 = NULL)
# or only if they share positive area (min_overlap_m2 = 1)
build_baseline <- function(polys, min_overlap_m2 = NULL) {
  P <- as_nodes(polys)
  comp <- components_from_edges(P$key, poly_edges(P, min_overlap_m2))
  dissolve_nodes(P, comp, "SU_") %>% finalize_ids()
}

# One merging phase between node sets A and B
run_phase <- function(A, B, thr, metric, mode, prefix) {
  pairs <- cross_pairs(A, B)
  q     <- pairs[passes(pairs, thr, metric), ]
  P     <- bind_rows(A, B)
  edges <- if (mode == "new") {
    q[, c("key_a", "key_b")]                       # edges ONLY from qualifying pairs
    # does not merge qualifying polygons that do not overlap more than 15%
  } else {                                         # old way: re-run st_intersects (touching merges)
    # qualifying pairs poligons overlapping by < 15% are still merged into a single unit
    keep <- unique(c(q$key_a, q$key_b))
    poly_edges(P[P$key %in% keep, ])
  }
  dissolve_nodes(P, components_from_edges(P$key, edges), prefix)
}

# Phase 1: ZUS x QPV2015, Phase 2: (phase-1 units) x QPV2024
build_15 <- function(zus, q15, q24, thr = THR, metric = METRIC, mode = c("new", "old")) {
  mode <- match.arg(mode)
  p1 <- run_phase(as_nodes(zus), as_nodes(q15), thr, metric, mode, "P1_")
  run_phase(p1, as_nodes(q24), thr, metric, mode, "SU_") %>% finalize_ids()
}


### 4. Creating datasets

versions <- list(
  baseline_touch   = build_baseline(polys, NULL),   # original: touching merges
  baseline_overlap = build_baseline(polys, 1),      # positive-area overlap only
  old15            = build_15(zus, q15, q24, mode = "old"),
  new15            = build_15(zus, q15, q24, mode = "new")
) %>% map(classify_units, lookup = treat_lookup)


### 5. Checks

version_summary <- imap_dfr(versions, function(u, nm) {
  d <- st_drop_geometry(u)
  tibble(version = nm, n_units = nrow(d), n_polygons = sum(d$n_poly),
         n_multi_polygon_units = sum(d$n_poly > 1), max_polygons_in_unit = max(d$n_poly),
         n_treated_units = sum(d$treatment_category != "Control"),
         n_control_units = sum(d$treatment_category == "Control"))
})
print(version_summary)
# new version should have 43 stable units more
cat("Input polygons:", nrow(polys), "\n")
stopifnot(all(version_summary$n_polygons == nrow(polys)))


# ==============================================================================
# Export sur disque pour ne plus avoir à le recalculer
# ==============================================================================
saveRDS(versions$new15, file.path(output_path, "versions_new15.rds"))
st_write(versions$new15, file.path(output_path, "versions_new15.gpkg"), delete_dsn = TRUE)
