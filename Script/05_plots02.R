library(cowplot)
library(ggplot2)
library(dplyr)
library(grid)
library(sf)

# Load datasets ----
infection_county   <- readRDS("ProcessedData/map_county_infection_proportion.rds")
infection_district <- readRDS("ProcessedData/map_district_infection_proportion.rds")

tx_counties <- counties(
  state = "TX",
  year  = 2024
) %>%
  mutate(county = toupper(NAME)) %>%
  select(county)


# Functions -----

compute_indirect_risk_county_new <- function(trans_mat, county_name = "GAINES", threshold = 0.5) {
  county_upper <- toupper(county_name)
  county_names <- rownames(trans_mat)
  
  if (!(county_upper %in% county_names))
    stop("County not found in matrix")
  
  # Direct transmission vector P_j→i
  Pji_vec <- trans_mat[county_upper, ]
  
  # Identify candidate first-generation counties k
  Pjk_vec  <- trans_mat[county_upper, ]          # Same row: j → k
  k_indices <- which(Pjk_vec >= threshold & !is.na(Pjk_vec))
  K         <- names(k_indices)
  
  # Output vector
  SG_vec <- rep(NA_real_, length(Pji_vec))
  names(SG_vec) <- names(Pji_vec)
  
  for (i in names(Pji_vec)) {
    
    Pji <- Pji_vec[i]  # direct probability j → i
    
    # If no k, then SG probability = Pji
    if (length(K) == 0) {
      SG_vec[i] <- Pji
      next
    }
    
    # Part 1: (1 - Pji)
    term1 <- (1 - Pji)
    
    # Part 2: product over k of (1 - Pjk * Pki)
    prod_term <- 1
    for (k in K) {
      Pjk <- Pjk_vec[k]  # j → k
      Pki <- trans_mat[k, i]  # k → i
      if (is.na(Pki)) next
      prod_term <- prod_term * (1 - Pjk * Pki)
    }
    
    # Final formula:
    # P_SG_ij = 1 - term1 * prod_term
    SG_vec[i] <- 1 - term1 * prod_term
  }
  
  # Remove self-transmission
  SG_vec[county_upper] <- NA
  
  return(SG_vec)
}


compute_indirect_risk_district_new <- function(district_mat,
                                               county_mat,
                                               county_name = "GAINES",
                                               threshold = 0.5) {
  # j0 = origin county
  j0 <- toupper(county_name)
  
  district_names <- rownames(district_mat)   # i: districts
  county_cols    <- colnames(district_mat)   # j,k: counties (as seen by districts)
  county_rows    <- rownames(county_mat)     # j: counties (origins in county_mat)
  county_cols_cc <- colnames(county_mat)     # k: counties (destinations in county_mat)
  
  
  # Direct probability from county j0 to each district i: P(j0 -> i)
  Pji_vec <- district_mat[, j0]   # names = districts
  # First-generation counties k from j0: P(j0 -> k)
  Pjk_vec <- county_mat[j0, ]     # row = j0, cols = counties k
  
  # Candidate first-generation counties K
  k_indices <- which(Pjk_vec >= threshold & !is.na(Pjk_vec))
  K <- names(k_indices)          # county names k
  
  # Output vector over districts
  SG_vec <- rep(NA_real_, length(Pji_vec))
  names(SG_vec) <- names(Pji_vec)   # district names
  
  # If no k above threshold: SG = direct
  if (length(K) == 0) {
    SG_vec <- Pji_vec
    return(SG_vec)
  }
  
  # Loop over districts i
  for (i in names(Pji_vec)) {
    Pji <- Pji_vec[i]   # j0 -> i
    
    # If no direct probability, keep NA
    if (is.na(Pji)) {
      SG_vec[i] <- NA_real_
      next
    }
    
    # (1 - Pji)
    term1 <- 1 - Pji
    
    # product over k of (1 - Pjk * Pki)
    prod_term <- 1
    for (k in K) {
      Pjk <- Pjk_vec[k]            # j0 -> k (county)
      Pki <- district_mat[i, k]    # k -> i (county k -> district i)
      
      if (is.na(Pjk) || is.na(Pki)) next
      prod_term <- prod_term * (1 - Pjk * Pki)
    }
    
    # Second-generation probability
    SG_vec[i] <- 1 - term1 * prod_term
  }
  
  return(SG_vec)
}








figure_MMR <- function() {
  
  # Plot A: County MMR
  p1 <- ggplot(infection_county) +
    geom_sf(aes(fill = MMR), color = "gray30", size = 0.1) +
    scale_fill_gradientn(
      colors = c("#d73027", "#fee08b", "#1a9850"),
      limits = c(0.6, 1.0),
      na.value = "lightgray",
      labels = scales::percent_format(accuracy = 1)
    ) +
    theme_minimal() +
    labs(
      title = "A: Weighted MMR Coverage by Texas County",
      fill  = "MMR"
    ) +
    theme_void() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 13),
      legend.position = "none"
    )
  
  # Plot B: District MMR
  p2 <- ggplot(infection_district) +
    geom_sf(aes(fill = MMR), color = "gray50", size = 0.1) +
    geom_sf(
      data  = tx_counties,
      fill  = NA,
      color = "gray20",
      size  = 0.4
    ) +
    scale_fill_gradientn(
      colors = c("#d73027", "#fee08b", "#1a9850"),
      limits = c(0.6, 1.0),
      na.value = "lightgray",
      labels   = scales::percent_format(accuracy = 1)
    ) +
    theme_minimal() +
    labs(
      title = "B: Weighted MMR Coverage by Texas School Districts",
      fill  = "MMR"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13)
    )
  
  legend   <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  
  row_plot   <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1, 1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2, 0.4))
  
  file_name <- "Figures/figure_B01.png"
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
}

figure_MMR()


figure_outbreaks_1st_2nd_county <- function(
    method    = 7,
    strategy  = 0,
    county    = "Gaines",
    threshold = 0.5,
    map_data  = infection_county,
    out_dir   = "Figures/"
) {
  
  
  # Define breaks and labels
  breaks     <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
  
  mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", strategy, ".rds"))
  
  county_upper <- toupper(county)
  county_names <- toupper(map_data$County)
  
  direct_vec <- mat[county_upper, match(county_names, colnames(mat))]
  
  indirect_vec <- compute_indirect_risk_county_new(
    trans_mat   = mat,
    county_name = county_upper,
    threshold   = threshold
  )
  
  highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
  
  # Cut values into bins
  map_data <- map_data %>%
    mutate(
      direct_bin = cut(
        direct_vec[match(toupper(County), names(direct_vec))],
        breaks = breaks,
        labels = bin_labels,
        include.lowest = TRUE,
        right = FALSE
      ),
      indirect_bin = cut(
        indirect_vec[match(toupper(County), names(indirect_vec))],
        breaks = breaks,
        labels = bin_labels,
        include.lowest = TRUE,
        right = FALSE
      )
    )
  
  # Plot A: Direct
  p1 <- ggplot(map_data) +
    geom_sf(aes(fill = direct_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "A: Probability of First-generation Outbreak") +
    theme_void() +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 13),
      legend.position = "none"
    )
  
  # Plot B: Indirect
  p2 <- ggplot(map_data) +
    geom_sf(aes(fill = indirect_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "B: Probability of First & Second-generation Outbreak") +
    theme_void() +
    theme(
      panel.grid  = element_blank(),
      axis.text   = element_blank(),
      axis.ticks  = element_blank(),
      axis.title  = element_blank()
    )
  
  legend   <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  
  row_plot   <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1, 1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2, 0.4))
  
  file_name <- paste0(out_dir)
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
}


figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "Gaines", out_dir="Figures/figure_B02_Gaines.png")
figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "KING", out_dir="Figures/figure_B02_KING.png")
figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "HALL", out_dir="Figures/figure_B02_HALL.png")
figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "THROCKMORTON", out_dir="Figures/figure_B02_THROCKMORTON.png")



figure_outbreaks_1st_2nd_district <- function(
    method    = 7,
    strategy  = 0,
    county    = "Gaines",       # outbreak source county j0
    threshold = 0.5,
    map_data  = infection_district,
    out_dir   = "Figures/"
) {
  
  
  # --- Load matrices ---
  district_mat <- readRDS(paste0("ProcessedData/district_pij_M", method, "_S", strategy, ".rds"))
  county_mat   <- readRDS(paste0("ProcessedData/county_pij_M",   method, "_S", strategy, ".rds"))
  
  county_upper <- toupper(county)
  
  # district names (map order)
  district_names <- map_data$district
  
  valid_districts <- rownames(district_mat)
  map_data <- map_data %>% filter(district %in% valid_districts)
  district_names <- map_data$district
  
  # DIRECT: P(j0 -> i)
  direct_vec <- district_mat[district_names, county_upper]
  
  # INDIRECT: first + second gen
  indirect_vec <- compute_indirect_risk_district_new(
    district_mat = district_mat,
    county_mat   = county_mat,
    county_name  = county_upper,
    threshold    = threshold
  )
  
  # --- Set breaks and color bins ---
  breaks     <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2","0.2–0.4","0.4–0.6","0.6–0.8","0.8–1.0")
  color_palette <- rev(c("#d73027","#fc8d59","#fee08b","#91cf60","#1a9850"))
  
  # --- Prepare map_data with bins ---
  map_data <- map_data %>%
    mutate(
      direct_bin = cut(
        direct_vec[match(district, names(direct_vec))],
        breaks = breaks, labels = bin_labels,
        include.lowest = TRUE, right = FALSE
      ),
      indirect_bin = cut(
        indirect_vec[match(district, names(indirect_vec))],
        breaks = breaks, labels = bin_labels,
        include.lowest = TRUE, right = FALSE
      )
    )
  
  
  # Highlight the source county’s districts
  highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
  
  # --- Plot A: Direct ---
  p1 <- ggplot(map_data) +
    geom_sf(aes(fill = direct_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "A: Probability of First-generation Outbreak") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none")
  
  # --- Plot B: Indirect ---
  p2 <- ggplot(map_data) +
    geom_sf(aes(fill = indirect_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "B: Probability of First & Second-generation Outbreak") +
    theme_void()
  
  legend <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  
  row_plot   <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1,1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2,0.4))
  
  file_name <- paste0(out_dir)
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
}


figure_outbreaks_1st_2nd_district(method=7, strategy=0, county="Gaines", threshold=0.5, out_dir="Figures/figure_B03_Gaines.png")
figure_outbreaks_1st_2nd_district(method=7, strategy=0, county="KING", threshold=0.5, out_dir="Figures/figure_B03_KING.png")
figure_outbreaks_1st_2nd_district(method=7, strategy=0, county="HALL", threshold=0.5, out_dir="Figures/figure_B03_HALL.png")
figure_outbreaks_1st_2nd_district(method=7, strategy=0, county="THROCKMORTON", threshold=0.5, out_dir="Figures/figure_B03_THROCKMORTON.png")


figure_strategies_county <- function(
    method = 7,
    counties = c("Gaines"),
    strategies = c(1, 2, 3),
    map_data = infection_county,
    out_dir = "Figures/"
) {
  
  color_palette <- c("white", "#1a9850")  # white to green
  
  county_names <- toupper(map_data$County)
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR < 95% → 95%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  # Load baseline
  baseline_mat <- readRDS(
    paste0("ProcessedData/county_pij_M", method, "_S0.rds")
  )
  
  for (county in counties) {
    
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>%
      filter(toupper(County) == county_upper)
    
    # Baseline vector
    baseline_vec <- baseline_mat[county_upper, match(county_names, colnames(baseline_mat))]
    
    # ----------------------------------------------------------
    # GLOBAL MAX across all strategies (for fixed color range)
    # ----------------------------------------------------------
    global_max <- max(baseline_vec, na.rm = TRUE)
    
    for (s in strategies) {
      mat_s <- readRDS(
        paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds")
      )
      vec_s <- mat_s[county_upper, match(county_names, colnames(mat_s))]
      global_max <- max(global_max, max(vec_s, na.rm = TRUE))
    }
    
    maps <- list()
    
    # ----------------------------------------------------------
    # BUILD ALL PANELS
    # ----------------------------------------------------------
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      
      mat_s <- readRDS(
        paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds")
      )
      vec_s <- mat_s[county_upper, match(county_names, colnames(mat_s))]
      
      map_data$outbreak <- vec_s[match(toupper(map_data$County), names(vec_s))]
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = outbreak), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom,
                fill = "blue", color = "black", size = 0.3) +
        scale_fill_gradient(
          low = "white",
          high = "#d73027",
          limits = c(0, global_max),     # FIXED RANGE ACROSS PANELS
          name = "Outbreak Probability",
          na.value = "gray80"
        ) +
        labs(title = paste0(label_prefix[i], ": ", strategy_labels[i])) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    # Legend from last plot
    legend <- get_legend(
      maps[[length(maps)]] +
        theme(legend.position = "right",
              legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(
      plotlist = lapply(maps, \(p) p + theme(legend.position = "none")),
      nrow = 1
    )
    
    final_plot <- plot_grid(
      row_plot, legend, nrow = 1,
      rel_widths = c(length(maps), 0.4)
    )
    
    file_name <- paste0(
      out_dir, "figure02_original_",
      tolower(county),
      "_method", method,
      "_S", paste(strategies, collapse = "_"),
      ".png"
    )
    
    ggsave(file_name, final_plot,
           width = 6 * length(maps),
           height = 6,
           dpi = 400)
  }
}

figure_strategies_county(method=7, counties=c("Gaines"), strategies=c(1, 2, 3), map_data=infection_county)





figure_strategies_district <- function(
    method     = 7,
    counties   = c("Gaines"),    # source counties to highlight
    strategies = c(1, 2, 3),
    map_data   = infection_district,
    out_dir    = "Figures/"
) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  color_palette <- c("white", "#d73027")   # white → red
  
  # -------------------------------------------------------
  # LOAD BASELINE FIRST — required for correct filtering
  # -------------------------------------------------------
  baseline_mat <- readRDS(
    paste0("ProcessedData/district_pij_M", method, "_S0.rds")
  )
  
  # -------------------------------------------------------
  # FIX: Remove districts missing from matrix
  # -------------------------------------------------------
  map_data <- map_data %>% 
    filter(district %in% rownames(baseline_mat))
  
  district_names <- map_data$district
  county_of_district <- toupper(map_data$County)
  
  # Labels
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR < 95% → 95%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix    <- LETTERS[seq_along(strategies)]
  
  for (county in counties) {
    
    county_upper <- toupper(county)
    
    # highlight all districts belonging to this county
    highlighted_geom <- map_data %>%
      filter(toupper(County) == county_upper)
    
    # baseline vector: P0(county → district)
    baseline_vec <- baseline_mat[district_names, county_upper]
    
    # -----------------------------------------------------
    # GLOBAL MAX across baseline and all strategies
    # -----------------------------------------------------
    global_max <- max(baseline_vec, na.rm = TRUE)
    
    for (s in strategies) {
      mat_s <- readRDS(
        paste0("ProcessedData/district_pij_M", method, "_S", s, ".rds")
      )
      vec_s <- mat_s[district_names, county_upper]
      global_max <- max(global_max, max(vec_s, na.rm = TRUE))
    }
    
    maps <- list()
    
    # -----------------------------------------------------
    # BUILD PANELS FOR EACH STRATEGY
    # -----------------------------------------------------
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      
      mat_s <- readRDS(
        paste0("ProcessedData/district_pij_M", method, "_S", s, ".rds")
      )
      vec_s <- mat_s[district_names, county_upper]
      
      # put values into map_data
      map_data$outbreak <- vec_s[match(map_data$district, names(vec_s))]
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = outbreak), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom,
                fill = "blue", color = "black", size = 0.3) +
        scale_fill_gradient(
          low  = "white",
          high = "#d73027",
          limits = c(0, global_max),
          name   = "Outbreak Probability",
          na.value = "gray80"
        ) +
        labs(title = paste0(label_prefix[i], ": ", strategy_labels[i])) +
        theme_void() +
        theme(plot.title = element_text(hjust = 0.5, size = 13),
              legend.position = "none")
      
      maps[[i]] <- p
    }
    
    # Legend from last plot
    legend <- get_legend(
      maps[[length(maps)]] +
        theme(legend.position = "right",
              legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(
      plotlist = lapply(maps, \(p) p + theme(legend.position = "none")),
      nrow = 1
    )
    
    final_plot <- plot_grid(
      row_plot, legend, nrow = 1,
      rel_widths = c(length(maps), 0.4)
    )
    
    file_name <- paste0(
      out_dir, "figureD02_original_",
      tolower(county),
      "_method", method,
      "_S", paste(strategies, collapse = "_"),
      ".png"
    )
    
    ggsave(file_name, final_plot,
           width = 6 * length(maps),
           height = 6,
           dpi = 400)
  }
}



figure_strategies_district(
  method = 7,
  counties = c("Gaines"),
  strategies = c(1, 2, 3),
  map_data = infection_district
)
