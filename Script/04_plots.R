library(cowplot)
library(ggplot2)
library(dplyr)
library(grid)
library(sf)

# Load datasets ----
infection_county   <- readRDS("ProcessedData/map_county_infection_proportion_ratio0.5_efficacy01.rds")


tx_counties <- counties(
  state = "TX",
  year  = 2024
) %>%
  mutate(county = toupper(NAME)) %>%
  select(county)


# Functions -----

compute_indirect_risk_county_new <- function(
    trans_mat,
    county_name, 
    threshold = -1
    ) {
  # Paper notation:
  # P_{ij} = trans_mat[i, j]
  # i = destination
  # j = origin (index county)
  
  j <- toupper(county_name) # origin county
  counties <- rownames(trans_mat)
  
  if (!(j %in% counties))
    stop("County j not found in matrix.")
  
  # Extract the j-th column: P_{ij} for all i
  Pij <- trans_mat[, j]
  
  # Output vector
  P_SG <- rep(NA_real_, length(Pij))
  names(P_SG) <- counties
  
  # Loop over destination counties i
  for (i in counties) {
    
    P_ij <- Pij[i]   # direct probability j → i
    
    if (is.na(P_ij)) {
      P_SG[i] <- NA
      next
    }
    
    # Identify intermediates k county
    # intermediates_county = all k such that k != j, and P_{kj} >= threshold
    Pkj <- trans_mat[, j]          # P_{kj} = j → k
    #intermediates_county <- names(Pkj)[!is.na(Pkj) & Pkj >= threshold & names(Pkj) != j]
    intermediates_county <- names(Pkj)[!is.na(Pkj) & names(Pkj) != j]
    # Compute the product term
    prod_term <- 1
    for (k in intermediates_county) {
      P_ik <- trans_mat[i, k]   # k → i
      P_kj <- trans_mat[k, j]   # j → k
      if (is.na(P_ik)) next
      prod_term <- prod_term * (1 - P_ik * P_kj)
    }
    
    # Paper formula:
    # P^{SG}_{ij} = 1 - (1 - P_{ij}) * product_k (1 - P_{ik} P_{kj})
    P_SG[i] <- 1 - (1 - P_ij) * prod_term
  }
  
  # Self infection is undefined
  P_SG[j] <- NA
  
  return(P_SG)
}



# Figure 01 ------

figure_outbreaks_1st_2nd_county <- function(
    method    = 7,
    strategy  = 0,
    county    = "Gaines",
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
  
  direct_vec <- mat[match(county_names, rownames(mat)), county_upper]
  
  indirect_vec <- compute_indirect_risk_county_new(
    trans_mat   = mat,
    county_name = county_upper
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
    labs(title = "A: Probability of First generation Outbreak") +
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
    labs(title = "B: Probability of First & Second generation Outbreak") +
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


figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "Gaines",       out_dir="Figures/figure_01_Gaines.png")
# figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "KING",         out_dir="Figures/figure_01_KING.png")
# figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "HALL",         out_dir="Figures/figure_01_HALL.png")
# figure_outbreaks_1st_2nd_county(method = 7, strategy = 0, county = "THROCKMORTON", out_dir="Figures/figure_01_THROCKMORTON.png")



# Figure 02 -----

figure_reduction_strategies <- function(
    method = 7,
    counties = c("Gaines"), 
    strategies = c(1, 2, 3),
    map_data = infection_county, 
    out_dir = "Figures/"
    ) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  color_palette <- c("white", "#1a9850")  # white to green
  county_names <- toupper(map_data$County)
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR < 95% → 95%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  baseline_mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S0.rds"))
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    baseline_vec <- baseline_mat[ match(county_names, rownames(baseline_mat)), county_upper]
    
    maps <- list()
    
    # Compute global max reduction across all strategies for this county
    global_max <- 0
    for (s in strategies) {
      mat_s <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds"))
      vec_s <- mat_s[match(county_names, rownames(mat_s)), county_upper]
      red_s <- (baseline_vec - vec_s)*100
      global_max <- max(global_max, max(red_s, na.rm = TRUE))
    }
    
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds"))
      vec <- mat[match(county_names, rownames(mat)), county_upper]
      reduction <- (baseline_vec - vec)*100
      #reduction[reduction < 0] <- 0
      
      map_data$reduction <- reduction[match(toupper(map_data$County), names(reduction))]
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = reduction), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        scale_fill_gradient(
          low = "white",
          high = "#1a9850", 
          #limits = c(0, max(reduction, na.rm = TRUE)),
          limits = c(0, global_max),
          name = "Reduction in\nOutbreak Probability (%)", na.value = "gray80"
        ) +
        labs(title = paste0(label_prefix[i], ": ", strategy_labels[i])) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    legend <- get_legend(
      maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
    
    file_name <- paste0(out_dir, "figure_02_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(maps), height = 6, dpi = 400)
  }
}

figure_reduction_strategies(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3), map_data = infection_county)



# Figure 03 -----

Figure_indirect_transmission_multiple_county <- function(
    method,
    map_data, 
    counties,
    strategy = 0, 
    out_dir = "Figures/"
    ) {
  
  # Define breaks and labels
  breaks     <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
  
  
  county_names <- toupper(map_data$County)
  
  mat <- readRDS(paste0("ProcessedData/county_pij_M", method, "_S", strategy, ".rds"))
  maps <- list()
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    indirect_vec <- compute_indirect_risk_county_new(mat, county_name = county_upper)
    
    map_data$indirect_risk <- indirect_vec[match(county_names, names(indirect_vec))]
    
    # Cut values into bins
    map_data <- map_data %>%
      mutate(
        indirect_bin = cut(
          indirect_vec[match(toupper(County), names(indirect_vec))],
          breaks = breaks,
          labels = bin_labels,
          include.lowest = TRUE,
          right = FALSE
        )
      )
    
    p <- ggplot(map_data) +
      geom_sf(aes(fill = indirect_bin), color = "gray40", size = 0.1) +
      geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
      scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
      labs(
        title = str_to_title(county),
        fill = "Outbreak Probability"
      ) +
      theme_minimal() +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = "none"
      )
    
    maps[[length(maps) + 1]] <- p
  }
  
  legend <- get_legend(
    maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
  )
  
  row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
  
  file_name <- paste0(out_dir, "figure_03_indirect_risk_multiple_method", method,
                      "_S", strategy, "_", paste(tolower(counties), collapse = "_"), ".png")
  ggsave(file_name, final_plot, width = 6 * length(counties), height = 6, dpi = 400)
  
}



Figure_indirect_transmission_multiple_county(
  method = 7,
  map_data = infection_county,
  counties = c("POLK", "MONTAGUE", "LIMESTONE"),
  #counties = c("KING", "HALL", "THROCKMORTON"),
  #counties = c("DALLAS", "TARRANT", "COLLIN"),
  strategy = 0
)







# Figure S01 -----

figure_MMR_strategies <- function(
    map_data = map_probability,
    out_dir = "Figures/"
    ) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  color_mmr <- c("#d73027", "#fee08b", "#1a9850")     # for baseline MMR
  color_increase <- c("white", "#1a9850")             # for MMR increases
  
  # Compute increases (clamped)
  map_data <- map_data %>%
    mutate(
      delta_mmr1 = pmin(1, pmax(0, MMR1 - MMR)),
      delta_mmr2 = pmin(1, pmax(0, MMR2 - MMR)),
      delta_mmr3 = pmin(1, pmax(0, MMR3 - MMR))
    )
  
  plots <- list(
    ggplot(map_data) +
      geom_sf(aes(fill = MMR), color = "gray40", size = 0.1) +
      scale_fill_gradientn(colors = color_mmr, limits = c(0.8, 1), name = "MMR") +
      labs(title = "A: Baseline MMR") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 13)),
    
    ggplot(map_data) +
      geom_sf(aes(fill = delta_mmr1), color = "gray40", size = 0.1) +
      scale_fill_gradientn(colors = color_increase, limits = c(0, 0.2), name = "Increase in MMR") +
      labs(title = "B: MMR Vaccine Rate Increase by Strategy 01") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 13)),
    
    ggplot(map_data) +
      geom_sf(aes(fill = delta_mmr2), color = "gray40", size = 0.1) +
      scale_fill_gradientn(colors = color_increase, limits = c(0, 0.2), name = "Increase in MMR") +
      labs(title = "C: MMR Vaccine Rate Increase by Strategy 02") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 13)),
    
    ggplot(map_data) +
      geom_sf(aes(fill = delta_mmr3), color = "gray40", size = 0.1) +
      scale_fill_gradientn(colors = color_increase, limits = c(0, 0.2), name = "Increase in MMR") +
      labs(title = "D: MMR Vaccine Rate Increase by Strategy 03") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, size = 13))
  )
  
  grid <- plot_grid(plotlist = plots, ncol = 2, align = "hv")
  
  file_name <- paste0(out_dir, "figure_S01_MMR_increase_strategies.png")
  ggsave(file_name, grid, width = 12, height = 10, dpi = 400)
  
  return(grid)
}

figure_MMR_strategies(map_data = infection_county)



# Figure S02 ------
figure_second_generation_strategies <- function(
    method = 7,
    counties = c("Gaines"), 
    strategies = c(1, 2, 3),
    map_data = map_probability, 
    out_dir = "Figures/"
    ) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Discrete bins
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
  names(color_palette) <- bin_labels
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR < 95% → 95%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  county_names <- toupper(map_data$County)
  mats <- lapply(strategies, function(s) readRDS(paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds")))
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- mats[[i]]
      indirect_vec <- compute_indirect_risk_county_new(mat, county_name = county_upper)
      
      bin_vec <- cut(indirect_vec, breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
      bin_fact <- factor(bin_vec[match(toupper(map_data$County), names(indirect_vec))], levels = bin_labels)
      map_data$bin <- bin_fact
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = bin), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        scale_fill_manual(
          values = color_palette,
          name = "Outbreak Probability",
          drop = FALSE
        ) +
        labs(title = paste0(label_prefix[i], ": ", strategy_labels[i])) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    legend <- get_legend(
      maps[[1]] + theme(legend.position = "right", legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
    
    file_name <- paste0(out_dir, "figure_S02_Indirect_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(strategies), height = 6, dpi = 400)
  }
}

figure_second_generation_strategies(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3) ,map_data=infection_county)



























# 
# 
# 
# figure_strategies_county <- function(
#     method = 7,
#     counties = c("Gaines"),
#     strategies = c(1, 2, 3),
#     map_data = infection_county,
#     out_dir = "Figures/"
# ) {
#   
#   color_palette <- c("white", "#1a9850")  # white to green
#   
#   county_names <- toupper(map_data$County)
#   
#   full_labels <- c(
#     "Strategy 0: Baseline",
#     "Strategy 1: MMR < 90% → 90%", 
#     "Strategy 2: MMR < 92% → 92%", 
#     "Strategy 3: MMR < 95% → 95%"
#   )
#   
#   strategy_labels <- full_labels[strategies + 1]
#   label_prefix <- LETTERS[seq_along(strategies)]
#   
#   # Load baseline
#   baseline_mat <- readRDS(
#     paste0("ProcessedData/county_pij_M", method, "_S0.rds")
#   )
#   
#   for (county in counties) {
#     
#     county_upper <- toupper(county)
#     highlighted_geom <- map_data %>%
#       filter(toupper(County) == county_upper)
#     
#     # Baseline vector
#     baseline_vec <- baseline_mat[county_upper, match(county_names, colnames(baseline_mat))]
#     
#     # ----------------------------------------------------------
#     # GLOBAL MAX across all strategies (for fixed color range)
#     # ----------------------------------------------------------
#     global_max <- max(baseline_vec, na.rm = TRUE)
#     
#     for (s in strategies) {
#       mat_s <- readRDS(
#         paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds")
#       )
#       vec_s <- mat_s[county_upper, match(county_names, colnames(mat_s))]
#       global_max <- max(global_max, max(vec_s, na.rm = TRUE))
#     }
#     
#     maps <- list()
#     
#     # ----------------------------------------------------------
#     # BUILD ALL PANELS
#     # ----------------------------------------------------------
#     for (i in seq_along(strategies)) {
#       s <- strategies[i]
#       
#       mat_s <- readRDS(
#         paste0("ProcessedData/county_pij_M", method, "_S", s, ".rds")
#       )
#       vec_s <- mat_s[county_upper, match(county_names, colnames(mat_s))]
#       
#       map_data$outbreak <- vec_s[match(toupper(map_data$County), names(vec_s))]
#       
#       p <- ggplot(map_data) +
#         geom_sf(aes(fill = outbreak), color = "gray40", size = 0.1) +
#         geom_sf(data = highlighted_geom,
#                 fill = "blue", color = "black", size = 0.3) +
#         scale_fill_gradient(
#           low = "white",
#           high = "#d73027",
#           limits = c(0, global_max),     # FIXED RANGE ACROSS PANELS
#           name = "Outbreak Probability",
#           na.value = "gray80"
#         ) +
#         labs(title = paste0(label_prefix[i], ": ", strategy_labels[i])) +
#         theme_void() +
#         theme(
#           plot.title = element_text(hjust = 0.5, size = 13),
#           legend.position = "none"
#         )
#       
#       maps[[i]] <- p
#     }
#     
#     # Legend from last plot
#     legend <- get_legend(
#       maps[[length(maps)]] +
#         theme(legend.position = "right",
#               legend.title = element_text(size = 10))
#     )
#     
#     row_plot <- plot_grid(
#       plotlist = lapply(maps, \(p) p + theme(legend.position = "none")),
#       nrow = 1
#     )
#     
#     final_plot <- plot_grid(
#       row_plot, legend, nrow = 1,
#       rel_widths = c(length(maps), 0.4)
#     )
#     
#     file_name <- paste0(
#       out_dir, "figure02_original_",
#       tolower(county),
#       "_method", method,
#       "_S", paste(strategies, collapse = "_"),
#       ".png"
#     )
#     
#     ggsave(file_name, final_plot,
#            width = 6 * length(maps),
#            height = 6,
#            dpi = 400)
#   }
# }
# 
# figure_strategies_county(method=7, counties=c("Gaines"), strategies=c(1, 2, 3), map_data=infection_county)
# 
# 
