library(ggplot2)
library(sf)
library(cowplot)
library(stringr)
library(dplyr)
library(tidyr)

map_probability <- readRDS("ProcessedData/map_probability.rds")



compute_indirect_risk <- function(
    trans_mat, 
    county_name = "GAINES", 
    threshold = .5
) {
  
  j <- toupper(county_name)  # source county j
  if (!(j %in% colnames(trans_mat))) stop("Source county not found in colnames(trans_mat).")
  
  i_names <- rownames(trans_mat)
  P_ij <- trans_mat[, j]   # j -> i
  
  k_all <- names(P_ij)[!is.na(P_ij) & P_ij >= threshold & names(P_ij) != j]
  
  P_sg <- rep(NA_real_, length(P_ij))
  names(P_sg) <- i_names
  
  
  for (i in i_names) {
    # enforce k != i from paper equation
    k_use <- setdiff(k_all, i)
    
    if (length(k_use) == 0) {
      p2 <- 0
    } else {
      vals <- trans_mat[i, k_use] * trans_mat[k_use, j]  # P_ik * P_kj
      p2   <- mean(vals, na.rm = TRUE)
    }  
    
    P_sg[i] <- min(1, P_ij[i] + p2) # combined first + second gen
  }
  
  P_sg[j] <- NA_real_  # self county not shown
  P_sg
}






############### 

figure01 <- function(
    method = 7, 
    strategy = 0, 
    county = "Gaines",
    threshold = 0.5, 
    map_data = map_probability,
    out_dir = "Figures/"
) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define breaks and labels
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0-0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))  # red → yellow → green
  
  mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", strategy, ".rds"))
  county_upper <- toupper(county)
  county_names <- toupper(map_data$County)
  
  direct_vec <- mat[match(county_names, colnames(mat)), county_upper]
  indirect_vec <- compute_indirect_risk(trans_mat = mat, county_name = county_upper, threshold = threshold)
  highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
  
  # Cut values into bins
  map_data <- map_data %>%
    mutate(
      direct_bin = cut(direct_vec[match(toupper(County), names(direct_vec))], breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE),
      indirect_bin = cut(indirect_vec[match(toupper(County), names(indirect_vec))], breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
    )
  
  # Plot A: Direct
  p1 <- ggplot(map_data) +
    geom_sf(aes(fill = direct_bin), color = "gray40", size = 0.1) +
    geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
    scale_fill_manual(values = color_palette, drop = FALSE, name = "Outbreak Probability") +
    labs(title = "A: Probability of First-generation Outbreak") +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 13),
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
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  legend <- get_legend(p2)
  p2_clean <- p2 + theme(legend.position = "none")
  row_plot <- plot_grid(p1, p2_clean, nrow = 1, rel_widths = c(1, 1))
  final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(2, 0.4))
  
  file_name <- paste0(out_dir, "figure01_outbreak_2ndGen.png")
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
  
}


figure01()

