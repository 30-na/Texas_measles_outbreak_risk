
library(cowplot)
library(ggplot2)
library(dplyr)
library(grid)
library(ggplot2)
library(cowplot)
library(dplyr)
library(grid)
library(sf)
############### Figure 01

figure01 <- function(method = 7, strategy = 0, county = "Gaines",
                                        threshold = 0.5, map_data = map_probability,
                                        out_dir = "Figures/") {

  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define breaks and labels
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))  # red → yellow → green
  
  mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", strategy, ".rds"))
  county_upper <- toupper(county)
  county_names <- toupper(map_data$County)
  
  direct_vec <- mat[county_upper, match(county_names, colnames(mat))]
  indirect_vec <- compute_indirect_risk(mat, county_name = county_upper, threshold = threshold)
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
  
  file_name <- paste0(out_dir, "figure01.png")
  ggsave(file_name, final_plot, width = 12, height = 6, dpi = 400)
  
}


figure01()




########################### Figure 02 ###################


figure02 <- function(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3),
                     map_data = map_probability, out_dir = "Figures/") {
  
  # Discrete bins
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))  # red → green
  names(color_palette) <- bin_labels
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  county_names <- toupper(map_data$County)
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds"))
      vec <- mat[county_upper, match(county_names, colnames(mat))]
      bin_vec <- cut(vec, breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
      bin_fact <- factor(bin_vec[match(toupper(map_data$County), names(vec))], levels = bin_labels)
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
    
    # Extract legend from the last map
    legend <- get_legend(
      maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
    
    file_name <- paste0(out_dir, "figure02_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(strategies), height = 6, dpi = 400)
  }
}



figure02(method = 7, map_probability, counties = c("Gaines"), strategies = c(1, 2, 3))





####################### Figure 03
figure03 <- function(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3),
                     map_data = map_probability, threshold = 0.5, out_dir = "Figures/") {
  
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
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  county_names <- toupper(map_data$County)
  mats <- lapply(strategies, function(s) readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds")))
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- mats[[i]]
      indirect_vec <- compute_indirect_risk(mat, county_name = county_upper, threshold = threshold)
      
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
      maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.4))
    
    file_name <- paste0(out_dir, "figure03_Indirect_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(strategies), height = 6, dpi = 400)
  }
}



figure03(method = 7, map_probability, counties = c("Gaines"), strategies = c(1, 2, 3),threshold = .5)



############################ Figure04
figure04 <- function(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3),
                     map_data = map_probability, out_dir = "Figures/") {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  color_palette <- c("white", "#1a9850")  # white to green
  county_names <- toupper(map_data$County)
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  baseline_mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S0.rds"))
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    baseline_vec <- baseline_mat[county_upper, match(county_names, colnames(baseline_mat))]
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds"))
      vec <- mat[county_upper, match(county_names, colnames(mat))]
      reduction <- baseline_vec - vec
      reduction[reduction < 0] <- 0
      
      map_data$reduction <- reduction[match(toupper(map_data$County), names(reduction))]
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = reduction), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        scale_fill_gradient(
          low = "white", high = "#1a9850", limits = c(0, max(reduction, na.rm = TRUE)),
          name = "Reduction in\nOutbreak Probability", na.value = "gray80"
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
    
    file_name <- paste0(out_dir, "figure04_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(maps), height = 6, dpi = 400)
  }
}



figure04(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3), map_data = map_probability)



figure05 <- function(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3),
                     map_data = map_probability, threshold = 0.5, out_dir = "Figures/") {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  color_palette <- c("white", "#1a9850")  # white to green
  county_names <- toupper(map_data$County)
  
  full_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )
  
  strategy_labels <- full_labels[strategies + 1]
  label_prefix <- LETTERS[seq_along(strategies)]
  
  baseline_mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S0.rds"))
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    baseline_vec <- compute_indirect_risk(baseline_mat, county_name = county_upper, threshold = threshold)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds"))
      vec <- compute_indirect_risk(mat, county_name = county_upper, threshold = threshold)
      reduction <- baseline_vec - vec
      reduction[reduction < 0] <- 0
      
      map_data$reduction <- reduction[match(toupper(map_data$County), names(reduction))]
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = reduction), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        scale_fill_gradient(
          low = "white", high = "#1a9850", limits = c(0, max(reduction, na.rm = TRUE)),
          name = "Reduction in\nOutbreak Probability", na.value = "gray80"
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
    
    file_name <- paste0(out_dir, "figure05_Indirect_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(maps), height = 6, dpi = 400)
  }
}



figure05(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3), map_data = map_probability, threshold = 0.5)


############ Figure06
figure06 <- function(map_data = map_probability, out_dir = "Figures/") {
  
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
  
  file_name <- paste0(out_dir, "figure06_mmr_and_increase_grid.png")
  ggsave(file_name, grid, width = 12, height = 10, dpi = 400)
  
  return(grid)
}



figure06(map_data = map_probability)























#################################### TEST



figure02 <- function(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3),
                     map_data = map_probability, out_dir = "Figures/") {
  
  # Discrete bins and color palette
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))
  names(color_palette) <- bin_labels
  strategy_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )[strategies + 1]
  
  county_names <- toupper(map_data$County)
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds"))
      vec <- mat[county_upper, match(county_names, colnames(mat))]
      bin_vec <- cut(vec, breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
      bin_fact <- factor(bin_vec[match(toupper(map_data$County), names(vec))], levels = bin_labels)
      map_data$bin <- bin_fact
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = bin), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        scale_fill_manual(values = color_palette, limits = bin_labels, name = "Outbreak Probability") +
        labs(title = strategy_labels[i]) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    # Create dummy map to extract legend (5 groups + 1 NA)
    dummy_geom <- st_sf(
      bin = factor(c(bin_labels, NA), levels = bin_labels),
      geometry = st_sfc(
        st_polygon(list(rbind(c(0,0), c(1,0), c(1,1), c(0,1), c(0,0)))),
        st_polygon(list(rbind(c(1,0), c(2,0), c(2,1), c(1,1), c(1,0)))),
        st_polygon(list(rbind(c(2,0), c(3,0), c(3,1), c(2,1), c(2,0)))),
        st_polygon(list(rbind(c(3,0), c(4,0), c(4,1), c(3,1), c(3,0)))),
        st_polygon(list(rbind(c(4,0), c(5,0), c(5,1), c(4,1), c(4,0)))),
        st_polygon(list(rbind(c(5,0), c(6,0), c(6,1), c(5,1), c(5,0))))
      )
    )
    
    p_dummy <- ggplot(dummy_geom) +
      geom_sf(aes(fill = bin), color = "black") +
      scale_fill_manual(
        values = color_palette,
        limits = bin_labels,
        drop = FALSE,
        na.value = "gray80",
        name = "Outbreak Probability"
      ) +
      theme_void() +
      theme(
        legend.position = "right",
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.key.size = unit(2, "lines"),         # controls both width & height (if key.width/key.height not set)
        legend.key.height = unit(2, "lines"),       # vertical size of color boxes
        legend.key.width = unit(2, "lines"),        # horizontal size of color boxes
        legend.spacing.y = unit(0.5, "lines")
      )
    
    legend <- get_legend(p_dummy)
    
    row_plot <- plot_grid(plotlist = maps, nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.5))
    
    file_name <- paste0(out_dir, "figure02_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(strategies), height = 6, dpi = 400)
  }
}



figure02(method = 7, map_probability, counties = c("Gaines"), strategies = c(1, 2, 3))

##################################

figure02 <- function(method = 7, counties = c("Gaines"), strategies = c(1, 2, 3),
                     map_data = map_probability, out_dir = "Figures/") {
  
  
  # Discrete bins
  breaks <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  bin_labels <- c("0–0.2", "0.2–0.4", "0.4–0.6", "0.6–0.8", "0.8–1.0")
  color_palette <- rev(c("#d73027", "#fc8d59", "#fee08b", "#91cf60", "#1a9850"))  # red → green
  names(color_palette) <- bin_labels
  strategy_labels <- c(
    "Strategy 0: Baseline",
    "Strategy 1: MMR < 90% → 90%", 
    "Strategy 2: MMR < 92% → 92%", 
    "Strategy 3: MMR +5%",
    "Strategy 4: Gaines → 90%",
    "Strategy 5: Gaines → 92%"
  )[strategies + 1]
  
  county_names <- toupper(map_data$County)
  
  for (county in counties) {
    county_upper <- toupper(county)
    highlighted_geom <- map_data %>% filter(toupper(County) == county_upper)
    
    maps <- list()
    
    for (i in seq_along(strategies)) {
      s <- strategies[i]
      mat <- readRDS(paste0("ProcessedData/pij_M", method, "_S", s, ".rds"))
      vec <- mat[county_upper, match(county_names, colnames(mat))]
      bin_vec <- cut(vec, breaks = breaks, labels = bin_labels, include.lowest = TRUE, right = FALSE)
      bin_fact <- factor(bin_vec[match(toupper(map_data$County), names(vec))], levels = bin_labels)
      map_data$bin <- bin_fact
      
      
      p <- ggplot(map_data) +
        geom_sf(aes(fill = bin), color = "gray40", size = 0.1) +
        geom_sf(data = highlighted_geom, fill = "blue", color = "black", size = 0.3) +
        
        scale_fill_manual(values = color_palette, limits = bin_labels, name = "Outbreak Probability")+
        
        labs(title = strategy_labels[i]) +
        theme_void() +
        theme(
          plot.title = element_text(hjust = 0.5, size = 13),
          legend.position = "none"
        )
      
      maps[[i]] <- p
    }
    
    # Extract legend from the last map
    legend <- get_legend(
      maps[[2]] + theme(legend.position = "right", legend.title = element_text(size = 10))
    )
    
    row_plot <- plot_grid(plotlist = lapply(maps, \(p) p + theme(legend.position = "none")), nrow = 1)
    final_plot <- plot_grid(row_plot, legend, nrow = 1, rel_widths = c(length(maps), 0.3))
    
    file_name <- paste0(out_dir, "figure02_", tolower(county), "_method", method, "_S", paste(strategies, collapse = "_"), ".png")
    ggsave(file_name, final_plot, width = 6 * length(strategies), height = 6, dpi = 400)
    
  }
}


figure02(method = 7, map_probability, counties = c("Gaines"), strategies = c(1, 2, 3))
