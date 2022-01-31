

# Cancer Resistance Selective Killing Simulation
# RSK.
# Freddie Whiting - 2022

# All plotting functions for VAF and 1/f comparisons. 

# Load libraries
library(plyr)
library(dplyr)
library(scales)
library(reldist)
library(tidyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(purrr)

# Need to change for plotting on cluster. 
options(bitmapType='cairo')

# Custom theme for plotting.

theme_FW <- function(text_size) {
  library(grid)
  library(ggthemes)
  theme_minimal() + 
    theme(panel.background = element_rect(colour = "grey80"),
          text = element_text(size = text_size), element_line(size = 0.8),
          panel.border = element_rect(colour = "black", fill=NA, size=2)
    )
}

####################
# Analysis Functions
####################


# Take a vector of raw mutation frequencies and convert into a table with: 
# [mutation relative frequency, count, 1/mutation relative frequency, 
# cumulative count]. 
# Can also be used on post-simulated sequencing VAF vectors. 

get_inv_mrf <- function(mut_rfs){
  
  # Get the counts of each relative frequency
  mut_rf_tbl <- data.frame(table(mut_rfs))
  colnames(mut_rf_tbl)[2] <- "count"
  mut_rf_tbl$mut_rfs <- as.numeric(as.character(mut_rf_tbl$mut_rfs))
  mut_rf_tbl <- mut_rf_tbl[nrow(mut_rf_tbl):1 ,]
  # Add 1/f column 
  mut_rf_tbl["inv_freq"] <- 1/mut_rf_tbl$mut_rfs
  # Add cumulative count column 
  mut_rf_tbl["cum_count"] <- cumsum(mut_rf_tbl$count)
  # Return the finished table
  return(mut_rf_tbl)
    
}




####################
# Plotting Functions
####################


# Plot the VAF distribution for: 
#   i) the original population 
#  ii) the population post genetic-bottleneck killing + subsequent expansion. 
# iii) the population post non-genetic bottleneck killing + subsequent
#      expansion. 

# Plot all as combined plot. Take histogram bin number and frequency of 
# 'resistant' allele (r) as arguments. 

comb_hist_plot <- function(mut_df, bin_n, min_freq, max_freq, USE_VAF=F, 
                           log_y=F){
  
  # Check only working with a single r and depth value: 
  sim_r <- unique(c(mut_df$r))
  sim_de <- unique(c(mut_df$depth))
  
  try(if(length(sim_r) != 1) stop("There should only be 1 unique r value in this dataframe."))
  try(if(length(sim_de) != 1) stop("There should only be 1 unique depth value in this dataframe."))
  
  # Subset the mutation dataframes by min/max freq, according to USE_VAF:
  if(USE_VAF==F){
    orig_df <- subset(mut_df, mut_rf_pre >= min_freq)
    btnk_gb_df <- subset(mut_df, mut_rf_gb >= min_freq)
    btnk_ngb_df <- subset(mut_df, mut_rf_ngb >= min_freq)
    # Set given frequency colname to 'x_freq' for plotting. 
    colnames(orig_df)[grep("mut_rf_pre", colnames(orig_df))] <- "x_freq"
    colnames(btnk_gb_df)[grep("mut_rf_gb", colnames(btnk_gb_df))] <-  "x_freq"
    colnames(btnk_ngb_df)[grep("mut_rf_ngb", colnames(btnk_ngb_df))] <- "x_freq"
    
  } else {
    orig_df <- subset(mut_df, pre_VAF >= min_freq)
    btnk_gb_df <- subset(mut_df, gb_VAF >= min_freq)
    btnk_ngb_df <- subset(mut_df, ngb_VAF >= min_freq)
    # Set given frequency colname to 'x_freq' for plotting. 
    colnames(orig_df)[grep("\\<pre_VAF", colnames(orig_df))] <- "x_freq"
    colnames(btnk_gb_df)[grep("\\<gb_VAF", colnames(btnk_gb_df))] <-  "x_freq"
    colnames(btnk_ngb_df)[grep("\\<ngb_VAF", colnames(btnk_ngb_df))] <- "x_freq"
  }
  
  p1 <- ggplot(data = orig_df, aes(x = x_freq)) + #, y = ..density..)) + 
    geom_histogram(bins = bin_n, alpha = 0.6, colour = "black") + 
    scale_x_continuous(limits = c(min_freq, max_freq)) + 
    ggtitle(paste0("r = ", sim_r, " depth = ", sim_de, ": \nOriginal VAF")) + 
    theme_FW(12)
  
  p2 <- ggplot(data = btnk_gb_df, aes(x = x_freq)) +#, y = ..density..)) +
    geom_histogram(bins = bin_n, alpha = 0.6, colour = "black") + 
    scale_x_continuous(limits = c(min_freq, max_freq)) + 
    ggtitle(paste0("r = ", sim_r, " depth = ", sim_de, ": \nGenetic-Bottleneck VAF")) + 
    theme_FW(12)
  
  p3 <- ggplot(data = btnk_ngb_df, aes(x = x_freq)) + #, y = ..density..)) + 
    geom_histogram(bins = bin_n, alpha = 0.6, colour = "black") + 
    scale_x_continuous(limits = c(min_freq, max_freq)) + 
    ggtitle(paste0("r = ", sim_r, " depth = ", sim_de, ": \nNon-Genetic-Bottleneck VAF")) + 
    theme_FW(12)
  
  if(log_y == T){
    p1 <- p1 + scale_y_log10()
    p2 <- p2 + scale_y_log10()
    p3 <- p3 + scale_y_log10()
  }

  if(USE_VAF==T){
    p1 <- p1 + xlab("VAF") + ylab("Count")
    p2 <- p2 + xlab("VAF") + ylab("Count")
    p3 <- p3 + xlab("VAF") + ylab("Count")
  } else {
    p1 <- p1 + xlab("Mutation Freq.") + ylab("Count")
    p2 <- p2 + xlab("Mutation Freq.") + ylab("Count")
    p3 <- p3 + xlab("Mutation Freq.") + ylab("Count")
  }
  
  fin_plot <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
  
  return(fin_plot)
  
}

# Plot an analagous plot of the 1/f vs cumulative number of mutations. 
# Can plot these all 

comb_1of_plot <- function(mut_df, min_freq, max_freq, ymax, USE_VAF=F, 
                          COMB_SIMS=F){
  
  # Check only working with a single r and depth value: 
  sim_r <- unique(c(mut_df$r))
  sim_de <- unique(c(mut_df$depth))
  sim_n <- unique(c(mut_df$Sim))
  # Only allow > 1 simualtion iterations if COMB_SIMS=T
  if(COMB_SIMS==F){
    try(if(length(sim_n) != 1) stop("There should only be a single simulation iteration when COMB_SIMS=F."))
  }
  
  try(if(length(sim_r) != 1) stop("There should only be 1 unique r value in this dataframe."))
  try(if(length(sim_de) != 1) stop("There should only be 1 unique depth value in this dataframe."))
  
  # Convert to 1/f tables, using columns according to whether USE_VAF=T
  if(COMB_SIMS==F){
    if(USE_VAF == F){
      orig_1of <- get_inv_mrf(mut_df$mut_rf_pre)
      btnk_gb_1of <- get_inv_mrf(mut_df$mut_rf_gb)
      btnk_ngb_1of <- get_inv_mrf(mut_df$mut_rf_ngb)
    } else {
      orig_1of <- get_inv_mrf(mut_df$pre_VAF)
      btnk_gb_1of <- get_inv_mrf(mut_df$gb_VAF)
      btnk_ngb_1of <- get_inv_mrf(mut_df$ngb_VAF)
    }  
  } else {
    if(USE_VAF == F){
      orig_1of <- plyr::ddply(mut_df, .(Sim), function(x){get_inv_mrf(x$mut_rf_pre)})
      btnk_gb_1of <- plyr::ddply(mut_df, .(Sim), function(x){get_inv_mrf(x$mut_rf_gb)})
      btnk_ngb_1of <- plyr::ddply(mut_df, .(Sim), function(x){get_inv_mrf(x$mut_rf_ngb)})
    } else {
      orig_1of <- plyr::ddply(mut_df, .(Sim), function(x){get_inv_mrf(x$pre_VAF)})
      btnk_gb_1of <- plyr::ddply(mut_df, .(Sim), function(x){get_inv_mrf(x$gb_VAF)})
      btnk_ngb_1of <- plyr::ddply(mut_df, .(Sim), function(x){get_inv_mrf(x$ngb_VAF)})
    }
  }
  
  if(COMB_SIMS==F){
    
    ggplot() + 
      geom_point(data = orig_1of, aes(x = inv_freq, y = cum_count, 
                                      fill = "original"),
                 shape = 21, colour = "black", alpha = 0.8, size = 2) + 
      geom_point(data = btnk_gb_1of, aes(x = inv_freq, y = cum_count, 
                                         fill = "genetic"),
                 shape = 21, colour = "black", alpha = 0.8, size = 2) + 
      geom_point(data = btnk_ngb_1of, aes(x = inv_freq, y = cum_count, 
                                          fill = "non-genetic"),
                 shape = 21, colour = "black", alpha = 0.8, size = 2) + 
      ggtitle(paste0("r = ", sim_r, " depth = ", sim_de, ": 1/f")) + 
      theme_FW(14) + 
      scale_colour_manual(values = c("navyblue", "red", "limegreen")) + 
      scale_x_continuous(limits = c(1/max_freq, 1/min_freq)) + 
      scale_y_continuous(limits = c(0, ymax)) + 
      xlab("1/f") + 
      ylab("M(f)")
    
  } else {
    
    ggplot() + 
      geom_point(data = orig_1of, aes(x = inv_freq, y = cum_count, 
                                      group = Sim, 
                                      fill = "original"),
                 shape = 21, colour = "black", alpha = 0.4, size = 2) + 
      geom_point(data = btnk_gb_1of, aes(x = inv_freq, y = cum_count, 
                                         group = Sim,
                                         fill = "genetic"),
                 shape = 21, colour = "black", alpha = 0.4, size = 2) + 
      geom_point(data = btnk_ngb_1of, aes(x = inv_freq, y = cum_count,
                                          group = Sim, 
                                          fill = "non-genetic"),
                 shape = 21, colour = "black", alpha = 0.4, size = 2) + 
      ggtitle(paste0("r = ", sim_r, " depth = ", sim_de, ": 1/f")) + 
      theme_FW(14) + 
      scale_colour_manual(values = c("navyblue", "red", "limegreen")) + 
      scale_x_continuous(limits = c(1/max_freq, 1/min_freq)) + 
      scale_y_continuous(limits = c(0, ymax)) + 
      xlab("1/f") + 
      ylab("M(f)")
    
  }
  
}


# Take a dataframe of simulation iterations and plot the combined VAF densities 
# as well as their combined mean and s.d. on the same plot. Return these plots
# for the original VAF, the genetic and non-genetic bottlenecks. 

# First, have a function that turns a list of dataframes and a given 
# VAF column into: 

# i) A combined density dataframe. 

comb_VAF_dens_df <- function(mut_df, VAF_col, min_VAF){
  
  # Split the mut_df by simulation ID. 
  split_dfs <- split(mut_df, mut_df$Sim)
  # Subset so all of the given VAFs are > min_VAF 
  split_dfs <- lapply(split_dfs, function(i){i[i[VAF_col] > min_VAF ,]})
  # Get the density for each.
  sim_denss <- lapply(split_dfs, function(i){density(i[[VAF_col]])})
  # Interploate the densities so they all have the same x coords
  approx_sim_denss <- lapply(sim_denss, approx, xout=seq(0.0, 1.0, by = 0.01))
  # Now create a matrix and populate with these values: 
  approx_sim_mat <- matrix(0, ncol = length(approx_sim_denss), nrow = length(approx_sim_denss[[1]]$y))
  
  for(i in 1:length(approx_sim_denss)){
    approx_sim_mat[,i] <- approx_sim_denss[[i]]$y
  }
  
  fin_df <- data.frame(approx_sim_mat)
  colnames(fin_df) <- paste0("Sim_", 1:ncol(fin_df))
  fin_df["x"] <- seq(0.0, 1.0, by = 0.01)
  
  # Turn from wide to long, by simulation number. 
  fin_df <- gather(fin_df, key = "Sim", value = y, grep("Sim_", colnames(fin_df)))
  
  return(fin_df)
  
}

# ii) A summary dataframe of the combined mean and sd VAF.

summ_VAF_dens_df <- function(mut_df, VAF_col, min_VAF){
  
  # Split the mut_df by simulation ID. 
  split_dfs <- split(mut_df, mut_df$Sim)
  # Subset so all of the given VAFs are > min_VAF 
  split_dfs <- lapply(split_dfs, function(i){i[i[VAF_col] > min_VAF ,]})
  # Get the density for each.
  sim_denss <- lapply(split_dfs, function(i){density(i[[VAF_col]])})
  # Interploate the densities so they all have the same x coords
  approx_sim_denss <- lapply(sim_denss, approx, xout=seq(0.0, 1.0, by = 0.01))
  # Now create a matrix and populate with these values: 
  approx_sim_mat <- matrix(0, ncol = length(approx_sim_denss), nrow = length(approx_sim_denss[[1]]$y))
  
  for(i in 1:length(approx_sim_denss)){
    approx_sim_mat[,i] <- approx_sim_denss[[i]]$y
  }
  
  # Now extract the mean and sd of density values for a given x value. 
  sim_mns <- rowMeans(approx_sim_mat, na.rm = T)
  sim_sds <- apply(approx_sim_mat, 1, sd, na.rm = T)
  
  # Now create an amalgamated dataframe with the: 
  # interpolated densities per simulation, 
  # amalgamated mean and sd for densitiies,
  # x axes used to standardies densitites, 
  
  summ_df <- data.frame(x = seq(0.0, 1.0, by = 0.01), 
                        sim_mn = sim_mns, sim_sd = sim_sds)
  
  return(summ_df)
  
}

# Have a function that plots the combined densities and summary (mean + sd) 
# density on the same panel: 

comb_VAF_dens_plot <- function(comb_df, summ_df){
  
  p <- ggplot(data = summ_df, aes(x = x, y = sim_mn)) + 
    geom_line(data = comb_df, aes(x = x, y = y, group = Sim), alpha = 0.7) +
    geom_ribbon(aes(ymin = sim_mn-sim_sd, ymax = sim_mn+sim_sd), alpha = 0.4, 
                colour = "black") +
    geom_line(size = 1.2, colour = "red") + 
    xlab("VAF") + 
    ylab("Density") +
    theme_FW(12)
  
  return(p)  
  
}


# 

all_comb_VAF_dens_plot <- function(mut_df, min_VAF, max_VAF){
  
  # Check only working with a single r and depth value. 
  sim_r <- unique(mut_df$r)
  sim_de <- unique(mut_df$depth)
  
  try(if(length(sim_r) != 1) stop("There should only be 1 unique r value in this dataframe."))
  try(if(length(sim_de) != 1) stop("There should only be 1 unique depth value in this dataframe."))
  
  # Get the combined and summary dataframes: 
  # Original VAF:  
  orig_comb_VAF_df <- comb_VAF_dens_df(mut_df, "pre_VAF", min_VAF)
  orig_summ_VAF_df <- summ_VAF_dens_df(mut_df, "pre_VAF", min_VAF)
  # Genetic bottleneck:
  btnk_comb_VAF_df_gb <- comb_VAF_dens_df(mut_df, "gb_VAF", min_VAF)
  btnk_summ_VAF_df_gb <- summ_VAF_dens_df(mut_df, "gb_VAF", min_VAF)
  # Non-genetic bottleneck:
  btnk_comb_VAF_df_ngb <- comb_VAF_dens_df(mut_df, "ngb_VAF", min_VAF)
  btnk_summ_VAF_df_ngb <- summ_VAF_dens_df(mut_df, "ngb_VAF", min_VAF)
  
  # Extract the maximum density to normalise plot axes: 
  omax <- max(orig_comb_VAF_df$y, na.rm = T)
  gmax <- max(btnk_comb_VAF_df_gb$y, na.rm = T)
  nmax <- max(btnk_comb_VAF_df_ngb$y, na.rm = T)
  ymax <- max(omax, gmax, nmax)
  
  # Create the plots: 
  p1 <- comb_VAF_dens_plot(orig_comb_VAF_df, orig_summ_VAF_df) + 
    scale_x_continuous(limits = c(min_VAF, max_VAF)) + 
    scale_y_continuous(limits = c(0, ymax+0.2)) + 
    ggtitle(paste0("r = ", sim_r, " depth = ", sim_de, ": Original VAF"))
  p2 <- comb_VAF_dens_plot(btnk_comb_VAF_df_gb, btnk_summ_VAF_df_gb) + 
    scale_x_continuous(limits = c(min_VAF, max_VAF)) + 
    scale_y_continuous(limits = c(0, ymax+0.2)) + 
    ggtitle(paste0("r = ", sim_r, " depth = ", sim_de, ": Genetic-Bottleneck VAF"))
  p3 <- comb_VAF_dens_plot(btnk_comb_VAF_df_ngb, btnk_summ_VAF_df_ngb) + 
    scale_x_continuous(limits = c(min_VAF, max_VAF)) + 
    scale_y_continuous(limits = c(0, ymax+0.2)) +
    ggtitle(paste0("r = ", sim_r, " depth = ", sim_de, ": Non-Genetic-Bottleneck VAF"))
  
  fin_p <- cowplot::plot_grid(p1, p2, p3, ncol = 1)
  
  return(fin_p)
  
  
}


################################################################################