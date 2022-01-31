

# Cancer Resistance Selective Killing Simulation
# RSK.
# Freddie Whiting - 2021

rm(list = ls())

# Load libraries
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


################################################################################


# Load in command line arguments: 

args <- commandArgs(trailingOnly = TRUE)

# ARG1 : output directory 

out_dir <- args[1]

# ARG2: plotting script location

plot_scr <- args[2]

# Need to already be in the Outputs directory.

setwd(out_dir)

source(plot_scr)

# Get the mutation and (if there are any) simulated VAF dataframes
mut_dfs <- grep("mut_df", list.files(), value = T)
VAF_dfs <- grep("VAF_df", list.files(), value = T)  

# Make plots directory if it doesn't exist. 
if(dir.exists("Plots")==F){
  dir.create("Plots")
}


# Functions to add a column to each output dataframe that distinguishes 
# simulation iteration

pop_out_df <- function(out_csv){
  
  out_df <- read.csv(out_csv, stringsAsFactors = F)
  sim_n <- sub(".*_(\\d{1,}).csv", out_csv, replacement = "\\1")
  out_df["Sim"] <- sim_n
  return(out_df)

}

# If there are VAF dfs, read in these and use their mutation frequencies, 
# otherwise read in the mutation dfs. 

if(length(VAF_dfs) > 0){
  VAF_pres <- TRUE
}
  

if(VAF_pres == T){
  
  out_dfs <- lapply(VAF_dfs, pop_out_df)
  
} else {
  
  out_dfs <- lapply(VAF_dfs, pop_out_df)
  
}

# Merge the output dataframes. 

out_df <- bind_rows(out_dfs)


# Want to change the min and max VAF depending on the depth used. Create a 
# dataframe that controls the min and max VAF for a given depth.  

depth_df <- data.frame(depth = c(10, 100, 1000), 
                       min_VAF = c(0.1, 0.1, 0.01),
                       max_VAF = c(0.3, 0.3, 0.3),
                       ymax = c(10, 100, 1000))


##############################################################
# Functions to save the plots given the chosen plot parameters. 
##############################################################

# Combined histogram plot. 

save_comb_hist_plot <- function(mut_df, bin_n, min_freq, max_freq, 
                                USE_VAF=F, log_y=F){
  
  setwd("Plots")
  # Check plotting parameter types.
  try(if(typeof(USE_VAF) != "logical") stop("USE_VAF should be TRUE or FALSE."))
  try(if(typeof(log_y) != "logical") stop("log_y should be TRUE or FALSE."))
  
  # Extract the simulation r and depth values. 
  sim_r <- as.numeric(unique(mut_df$r))
  sim_de <- as.numeric(unique(mut_df$depth))
  sim_n <- as.numeric(unique(mut_df$Sim))
  # Should only have one simulation for this plot. 
  try(if(length(sim_n) > 1) stop("There should only be a single simulation output in this DataFrame."))
  
  chp <- comb_hist_plot(mut_df, bin_n, min_freq, max_freq, 
                        USE_VAF=USE_VAF, log_y=log_y)
  
  if(USE_VAF == T){
    plot_name <- paste0("Sim_", sim_n, "_r-", sim_r, "_de-", sim_de,
                        "_VAF_chplot.jpeg")
  } else {
    plot_name <- paste0("Sim_", sim_n, "_r-", sim_r, "_de-", sim_de,
                        "_mut_chplot.jpeg")
  }
  
  # Save the plot.
  ggsave(plot_name, chp, width = 6, height = 12)
  setwd("../")
}

# Combined 1/f plot. 

save_comb_1of_plot <- function(mut_df, min_freq, max_freq, ymax, 
                               USE_VAF=F, COMB_SIMS=F){
  
  setwd("Plots")
  # Check plotting parameter types.
  try(if(typeof(USE_VAF) != "logical") stop("USE_VAF should be TRUE or FALSE."))
  try(if(typeof(COMB_SIMS) != "logical") stop("COMB_SIMS should be TRUE or FALSE."))
  
  # Extract the simulation r and depth values. 
  sim_r <- as.numeric(unique(mut_df$r))
  sim_de <- as.numeric(unique(mut_df$depth))
  sim_n <- as.numeric(unique(mut_df$Sim))
  # If COMB_SIMS is FALSE, should only have one simulation for this plot.
  if(COMB_SIMS == F){
    try(if(length(sim_n) > 1) stop("There should only be a single simulation output in this dataframe."))  
  }
  
  cfp <- comb_1of_plot(mut_df, min_freq, max_freq, ymax, USE_VAF=USE_VAF,
                       COMB_SIMS=COMB_SIMS)
  
  if(COMB_SIMS == F){
    if(USE_VAF == T){
      plot_name <- paste0("Sim_", sim_n, "_r-", sim_r, "_de-", sim_de,
                          "_VAF_cfplot.jpeg")
    } else {
      plot_name <- paste0("Sim_", sim_n, "_r-", sim_r, "_de-", sim_de,
                          "_mut_cfplot.jpeg")
    }

  } else {
    if(USE_VAF == T){
      plot_name <- paste0("All_Sims_r-", sim_r, "_de-", sim_de, 
                          "_VAF_cfplot.jpeg")
    } else {
      plot_name <- paste0("All_Sims_r-", sim_r, "_de-", sim_de, 
                          "_mut_cfplot.jpeg")
    }
  }
  
  # Save the plot. 
  ggsave(plot_name, cfp, width = 8, height = 6)
  setwd("../")
  

}


# Combined VAF density plot

save_all_comb_VAF_dens_plot <- function(mut_df, min_VAF, max_VAF){
  
  setwd("Plots")
  # Extract the simulation r and depth values. 
  sim_r <- as.numeric(unique(mut_df$r))
  sim_de <- as.numeric(unique(mut_df$depth))
  
  # There should be multiple simulation iterations for this plot. 
  sim_n <- as.numeric(unique(mut_df$Sim))
  if(length(sim_n) == 1){
    warning("There is only a single simulation iteration in this dataframe.")
  }
  
  avp <- all_comb_VAF_dens_plot(mut_df, min_VAF, max_VAF)
  
  plot_name <- paste0("All_Sims_r-", sim_r, "_de-", sim_de, "_VAF_avplot.jpeg")
  
  # Save the plot. 
  ggsave(plot_name, avp, width = 7, height = 16)
  setwd("../")
  
}
  
################################################################################

# Split the output dataframes, first by r, and then by depth. 

split_dfs <- split(out_df, f = out_df$r)

split_dfs <- lapply(split_dfs, function(x){split(x, f = x$depth)})

# Go through each parameter combination and save the plots. 

for(i in seq_along(split_dfs)){
  
  for(j in seq_along(split_dfs[[i]])){
    
    # Simulation depth and r...
    sim_r <- as.numeric(unique(split_dfs[[i]][[j]]$r))
    sim_de <- as.numeric(unique(split_dfs[[i]][[j]]$depth))
    # min_VAF, max_VAF and ymax for plotting given the depth...
    de_min_VAF <- depth_df[depth_df$depth == sim_de ,]$min_VAF
    de_max_VAF <- depth_df[depth_df$depth == sim_de ,]$max_VAF
    de_ymax <- depth_df[depth_df$depth == sim_de ,]$ymax
    
    
    if(VAF_pres==T){
      
      # Histograms:
      # - raw mutation frequencies for each sim. 
      plyr::ddply(split_dfs[[i]][[j]], .(Sim), function(x){
        save_comb_hist_plot(x, 100, 0, 1.01, USE_VAF=F, log_y=T)
      })
      # - simulated VAFs for each sim. 
      plyr::ddply(split_dfs[[i]][[j]], .(Sim), function(x){
        save_comb_hist_plot(x, 100, 0.05, 0.7, USE_VAF=T, log_y=F)
      })
      
      # 1/f plots:
      # - raw mutation frequencies for all sims. 
      save_comb_1of_plot(split_dfs[[i]][[j]], 0.001, 1.0, 10000, 
                         USE_VAF=F, COMB_SIMS=T)
      # - simulated VAFs for all sims. 
      save_comb_1of_plot(split_dfs[[i]][[j]], de_min_VAF, de_max_VAF, de_ymax, 
                         USE_VAF=T, COMB_SIMS=T)
      
      # Combined VAF densities:
      # - simulated VAFs
      save_all_comb_VAF_dens_plot(split_dfs[[i]][[j]], 0.1, 0.7)

        
    } else {
      
      # Histograms:
      # - raw mutation frequencies
      save_comb_hist_plot(split_dfs[[i]][[j]], 100, 0, 1.0, 
                          USE_VAF=F, log_y=F)

      # 1/f plots:
      # - raw mutation frequencies
      save_comb_1of_plot(split_dfs[[i]][[j]], 0.0001, 1.0, 10000, 
                         USE_VAF=F, COMB_SIMS=T)

    }
    
  }
  
}







