
library(rsofun)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(devtools)

load("data/biomee_gs_leuning_drivers.rda")
biomee_gs_leuning_drivers$params_siml
head(biomee_gs_leuning_drivers$forcing)

load("data/biomee_p_model_drivers.rda")
biomee_p_model_drivers$params_siml
head(biomee_p_model_drivers$forcing)

set.seed(2023)

# run the model
# build()
# install()
out <- runread_biomee_f(
  biomee_p_model_drivers,
  makecheck = TRUE,
  parallel = FALSE
)

# split out the annual data
biomee_pmodel_output <- out$data[[1]]$output_annual_tile

# plot
cowplot::plot_grid(
  biomee_pmodel_output |>
    ggplot() +
    geom_line(aes(x = year, y = GPP)) +
    theme_classic()+labs(x = "Year", y = "GPP"),
  biomee_pmodel_output |>
    ggplot() +
    geom_line(aes(x = year, y = plantC)) +
    theme_classic()+labs(x = "Year", y = "plantC")
)

# function to run biomee multiple times and detect the segmentation fault error
biomeextimes <- function(drivers,n){
  results <-  list()
  for(i in 1:n) {
    set.seed(2023)
    out <- runread_biomee_f(
    drivers,
    makecheck = TRUE,
    parallel = FALSE
  )$data[[1]]$output_annual_tile['plantC'] %>% 
    slice(tail(row_number(), 1)) 
  if(is.na(out)|out==0) message("Error: Simulation failed")
  results <- append(results, out)
  }
  save(results, file= "data/results.rda")
  return(results)
}

# run the model n times 
n = 500
biomeextimes(biomee_p_model_drivers,n) 
load("data/results.rda")
results
