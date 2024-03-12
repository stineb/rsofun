library(dplyr)
library(tidyr)
library(rsofun)
library(ggplot2)
library(patchwork)
library(readr)
library(lubridate)

## Parameters ------------------------
# in LT review paper v3
pars <- list(

  # P-model
  kphio                 = 0.04998,    # setup ORG in Stocker et al. 2020 GMD
  kphio_par_a           = 0.0,        # set to zero to disable temperature-dependence of kphio
  kphio_par_b           = 1.0,
  soilm_thetastar       = 0.6 * 240,  # to recover old setup with soil moisture stress
  soilm_betao           = 0.0,
  beta_unitcostratio    = 146.0,
  rd_to_vcmax           = 0.014,      # value from Atkin et al. 2015 for C3 herbaceous
  tau_acclim            = 30.0,
  kc_jmax               = 0.41,

  # Plant
  f_nretain             = 0.500000,
  fpc_tree_max          = 0.950000,
  growtheff             = 0.600000,
  r_root                = 2*0.913000,
  r_sapw                = 2*0.044000,
  exurate               = 0.003000,

  k_decay_leaf          = 1.90000,
  k_decay_root          = 1.90000,
  k_decay_labl          = 1.90000,
  k_decay_sapw          = 1.90000,

  r_cton_root           = 37.0000,
  r_cton_wood           = 100.000,
  r_cton_seed           = 15.0000,
  nv_vcmax25            = 0.02 * 13681.77, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.695; previously: 5000.0,
  ncw_min               = 0.08 * 1.116222, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.691; previously used: 0.056,
  r_n_cw_v              = 0, # assumed that LMA is independent of Vcmax25; previously: 0.1,
  r_ctostructn_leaf     = 1.3 * 45.84125, # see ln_cn_review/vignettes/analysis_leafn_vcmax_field.Rmd, l.699; previously used: 80.0000,
  kbeer                 = 0.500000,

  # Phenology (should be PFT-specific)
  gddbase               = 5.0,
  ramp                  = 0.0,
  phentype              = 2.0,

  # Soil physics (should be derived from params_soil, fsand, fclay, forg, fgravel)
  perc_k1               = 5.0,
  thdiff_wp             = 0.2,
  thdiff_whc15          = 0.8,
  thdiff_fc             = 0.4,
  forg                  = 0.01,
  wbwp                  = 0.029,
  por                   = 0.421,
  fsand                 = 0.82,
  fclay                 = 0.06,
  fsilt                 = 0.12,

  # Water and energy balance
  kA                    = 107,
  kalb_sw               = 0.17,
  kalb_vis              = 0.03,
  kb                    = 0.20,
  kc                    = 0.25,
  kCw                   = 1.05,
  kd                    = 0.50,
  ke                    = 0.0167,
  keps                  = 23.44,
  kWm                   = 220.0,
  kw                    = 0.26,
  komega                = 283.0,
  maxmeltrate           = 3.0,

  # Soil BGC
  klitt_af10            = 1.2,
  klitt_as10            = 0.35,
  klitt_bg10            = 0.35,
  kexu10                = 50.0,
  ksoil_fs10            = 0.021,
  ksoil_sl10            = 7.0e-04,
  ntoc_crit1            = 0.45,
  ntoc_crit2            = 0.76,
  cton_microb           = 10.0,
  cton_soil             = 9.77,
  fastfrac              = 0.985,

  # N uptake
  eff_nup               = 0.0001000,
  minimumcostfix        = 1.000000,
  fixoptimum            = 25.15000,
  a_param_fix           = -3.62000,
  b_param_fix           = 0.270000,

  # Inorganic N transformations (re-interpreted for simple ntransform model)
  maxnitr               =  0.00005,

  # Inorganic N transformations for full ntransform model (not used in simple model)
  non                   = 0.01,
  n2on                  = 0.0005,
  kn                    = 83.0,
  kdoc                  = 17.0,
  docmax                = 1.0,
  dnitr2n2o             = 0.01,

  # Additional parameters - previously forgotten
  frac_leaf             = 0.5,           # after wood allocation
  frac_wood             = 0,           # highest priority in allocation
  frac_avl_labl         = 0.1,

  # for development
  tmppar                = 9999,

  # simple N uptake module parameters
  nuptake_kc            = 250,
  nuptake_kv            = 5,
  nuptake_vmax          = 0.15

)

## CH-Oe1 climate forcing -------------------
filnam <- "data-raw/df_drivers_ch_oe1.rds"
overwrite <- FALSE

if (!file.exists(filnam) || overwrite){
  
  df_drivers_ch_oe1 <- read_rds("~/data/FluxDataKit/rsofun_driver_data_clean.rds") |> 
    filter(sitename == "CH-Oe1")
  
  df_ndep <- ingestr::ingest(
    df_drivers_ch_oe1 |>
      unnest(site_info) |>
      select(sitename, lon, lat) |>
      mutate(year_start = min(lubridate::year(df_drivers_ch_oe1$forcing[[1]]$date)),
             year_end = max(lubridate::year(df_drivers_ch_oe1$forcing[[1]]$date))),
    source    = "ndep",
    timescale = "y",
    dir       = "~/data/ndep_lamarque/",
    verbose   = FALSE
  )
  
  df_ndep_mean <- df_ndep |>
    unnest(data) |>
    summarise(noy = mean(noy), nhx = mean(nhx))
  
  ## add new required columns to forcing
  use_cseed <- 5
  cn_seed <- 20
  use_nseed <- use_cseed / cn_seed
  
  df_drivers_ch_oe1 <- df_drivers_ch_oe1 |>
    mutate(forcing = purrr::map(forcing, ~mutate(.,
                                                 fharv = 0.0,
                                                 dno3 = df_ndep_mean$noy / 365,
                                                 dnh4 = df_ndep_mean$nhx / 365)),
           forcing = purrr::map(forcing, ~mutate(.,
                                                 fharv = ifelse(month(date) == 7 & mday(date) == 15, 0.0, 0.0),
                                                 cseed = ifelse(month(date) == 2 & mday(date) == 15, use_cseed, 0.0),
                                                 nseed = ifelse(month(date) == 2 & mday(date) == 15, use_nseed, 0.0))))
  
  df_drivers_ch_oe1$params_siml[[1]]$spinupyears <- 2002
  df_drivers_ch_oe1$params_siml[[1]]$recycle <- 5
  
  write_rds(df_drivers_ch_oe1, file = filnam)
  
} else {
  
  df_drivers_ch_oe1 <- readRDS(filnam)
  
}

df_drivers <- df_drivers_ch_oe1

## Define whether to use interactive C-N cycling
df_drivers$params_siml[[1]]$c_only <- FALSE

### Synthetic forcing: Mean seasonal cycle -----------------------
df_drivers$forcing[[1]] <- df_drivers$forcing[[1]] |>
  filter(!(lubridate::month(date) == 2 & lubridate::mday(date) == 29))
df_meanann <- df_drivers$forcing[[1]] |>
  mutate(doy = lubridate::yday(date)) |>
  group_by(doy) |>
  summarise(across(where(is.double), .fns = mean)) |>
  filter(!(doy == 366))
nyears <- df_drivers$forcing[[1]] |>
  mutate(year = lubridate::year(date)) |>
  pull(year) |>
  unique() |>
  length()
df_drivers$forcing[[1]] <- purrr::map_dfr(
  as.list(seq(nyears)),
  ~{df_meanann}) |>
  mutate(date = df_drivers$forcing[[1]]$date)

### 100 years at elevated CO2 ---------------
# from 2009 onwards
df_meanann_ele <- df_meanann |>
  mutate(co2 = co2*2)

nyears_ele <- 100

dates <- tibble(
  date = seq(ymd("2009-01-01"), ymd("2108-12-31"), by = "days")) |>
  mutate(month = month(date), dom = mday(date)) |>
  filter(!(month == 2 & dom == 29))

df_drivers$forcing[[1]] <- bind_rows(
  df_drivers$forcing[[1]],
  purrr::map_dfr(
    as.list(seq(nyears_ele)),
    ~{df_meanann_ele}) |>
    mutate(date = dates$date)
    )

## Model run ------------------------
output <- runread_cnmodel_f(
  df_drivers,
  par = pars
)

adf <- output$data[[1]] |> 
  mutate(year = year(date)) |> 
  group_by(year) |> 
  summarise(
    lai = max(lai),
    gpp = sum(gpp),
    npp = sum(npp),
    npp_leaf = sum(npp_leaf),
    npp_root = sum(npp_root),
    npp_wood = sum(npp_wood),
    nup = sum(nup),
    netmin = sum(netmin),
    ninorg = mean(pno3 + pnh4),
    nloss = sum(nloss)
  )

# keep first 7 years for visualising seasonal course
ddf <- output$data[[1]][1:(7*365),]

## Visualisations  ------------------------
### Daily time series ---------------------------
# for visualising seasonality
#### LAI -----------------------------
gg1 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, lai)) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("LAI (m"^2, " m"^-2, ")")))

#### GPP -----------------------------
gg2 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, gpp)) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("GPP (gC m"^-2, " d"^-1, ")")))

#### NEE -----------------------------
gg3 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, gpp - rleaf - rwood - rroot - rhet)) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("NEE (gC m"^-2, " d"^-1, ")")))

#### cumulative NEE -----------------------------
gg4 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, cumsum(gpp - rleaf - rwood - rroot - rhet))) + 
  geom_line() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Date", y = expression(paste("Cumulative NEE by GPP (gC m"^-2, " d"^-1, ")")))

ddf |> 
  as_tibble() |> 
  ggplot(aes(date, cumsum(npp - rhet))) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("Cumulative NEE by NPP (gC m"^-2, " d"^-1, ")")))

gg1 / gg2 / gg3 / gg4

#### all pools --------------
ggtest <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, cleaf + croot + cwood + cseed + cresv + clabl + clitt + csoil)) + 
  geom_line() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Date", y = expression(paste("Total C (gC m"^-2, ")")))

#### NPP -----------------------------
gg5 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, npp)) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("NPP (gC m"^-2, " d"^-1, ")")))

#### NPP ----------------------------- fraction to leaves
gg6 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, npp_leaf/npp)) + 
  geom_line() +
  labs(x = "Date", y = "Fraction of leaf BP")

#### NPP ----------------------------- fraction to root
gg7 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, npp_root/npp)) + 
  geom_line() +
  labs(x = "Date", y = "Fraction of root BP")

#### NPP ----------------------------- fraction to wood
gg8 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, npp_wood/npp)) + 
  geom_line() +
  labs(x = "Date", y = "Fraction of wood BP")

#### BPE -----------------------------
gg9 <- ddf |> 
  mutate(year = lubridate::year(date)) |> 
  group_by(year) |> 
  summarise(npp = sum(npp),
            gpp = sum(gpp)) |> 
  ggplot(aes(year, npp/gpp)) + 
  geom_line() +
  labs(x = "Year", x = "BPE (unitless)")

#### Cleaf -----------------------------
gg10 <- ddf |> 
  ggplot(aes(date, cleaf)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Leaf C (gC m"^-2, ")")))

#### Croot -----------------------------
gg11 <- ddf |> 
  ggplot(aes(date, croot)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Root C (gC m"^-2, ")")))

#### Cwood -----------------------------
gg26 <- ddf |> 
  ggplot(aes(date, cwood)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Wood C (gC m"^-2, ")")))

#### Clabl -----------------------------
gg12 <- ddf |> 
  ggplot(aes(date, clabl)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Labile C (gC m"^-2, ")")))

#### Cresv -----------------------------
gg13 <- ddf |> 
  ggplot(aes(date, cresv)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Reserves C (gC m"^-2, ")")))

#### RMF -----------------------------
gg14 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, croot/(croot + cleaf + cwood))) + 
  geom_line() +
  labs(x = "Date", y = "Root mass fraction")

#### Clitt -----------------------------
gg15 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, clitt)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Litter C (gC m"^-2, ")")))

#### Csoil -----------------------------
gg16 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, csoil)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Soil C (gC m"^-2, ")")))

#### CNleaf -----------------------------
gg17 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, cleaf/nleaf)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Leaf C:N (gC gN"^-1, ")")))

#### CNlitt -----------------------------
gg18 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, clitt/nlitt)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Litter C:N (gC gN"^-1, ")")))

#### CNsoil -----------------------------
gg19 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, csoil/nsoil)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Soil C:N (gC gN"^-1, ")")))

#### Ninorg -----------------------------
gg20 <- ddf |> 
  ggplot(aes(date, ninorg)) + 
  geom_line() +
  labs(x = "Year", x = expression(paste("Soil inorganic N (gN m"^-2, ")")))

#### Netmin -----------------------------
gg21 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, netmin)) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("Net N mineralization (gN m"^-2, " d"^-1, ")")))

#### Nup -----------------------------
gg22 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, nup)) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("N uptake (gN m"^-2, " d"^-1, ")")))

#### Nloss ----------------------------- XXX problem: cannot be so high
gg23 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, nloss)) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("N loss (gN m"^-2, " d"^-1, ")")))

#### Nup -----------------------------
gg24 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, nfix)) + 
  geom_line() +
  labs(x = "Date", y = expression(paste("N fixation (gN m"^-2, " d"^-1, ")")))

#### tsoil -----------------------------
gg25 <- ddf |> 
  as_tibble() |> 
  ggplot(aes(date, tsoil)) + 
  geom_line() +
  labs(x = "Date", y = "°C")

#### All -----------------
ggout <- cowplot::plot_grid(
  gg1, 
  gg2, 
  gg3, 
  gg4, 
  gg5, 
  gg6, 
  gg7, 
  gg8, 
  gg9, 
  gg10, 
  gg11, 
  gg26,
  gg12, 
  gg13, 
  gg14, 
  gg15, 
  gg16, 
  gg17, 
  gg18, 
  gg19, 
  gg20, 
  gg21, 
  gg22, 
  gg23,
  gg24,
  ncol = 1
  )

ggsave(here::here("fig/tseries.pdf"), 
       plot = ggout,
       width = 8,
       height = 35 )

### Annual time series: -----------------------
# for visualising response to step increase
#### LAI -----------------------------
gg1 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, lai)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = expression(paste("LAI (m"^2, " m"^-2, ")")))

#### GPP -----------------------------
gg2 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, gpp)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = expression(paste("GPP (gC m"^-2, " d"^-1, ")")))

#### NPP -----------------------------
gg5 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, npp)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = expression(paste("NPP (gC m"^-2, " d"^-1, ")")))

#### NPP fraction to leaves ----------------------------- 
gg6 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, npp_leaf/npp)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = "Fraction of leaf BP")

#### NPP fraction to roots ----------------------------- 
gg7 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, npp_root/npp)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = "Fraction of root BP")

#### NPP fraction to wood ----------------------------- 
gg8 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, npp_wood/npp)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = "Fraction of wood BP")

#### N uptake ----------------------------- 
gg9 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, nup)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = "N uptake")

#### N inorganic pool ----------------------------- 
gg10 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, ninorg)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = "Soil mineral N")

#### N loss ----------------------------- 
gg11 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, nloss)) + 
  geom_vline(xintercept = 2009, linetype = "dotted") +
  geom_line() +
  labs(x = "Date", y = "N loss")

#### N mineralisation ----------------------------- 
gg12 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, netmin)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = "Net N mineralisation")

#### N cycle openness ----------------------------- 
gg13 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, nloss/netmin)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = "N cycle openness")

#### N cost ----------------------------- 
gg14 <- adf |> 
  as_tibble() |> 
  ggplot(aes(year, npp_root/nup)) + 
  geom_line() +
  geom_vline(xintercept = 2009, linetype = "dotted") +
  labs(x = "Date", y = "N cost")

#### All -----------------
ggout <- cowplot::plot_grid(
  gg1, 
  gg2, 
  gg5, 
  gg6, 
  gg7, 
  gg8, 
  gg9,
  gg10,
  gg11,
  gg12,
  gg13,
  gg14,
  ncol = 1
)

ggsave(here::here("fig/tseries_eco2.pdf"), 
       plot = ggout,
       width = 8,
       height = 16 )

### Annual total N budget -----------------
meanadf <- ddf |> 
  as_tibble() |> 
  left_join(
    df_drivers$forcing[[1]],
    by = "date"
  ) |> 
  mutate(year = lubridate::year(date)) |> 
  group_by(year) |> 
  summarise(across(where(is.numeric), sum)) |> 
  summarise(across(where(is.numeric), mean))
  
meanadf |> 
  mutate(ndep = dno3 + dnh4) |> 
  select(ndep, nloss, netmin, nup) |> 
  pivot_longer(
    cols = c(ndep, nloss, netmin, nup),
    names_to = "var",
    values_to = "val"
  ) |> 
  ggplot(aes(var, val)) +
  geom_bar(stat = "identity")

### Spinup-----------------
## read (experimental) files
aout <- read_fwf(file = "out/out_rsofun.a.csoil.txt", col_types = "in") |>
  setNames(c("year", "csoil")) |>
  left_join(
    read_fwf(file = "out/out_rsofun.a.nsoil.txt", col_types = "in") |>
      setNames(c("year", "nsoil")),
    by = "year"
  )

# soil C spinup
aout |>
  # slice(1000:2008) |> 
  ggplot(aes(year, csoil)) +
  
  # first soil equilibration year
  geom_vline(xintercept = 600, linetype = "dotted") +

  # start free allocation
  geom_vline(xintercept = 900, linetype = "dotted") +
  
  # second soil equilibration year
  geom_vline(xintercept = 1500, linetype = "dotted") +
  
  geom_line() + 
  theme_classic()

# soil N spinup
aout |>
  # slice(1000:2008) |> 
  ggplot(aes(year, nsoil)) +
  
  # first soil equilibration year
  geom_vline(xintercept = 600, linetype = "dotted") +

  # start free allocation
  geom_vline(xintercept = 900, linetype = "dotted") +
  
  # second soil equilibration year
  geom_vline(xintercept = 1500, linetype = "dotted") +
  
  geom_line() + 
  theme_classic()

### Response ratios ---------------------------
df_out <- output$data[[1]] |> 
  mutate(leaf_cn = cleaf/nleaf, 
         root_shoot = croot/cleaf, 
         n_inorg = pno3 + pnh4,
         anpp = npp_leaf + npp_wood, 
         bnpp = npp_root + cex) |> 
  select(date, asat, gpp, vcmax, jmax, gs = gs_accl, narea, leaf_cn, lai, cleaf, 
         croot, root_shoot, nup, n_inorg, anpp, bnpp)

df_amb <- df_out |> 
  filter(year(date) < 2010) |> 
  summarise(across(where(is.numeric), mean))

df_ele <- df_out |> 
  filter(year(date) %in% 2010:2015) |> 
  summarise(across(where(is.numeric), mean))

df_ele2 <- df_out |> 
  filter(year(date) %in% 2100:2110) |> 
  summarise(across(where(is.numeric), mean))

df_exp <- bind_rows(df_amb, df_ele)
df_rr  <- log(df_exp[2,]/df_exp[1,]) |> 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "response") |> 
  mutate(variable = factor(variable, 
                           levels = rev(c("gpp", "asat", "vcmax", "jmax", "gs", 
                                          "narea", 
                                          "leaf_cn", "lai", "cleaf", "anpp",
                                          "croot", "bnpp", "root_shoot", "nup", 
                                          "n_inorg"))))

df_exp2 <- bind_rows(df_amb, df_ele2)
df_rr2  <- log(df_exp2[2,]/df_exp2[1,]) |> 
  pivot_longer(cols = everything(), names_to = "variable", values_to = "response") |> 
  mutate(variable = factor(variable, 
                           levels = rev(c("gpp", "asat", "vcmax", "jmax", "gs", 
                                          "narea", 
                                          "leaf_cn", "lai", "cleaf", "anpp",
                                          "croot", "bnpp", "root_shoot", "nup", 
                                          "n_inorg"))))

ggrr <- ggplot() +
  geom_point(aes(variable, response, color = "Long-term"), data = df_rr2, size = 2) +
  geom_point(aes(variable, response, color = "First 5 years"), data = df_rr, size = 2) +
  geom_hline( yintercept = 0.0, linewidth = 0.5, linetype = "dotted" ) +
  labs(x = "Variable", y = "Log Response Ratio") +
  coord_flip() +
  labs(title = "cnmodel prediction", x = "Response to eCO2")

### Function visualisations -------------
ggplot() +
  geom_function(fun = calc_f_seed) +
  xlim(-0.02,0.02) +
  # geom_vline(xintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")

calc_ft_growth <- function(temp){
  yy <- 1 / (1 + exp(-(temp-5)))
  return(yy)
}

ggplot() +
  geom_function(fun = calc_ft_growth) +
  xlim(-10, 20) +
  geom_vline(xintercept = 0, linetype = "dotted")


calc_f_nup <- function(conc){
  jmax <- 2
  k_nup <- 1
  yy <- jmax / (1 + jmax/(k_nup * conc))
  return(yy)
}

ggplot() +
  geom_function(fun = calc_f_nup) +
  xlim(0, 100) +
  geom_vline(xintercept = 0, linetype = "dotted")


calc_ftemp <- function(temp){
  E0 = 308.56
  T0 = 227.13
  Tzero = 273.15
  ref_temp_local = 22
  ftemp = exp(E0 * ((1.0 / (ref_temp_local + Tzero - T0)) - (1.0 / (temp + Tzero - T0))))
  return(ftemp)
}

ggplot() +
  geom_function(fun = calc_ftemp) +
  xlim(0, 40) +
  geom_vline(xintercept = 0, linetype = "dotted")

calc_decay <- function(t){
  k = 1/100
  y = exp(-k * t)
  return(y)
}

ggplot() +
  geom_function(fun = calc_decay) +
  geom_point(aes(x = 10, y = 0.9), col = "red" ) +
  geom_point(aes(x = 20, y = 0.8), col = "red" ) +
  geom_point(aes(x = 50, y = 0.5), col = "red" ) +
  xlim(0, 500) +
  ylim(0, 1)

calc_dc <- function(t){
  k = 1/100
  dc = 1 - exp(-k * t)
  return(dc)
}

ggplot() +
  geom_function(fun = calc_dc) +
  geom_point(aes(x = 10, y = 0.1), col = "red" ) +
  geom_point(aes(x = 20, y = 0.2), col = "red" ) +
  geom_point(aes(x = 50, y = 0.5), col = "red" ) +
  xlim(0, 500) +
  ylim(0, 1)


## Write output to file --------------------
readr::write_csv(as_tibble(df_exp), file = paste0(here::here(), "/data/df_exp_co2.csv"))
readr::write_csv(as_tibble(output), file = paste0(here::here(), "/data/output_cnmodel_co2.csv"))
# readr::write_csv(as_tibble(output), file = "../data/output_cnmodel_nfert.csv")
