# ===============================================================
# R/Penman_Combination_ETo.R
#
# Open-water evaporation, Penman combination.
# Requires R/utils_physics.R to be sourced first.
#
# G = WATER heat flux (MJ m^-2 d^-1). Default 0 is a PLACEHOLDER.
#     For a deep reservoir this term is large and drives the seasonal
#     phase lag. It must be supplied for real work.
#
# NOT YET VERIFIED against a published worked example. Gate B.
# ===============================================================

Eopen_Penman <- function(df, lat_deg, elev_m, albedo = ALBEDO_WATER,
                         G = 0, Tw_C = NULL, Ld_MJ_m2_day = NULL,
                         verbose = TRUE) {
  
  if (!exists("prepare_forcing"))
    stop("Source R/utils_physics.R before this file.")
  
  miss <- setdiff(REQUIRED_FORCING, names(df))
  if (length(miss))
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  if (is.null(Tw_C) != is.null(Ld_MJ_m2_day))
    stop("Tw_C and Ld_MJ_m2_day must be supplied together.")
  
  if (identical(G, 0) && verbose)
    message("Eopen_Penman: G = 0 (no heat storage). ",
            "Placeholder only for a deep reservoir.")
  
  # --- sorting, with externally supplied vectors carried along ------
  # Previously df was reordered while G and Tw_C were not, silently
  # misaligning the heat storage term against the forcing whenever the
  # input arrived unsorted.
  ord <- order(df$date)
  if (is.numeric(G) && length(G) == nrow(df))       G    <- G[ord]
  if (!is.null(Tw_C) && length(Tw_C) == nrow(df))   Tw_C <- Tw_C[ord]
  if (!is.null(Tw_C) && !length(Tw_C) %in% c(1L, nrow(df)))
    stop("Tw_C must be length 1 or nrow(df). Got ", length(Tw_C), ".")
  if (is.numeric(G) && !length(G) %in% c(1L, nrow(df)))
    stop("G must be length 1 or nrow(df). Got ", length(G), ".")
  df <- df[ord, ]
  
  f <- prepare_forcing(df, lat_deg, elev_m, verbose = verbose)
  
  Rns <- (1 - albedo) * f$Rs
  
  # --- net longwave: branch is reported, never silent ---------------
  if (is.null(Tw_C)) {
    if (verbose)
      message("Eopen_Penman: longwave branch = FAO-56 Eq. 39 ",
              "(air temperature surrogate). Approximate for open water.")
    Rnl <- Rnl_fao56(df$Tmax_C, df$Tmin_C, f$ea, f$Rs, f$Rso)
  } else {
    if (!Ld_MJ_m2_day %in% names(df))
      stop("Column '", Ld_MJ_m2_day, "' not found in df.")
    if (verbose)
      message("Eopen_Penman: longwave branch = open water ",
              "(eps*sigma*Tw^4 - Ld), Ld from '", Ld_MJ_m2_day, "'.")
    Rnl <- Rnl_water(Tw_C, df[[Ld_MJ_m2_day]])
  }
  
  Rn <- Rns - Rnl
  
  # --- Penman combination, open water -------------------------------
  # lambda stays 2.45: the 6.43 wind constant is expressed in
  # MJ m^-2 d^-1 kPa^-1 and was derived on that basis.
  rad_term  <- (f$Delta / (f$Delta + f$gamma)) * ((Rn - G) / LAMBDA_MJ)
  wind_term <- 6.43 * (1 + 0.536 * df$u2)
  aero_term <- (f$gamma / (f$Delta + f$gamma)) * (wind_term / LAMBDA_MJ) * f$vpd
  
  cbind(df,
        Tmean_C = f$Tmean, es = f$es, ea = f$ea, vpd = f$vpd,
        Ra_MJ = f$Ra, Rso_MJ = f$Rso, Rns_MJ = Rns, Rnl_MJ = Rnl,
        Rn_MJ = Rn, G_MJ = G,
        rad_mm = rad_term, aero_mm = aero_term,
        Eopen_mm = pmax(rad_term + aero_term, 0))
}