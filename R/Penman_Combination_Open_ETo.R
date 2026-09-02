# ===============================================================
# R/Penman_Combination_Open_ETo.R
#
# Open-water evaporation, Penman combination.
# Requires R/utils_physics.R to be sourced first.
#
# G = WATER heat flux (MJ m^-2 d^-1). Default 0 is a PLACEHOLDER.
#     For a deep reservoir this term is large and drives the seasonal
#     phase lag. It must be supplied for real work.
#
# VERIFIED — Gate B passed. Reproduces McMahon et al. (2013), HESS 17,
# 1331-1363, Supplement Section S19, Worked Examples 1-3 (Alice Springs
# Airport, 20 Jul 1980, G = 0):
#   - Every psychrometric/radiation intermediate (Tmean, es, ea, vpd,
#     Ra, Rso, P, gamma, Delta) matches to ~1e-5, i.e. the shared
#     prepare_forcing() chain is independently confirmed against a
#     second published source, not just FAO-56 Uccle/Lyon.
#   - Rnl/Rn differ from the paper by ~0.06%, fully explained by FAO-56
#     Eq. 39's T + 273.16 vs McMahon's T + 273.2 (both are legitimate,
#     see utils_physics.R header note on this convention).
#   - With wind_fn = "penman1956" (McMahon's own choice, Eq. S4.3),
#     Eopen_mm = 2.9808 vs published EPenOW = 2.9797 (0.04%, within the
#     rounding cascade of a hand-worked example).
#
# WIND FUNCTION: the default, "fao56_shuttleworth" (6.43*(1+0.536*u2)
# MJ m^-2 d^-1 kPa^-1, i.e. Shuttleworth's (1993) recasting of the
# land-calibrated FAO/Penman wind function), is a DIFFERENT published
# parameterisation from Penman's own (1956) open-water wind function
# (f(u) = 1.313 + 1.381*u2, used by McMahon et al. throughout their
# worked examples). For the Alice Springs case the two give Eopen_mm of
# 3.548 vs 2.981 mm/day respectively — a 19% difference from wind
# function choice alone, at typical open-water u2. This is a real
# methodological choice, not a bug; report it as a sensitivity in
# Methods (see PROJECT_CONTEXT.md section 5, "Wind conversion choice").
# ===============================================================

Eopen_Penman <- function(df, lat_deg, elev_m, albedo = ALBEDO_WATER,
                         G = 0, Tw_C = NULL, Ld_MJ_m2_day = NULL,
                         wind_fn = c("fao56_shuttleworth", "penman1956"),
                         verbose = TRUE) {

  wind_fn <- match.arg(wind_fn)
  
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
  rad_term <- (f$Delta / (f$Delta + f$gamma)) * ((Rn - G) / LAMBDA_MJ)

  if (wind_fn == "fao56_shuttleworth") {
    # Shuttleworth (1993)-form wind function, calibrated on land-based
    # wind. MJ m^-2 d^-1 kPa^-1; divided by lambda to reach mm/day.
    aero_term <- (f$gamma / (f$Delta + f$gamma)) *
      (6.43 * (1 + 0.536 * df$u2) / LAMBDA_MJ) * f$vpd
  } else {
    # Penman (1956) open-water wind function, McMahon et al. (2013)
    # Eq. S4.3. Already in mm/day/kPa - no lambda division.
    aero_term <- (f$gamma / (f$Delta + f$gamma)) *
      (1.313 + 1.381 * df$u2) * f$vpd
  }
  
  cbind(df,
        Tmean_C = f$Tmean, es = f$es, ea = f$ea, vpd = f$vpd,
        Ra_MJ = f$Ra, Rso_MJ = f$Rso, Rns_MJ = Rns, Rnl_MJ = Rnl,
        Rn_MJ = Rn, G_MJ = G,
        rad_mm = rad_term, aero_mm = aero_term,
        Eopen_mm = pmax(rad_term + aero_term, 0))
}