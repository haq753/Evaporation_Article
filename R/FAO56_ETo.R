# ===============================================================
# R/FAO56_ETo.R
#
# FAO-56 reference evapotranspiration (grass reference surface).
# Requires R/utils_physics.R to be sourced first.
#
# Verified against FAO-56 Example 18 (3.8796 mm/d) and
# Example 20 (4.5604 mm/d).
# ===============================================================

ETo_FAO56 <- function(df, lat_deg, elev_m,
                      albedo = ALBEDO_GRASS, G = 0, verbose = TRUE) {
  
  if (!exists("prepare_forcing"))
    stop("Source R/utils_physics.R before this file.")
  
  if (albedo != ALBEDO_GRASS && verbose)
    warning("FAO-56 ETo is defined for albedo = 0.23. You passed ", albedo,
            ". For open water use Eopen_Penman() instead.")
  
  f <- prepare_forcing(df, lat_deg, elev_m, verbose = verbose)
  
  Rns <- (1 - albedo) * f$Rs
  Rnl <- Rnl_fao56(df$Tmax_C, df$Tmin_C, f$ea, f$Rs, f$Rso)
  Rn  <- Rns - Rnl
  
  # FAO-56 Eq. 6. Note: T + 273 here, not 273.16. As published.
  ETo <- (0.408 * f$Delta * (Rn - G) +
            f$gamma * (900 / (f$Tmean + 273)) * df$u2 * f$vpd) /
    (f$Delta + f$gamma * (1 + 0.34 * df$u2))
  
  cbind(df,
        Tmean_C = f$Tmean, es = f$es, ea = f$ea, vpd = f$vpd,
        Ra_MJ = f$Ra, Rso_MJ = f$Rso, Rns_MJ = Rns, Rnl_MJ = Rnl,
        Rn_MJ = Rn, ETo_mm = pmax(ETo, 0))
}