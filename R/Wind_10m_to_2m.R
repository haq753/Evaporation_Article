wind10_to_wind2 <- function(uz, z_m = 10, method = c("fao56", "log"), z0 = NULL) {
  
  method <- match.arg(method)
  
  if (any(uz < 0, na.rm = TRUE)) stop("Negative wind speeds in input.")
  
  ratio <- switch(method,
                  fao56 = 4.87 / log(67.8 * z_m - 5.42),
                  log   = {
                    if (is.null(z0)) stop("method = 'log' requires z0 (e.g. 0.0002 for open water).")
                    if (z0 <= 0 || z0 >= 2) stop("z0 must be > 0 and < 2 m.")
                    log(2 / z0) / log(z_m / z0)
                  }
  )
  
  uz * ratio
}