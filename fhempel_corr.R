################################################################################
# FUNCTION TO APPLY THE CALIBRATION PROCEDURE DESCRIBED IN:
#   Hempel et al. A trend-preserving bias correction - the ISI-MIP approach.
#      Earth Syst Dynam, 2013;4:219-236.
#
# ARGUMENTS:
#   - series: DATE + OPTIONALLY MULTIPLE SERIES
#   - correction: MATRIX WITH CORRECTION PARAMETERS (OR LIST OF)
#
# Author: Antonio Gasparrini - GNU General Public License (version 3)
# Update: 04 March 2019
################################################################################

fhempel_corr <- function(series, correction) {
  #
  # CHECK CONSISTENCY
  if (ncol(series) < 2) stop("series must have at least two variables (date and series")
  if (class(series[, 1]) != "Date") {
    stop("Class of first variable of 'series' must be 'Date'")
  }
  if (!is.list(correction)) correction <- list(correction)
  if (any(sapply(correction, nrow) != 12) | any(sapply(correction, ncol) != 2)) {
    stop("'correction' must be a (list of) matrix with 12 rows and 2 columns")
  }
  if (ncol(series) - 1 != length(correction)) {
    stop("number of series and correction matrices not consistent")
  }
  #
  # APPLY TO EACH MODELLED SERIES
  out <- lapply(seq(ncol(series) - 1), function(j) {
    #
    # EXTRACT MONTH AND YEAR
    day <- as.numeric(format(series[, 1], format = "%d"))
    month <- as.numeric(format(series[, 1], format = "%m"))
    if (length(unique(month)) != 12) {
      stop("some months are not reprensented in 'series'")
    }
    year <- as.numeric(format(series[, 1], format = "%Y"))
    monthyear <- factor(paste(year, month, sep = "-"))
    #
    # COMPUTE AVERAGES AND RESIDUALS OF THE SERIES
    mavg <- tapply(series[, j + 1], month, mean, na.rm = T)
    myavg <- tapply(series[, j + 1], monthyear, mean, na.rm = T)
    res <- series[, j + 1] - myavg[monthyear]
    #
    # DERIVE THE QUANTITIES TO REMOVE DISCONTINUITIES FROM CORRECTION
    # (SEE Hempel et al, PAGE 228)
    nday <- tapply(series[, j + 1], monthyear, length)
    d <- (day - 1) / (nday[monthyear] - 1) - 0.5
    dm <- 0.5 * (abs(d) - d)
    d0 <- 1 - abs(d)
    dp <- 0.5 * (abs(d) + d)
    #
    # DERIVE THE CORRECTIONS ACCOUNTING FOR DISCONTINUITIES
    C <- corr[[j]][, 1]
    B <- corr[[j]][, 2]
    Cm <- C[c(12, 1:11)]
    Bm <- B[c(12, 1:11)]
    Cp <- C[c(2:12, 1)]
    Bp <- B[c(2:12, 1)]
    C <- dm * Cm[month] + d0 * C[month] + dp * Cp[month]
    B <- dm * Bm[month] + d0 * B[month] + dp * Bp[month]
    #
    # RETURN THE CORRECTED SERIES
    return(myavg[monthyear] + C + res * B)
  })
  #
  # ADD DATE AND RENAME
  out <- cbind(series[1], as.data.frame(do.call(cbind, out)))
  dimnames(out) <- dimnames(series)
  #
  # RETURN
  return(out)
}
