################################################################################
# FUNCTION TO APPLY THE CALIBRATION PROCEDURE DESCRIBED IN:
#   Hempel et al. A trend-preserving bias correction - the ISI-MIP approach.
#      Earth Syst Dynam, 2013;4:219-236.
#
# ARGUMENTS:
#   - obs: DATE + SINGLE OBSERVED HISTORICAL SERIES
#   - mod: DATE + OPTIONALLY MULTIPLE MODELLED SERIES
#   - add: LOGICAL. IF TRUE, THE ADDITIVE CORRECTION IS APPLIED
#   - mult: LOGICAL. IF TRUE, THE MULTIPLICATIVE CORRECTION IS APPLIED
#   - output: RETURN THE SERIES ("series") OR THE PARAMETERS ("correction")
#
# Author: Antonio Gasparrini - GNU General Public License (version 3)
# Update: 04 March 2019
################################################################################

fhempel_main <- function(obs, mod, add = TRUE, mult = TRUE, output = "series") {
  #
  # CHECK THE OUTPUT
  output <- match.arg(output, c("series", "correction"))
  #
  # CHECK CONSISTENCY
  if (ncol(obs) != 2) stop("obs must have two variables (date and series")
  if (ncol(mod) < 2) stop("mod must have at least two variables (date and series")
  if (!base:::any(class(obs[, 1]) %in% "Date") | !base:::any(class(mod[, 1]) %in% "Date")) {
    stop("Class of first variable of 'obs' and 'mod' must be 'Date'")
  }
  #
  # APPLY TO EACH MODELLED SERIES
  out <- lapply(seq(ncol(mod) - 1), function(j) {
    #
    # IDENTIFY PERIOD WITH NO MISSING
    ind <- intersect(obs[, 1], mod[, 1])
    notna <- complete.cases(obs[obs[, 1] %in% ind, 2], mod[mod[, 1] %in% ind, j + 1])
    indobs <- seq(nrow(obs))[obs[, 1] %in% ind][notna]
    indmod <- seq(nrow(mod))[mod[, 1] %in% ind][notna]
    if (length(indobs) == 0) stop("no common non-missing days in 'obs' and 'mod'")
    #
    # EXTRACT MONTH AND YEAR
    month <- as.numeric(format(obs[indobs, 1], format = "%m"))
    if (length(unique(month)) != 12) {
      stop("some months are not reprensented in the overlapping period")
    }
    year <- as.numeric(format(obs[indobs, 1], format = "%Y"))
    monthyear <- factor(paste(year, month, sep = "-"))
    #
    # COMPUTE ADDITIVE CORRECTION
    mavgobs <- tapply(obs[indobs, 2], month, mean, na.rm = T)
    mavgmod <- tapply(mod[indmod, j + 1], month, mean, na.rm = T)
    C <- mavgobs - mavgmod
    if (!add) C[] <- 0
    #
    # RESIDUALS FROM MONTHLY/YEAR AVERAGES, THEN MULTIPLICATIVE CORRECTION
    myavgobs <- tapply(obs[indobs, 2], monthyear, mean, na.rm = T)
    myavgmod <- tapply(mod[indmod, j + 1], monthyear, mean, na.rm = T)
    resobs <- obs[indobs, 2] - myavgobs[monthyear]
    resmod <- mod[indmod, j + 1] - myavgmod[monthyear]
    B <- sapply(
      seq(12),
      function(m) coef(lm(sort(resobs[month == m]) ~ 0 + sort(resmod[month == m])))
    )
    if (!mult) B[] <- 1
    #
    # RETURN CORRECTION IF STATED
    if (output == "correction") {
      return(matrix(c(C, B),
        ncol = 2,
        dimnames = list(unique(months(mod[, 1], abbr = TRUE)), c("add", "mult"))
      ))
    }
    #
    # OBTAIN DAY, MONTH AND YEAR FROM DATE (ORIGINAL SERIES)
    day <- as.numeric(format(mod[, 1], format = "%d"))
    month <- as.numeric(format(mod[, 1], format = "%m"))
    year <- as.numeric(format(mod[, 1], format = "%Y"))
    monthyear <- factor(paste(year, month, sep = "-"))
    #
    # DERIVE THE QUANTITIES TO REMOVE DISCONTINUITIES FROM CORRECTION
    # (SEE Hempel et al, PAGE 228)
    nday <- tapply(mod[, j + 1], monthyear, length)
    d <- (day - 1) / (nday[monthyear] - 1) - 0.5
    dm <- 0.5 * (abs(d) - d)
    d0 <- 1 - abs(d)
    dp <- 0.5 * (abs(d) + d)
    #
    # DERIVE THE CORRECTIONS ACCOUNTING FOR DISCONTINUITIES
    Cm <- C[c(12, 1:11)]
    Bm <- B[c(12, 1:11)]
    Cp <- C[c(2:12, 1)]
    Bp <- B[c(2:12, 1)]
    C <- dm * Cm[month] + d0 * C[month] + dp * Cp[month]
    B <- dm * Bm[month] + d0 * B[month] + dp * Bp[month]
    #
    # DERIVE THE CORRECTED SERIES
    myavgout <- tapply(mod[, j + 1], monthyear, mean, na.rm = T)
    resout <- mod[, j + 1] - myavgout[monthyear]
    #
    # RETURN CORRECTED SERIES IF STATED
    return(myavgout[monthyear] + C + resout * B)
  })
  #
  # RETURN CORRECTION IF STATED
  if (output == "correction") {
    names(out) <- names(mod)[-1]
    return(out)
  }
  #
  # ADD DATE AND RENAME
  out <- cbind(mod[1], as.data.frame(do.call(cbind, out)))
  dimnames(out) <- dimnames(mod)
  #
  # RETURN
  return(out)
}
