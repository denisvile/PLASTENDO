# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
# Function for sigmoidal fitting ----
# €€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€€
fc.sigmoidalFitting <- function(data, plot=TRUE) { 
  idPot <- unique(data$idPot)
  # Test of the sigmoidal fitting
  fit <- NULL
  for(Btest in 1:10) for (X0test in 2:50) {
    try(fit <- nls(Area.total.mm2  ~ A / ( 1 + exp (-((day - X0) / B))), 
                   data = data, 
                   start = list(A = max(data$Area.total.mm2, na.rm = TRUE), B = Btest , X0 = X0test)),
        silent = T)
    if(!is.null(fit)) break
  }
  # Model parameters
  try(A <- summary(fit)$coef[1], silent = T)
  try(B <- summary(fit)$coef[2], silent = T)
  try(X0 <- summary(fit)$coef[3], silent = T)
  
  # Computation of the maximum expansion rate : ER = A / (4*B)
  try(ER <- 0.25 * A / B, silent = T)
  
  # Computation of the expansion duration : Duration = X0 - B * ln ((1/0.95)-1)
  try(Duration <- X0 - B * log ((1/0.95)-1), silent = T)
  
  if(is.null(fit)) {
    resultTemp <- data.frame(idPot = idPot, A = NA, B = NA, X0 = NA, ER = NA, Duration = NA)
  }
  
  if(!is.null(fit)) {
    resultTemp <- data.frame(idPot = idPot, A = A, B = B, X0 = X0, ER = ER, Duration = Duration)
  }
  
  fun.sigm <- function(x)  A / ( 1 + exp (-((x - X0) / B)))
  
  if(plot==TRUE) {
  gp.fit <- ggplot(data = data,
                   aes(x = day,
                       y = Area.total.mm2)) +
    geom_point() +
    ggtitle(idPot) +
    stat_function(fun = fun.sigm) +
    myTheme
  gp.fit
  return(list(resultTemp = resultTemp, gp.fit = gp.fit))
  }
  return(list(resultTemp = resultTemp))
}


