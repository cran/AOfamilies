### GLM part

## define the AO transformations (symmetric and asymmetric)
ao.sym <- function(phi, verbose = FALSE)
{
  ## parameter processing
  if(verbose && phi < 0) 
  warning("sign of phi ignored as ao1(phi) = ao1(-phi)")
  #phi <- abs(phi)
  phi <- phi [phi >= 0] # positive values for ao.sym
  if(phi == 0) {
    rval <- make.link("logit")
    rval$name <- "ao.sym"
    return(rval)
  }

  linkfun <- function(mu){(2/phi) * (mu^phi - (1 - mu)^phi)/
                                    (mu^phi + (1 - mu)^phi)
  }
 
  linkinv <- function(eta){
    etastar <- pmin(2/phi - .Machine$double.eps, 
                    pmax(-2/phi + .Machine$double.eps, eta))
    if(verbose && !isTRUE(all.equal(as.vector(eta), 
	                      as.vector(etastar))))
      warning("truncation in inverse link function")
    ((1 + 0.5 * phi * etastar)^(1/phi))/
    ((1 + 0.5 * phi * etastar)^(1/phi) + 
     (1 - 0.5 * phi * etastar)^(1/phi))
  }
      
  mu.eta <- function(eta){
    phieta1 <- 1 + 0.5 * phi * eta
    phieta2 <- 1 - 0.5 * phi * eta
    phieta1^((1/phi) - 1) * 0.5/(phieta1^(1/phi) + 
    phieta2^(1/phi)) - (phieta1^(1/phi)) * (phieta1^((1/phi) - 1) * 
    0.5 - phieta2^((1/phi) - 1) * (0.5))/(phieta1^(1/phi) + 
    phieta2^(1/phi))^2
  }
  
  valideta <- function(eta) {
    if(verbose && !all(abs(0.5 * phi * eta) < 1)) 
    warning("some of the current etas are out of range")
    TRUE
  }
  
  name <- "ao.sym" # "Aranda-Ordaz symmetric"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}

ao.asym <- function(phi, verbose = FALSE)
{
  ## parameter processing
  if(phi == 1) {
    rval <- make.link("logit")
    rval$name <- "ao.asym"
    return(rval)
  }
  if(phi == 0) {
    rval <- make.link("cloglog")
    rval$name <- "ao.asym"
    return(rval)
  }

  linkfun <- function(mu) log(((1 - mu)^(-phi) - 1)/phi)
  
  linkinv <- function(eta){
    if (phi * exp(sum(eta)) <= -1) 1 else
    (1 - (1 + phi * exp(eta))^(-1/phi))
  }  
      
  mu.eta <- function(eta){
    exp(eta) * (1 + phi * exp(eta))^(-(1 + phi)/phi)
  }					  

  valideta <- function(eta){
    if(verbose && !all(phi * exp(eta) > -1)) 
    warning("some of the current etas are out of range")
    TRUE
  }
  
  name <- "ao.asym" # "Aranda-Ordaz asymmetric"

  ## return link-glm object
  structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta,
    valideta = valideta, name = name), class = "link-glm")
}

## define the central fitting function of the package for GLM
ao.glm.fit <- function (x, y, link, phi, weights, maxit = 500, ...){					

## fit model repetitively for sequence of transformation parameters 
## and store log-likelihoods
logLikvector <- rep (-1/0, length(phi)) ## to begin, define a vector 
                                        ## of log-lik = -Infinity
                                        ## whose elements are to be 
                                        ## replaced progressively								

for (i in 1:length(phi)) {
fit <- try (glm.fit (
            x, y,  weights = weights, start = rep (0, ncol(x)),
            family = binomial (link = link (phi = phi[i])), 
            ...), silent = TRUE)

## keep -Infinity as log-lik if glm procedure did not converge,
## otherwise extract log-lik from fitted model 
if (class (fit) == "try-error") {
logLikvector[i] <- -1/0
} else {
logLikvector[i] <- (fit$aic - 2 * ncol(x))/(-2)
}
}

## find the position of the MLE of the transformation parameter
ok_logLik <- logLikvector > -1/0 # record positions of 
                                 # log-lik > -Infinity 
MLE.pos <- which.max (logLikvector [ok_logLik])

## store the range of valid values for lambda and 
## corresponding log-lik
valid.lambda <- phi [which (ok_logLik)]
valid.logLik <- logLikvector [ok_logLik]

## fit model with MLE of transformation parameter
fit.MLE <- glm.fit(x, y, weights = weights,
                   start = rep (0, ncol(x)),
                   family = binomial
                   (link = link (phi = valid.lambda[MLE.pos])), 
				   ...) 

## test of asymmetry
fit.logit <- glm.fit(x = x, y = y, 
                     weights = weights,
                     start = rep (0, ncol(x)),
                     family = binomial (link = "logit"), ...)

lin.pred <- rep (0, length (y))
for (i in 1:length (fit.logit$coef)){
lin.pred <- lin.pred + fit.logit$coef[i] * x[,i]
}

theta <- exp (lin.pred)/(1 + exp(lin.pred))

r <- y * weights

u.stat <- sum (((r - weights * theta)/theta) * 
               (theta + log (1 - theta)))

i.ll <- sum (weights * ((theta + log (1 - theta))^2)/
                       (exp (lin.pred)))

i.bl <- c()
for (i in 1:length (fit.logit$coef)){
i.bl <- c(i.bl, sum ((theta + log (1 - theta)) * weights * x[,i] * 
                             (1 - theta)))
} 

i.bb <- matrix (ncol = length (fit.logit$coef), 
                nrow = length (fit.logit$coef))

for (i in 1:length (fit.logit$coef)){
for (j in 1:length (fit.logit$coef)){
i.bb[i,j] <- sum (weights * x[,i] * x[,j] * theta * (1 - theta))}
}

var.u.stat <- i.ll - t (i.bl) %*% solve (i.bb) %*% i.bl

stand.u.stat <- u.stat/var.u.stat

## rename link with numbering code				   
if (fit.MLE$family[2] == "ao.sym") {
link.bis <- 1
} else {
link.bis <- 2
}

## return list as output				   
list (link.bis = link.bis, 
      MLE = valid.lambda[MLE.pos],
      MLE.pos = MLE.pos,
      fit.MLE = fit.MLE,
      fit.MLE.coef = fit.MLE$coefficients,
      logLik = valid.logLik[MLE.pos],
      valid.lambda = valid.lambda,
      valid.logLik = valid.logLik,
      stand.u.stat = stand.u.stat,
      family = fit.MLE$family,
      fitted.values = fit.MLE$fitted.values
	  )
}

## define the generic function of the package for GLM
ao.glm <- function (x, ...) UseMethod ("ao.glm")

## define the default method of the package for GLM
ao.glm.default <- function (x, y, link, phi, weights, ...){
x <- as.matrix (x)
y <- as.numeric (y)
phi <- as.numeric (phi)
weights <- as.numeric (weights)
fit <- ao.glm.fit (x, y, link, phi, weights, ...)
class (fit) <- c ("ao.glm") 
fit
}

## define the print method of the package for GLM
print.ao.glm <- function(x, ...){
cat("Call:\n")
print(x$call)
cat("\nMLE of lambda:\n")
print(x$MLE)
cat("\nLog-likelihood associated with MLE of lambda:\n")
print(x$logLik)
cat("\nMLE coefficients for Aranda-Ordaz regression:\n")
print(x$fit.MLE.coef)
}

## define the formula method of the package for GLM
ao.glm.formula <- function (formula, data = list(), link, 
                            phi = seq (-2, 2, 0.01), weights, 
                            plotit = "TRUE", 
                            plot.spline = "TRUE", ...){
							
## keep the arguments which should go into the model frame
mf <- match.call (expand.dots = TRUE)
m <- match (c ("formula", "data", "weights"), names (mf), 0)
mf <- mf[c (1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name ("model.frame")
mf <- eval.parent (mf)

## allow model.frame to update the terms object before saving it
mt <- attr (mf, "terms")

## define response variable
y <- model.response (mf, "numeric")

## define design matrix
x <- model.matrix (mt, mf, contrasts)

## retrieve the weights (or define them if not provided) 
if (is.null (model.weights (mf))) {
weights <- rep (1, length (y)) 
} else {
weights <- model.weights (mf)
}

fit <- ao.glm.default (x, y, link, phi, weights, ...)

fit$call <- match.call()
fit$formula <- formula
fit$plotit <- plotit
fit$plot.spline <- plot.spline
fit$model <- mf
fit$design.matrix <- x
fit
}

## define the summary method of the package for GLM
summary.ao.glm <- function (object, ...){
if ( # check if we have enough data with reasonable precision 
     # to calculate a CI
-(qchisq(0.95, 1)/2) + max (object$valid.logLik) > 
min (object$valid.logLik)                  
&& 
min (abs (object$valid.logLik[1:object$MLE.pos] -
          (-(qchisq(0.95, 1)/2) +
          max (object$valid.logLik)))) < 1
&& 
min (abs (
     object$valid.logLik[object$MLE.pos:length(object$valid.logLik)] -
     (-(qchisq(0.95, 1)/2) +
     max (object$valid.logLik)))) < 1){ 	

Lower.bound <- round(object$valid.lambda[
               which.min (abs (object$valid.logLik[1:object$MLE.pos] -
                          (-(qchisq(0.95, 1)/2) +
                          max (object$valid.logLik))))], 4)                    

Upper.bound <- round(object$valid.lambda[
               object$MLE.pos + 
               which.min (abs (
			              object$valid.logLik[object$MLE.pos:
                          length(object$valid.logLik)] -
                          (-(qchisq(0.95, 1)/2) +
                          max (object$valid.logLik))))], 4)

CI.table <- data.frame(Lower.bound, Upper.bound)                          
row.names(CI.table) <- c("")
}
else{
Lower.bound <- round(min (object$valid.lambda), 4)
Upper.bound <- round(max (object$valid.lambda), 4)
CI.table <- data.frame(Lower.bound, Upper.bound)
row.names(CI.table) <- c("") 
}

## get which transformations is included in CI

epsilon <- .Machine$double.eps

if (object$link.bis == 2) { 
possible.transfo <- c ("none") 
if (epsilon >= Lower.bound && epsilon <= Upper.bound
    && 
    1 >= Upper.bound){
possible.transfo <- c("cloglog")
}
if (epsilon >= Lower.bound && 1 <= Upper.bound){
possible.transfo <- c("cloglog", "logistic")
}
if (1 >= Lower.bound && 1 <= Upper.bound
    &&
	epsilon <= Lower.bound){
possible.transfo <- c("logistic")
}
}

if (object$link.bis == 1){ 
possible.transfo <- c("none")  
if (epsilon >= Lower.bound && epsilon <= Upper.bound
    && 
    0.3955 >= Upper.bound){
possible.transfo <- c("logistic")
}
if (epsilon >= Lower.bound && 0.3955 <= Upper.bound
    && 
    0.6755 >= Upper.bound){
possible.transfo <- c("logistic", "probit")
}
if (epsilon >= Lower.bound && 0.6755 <= Upper.bound
    && 
    1 >= Upper.bound){
possible.transfo <- c("logistic", "probit", "arcsine")
}
if (epsilon >= Lower.bound && 1 <= Upper.bound){
possible.transfo <- c("logistic", "probit", "arcsine", 
                      "linear")
}
if (epsilon <= Lower.bound && 0.3955 <= Upper.bound
    && 
    0.3955 >= Lower.bound && 0.6755 >= Upper.bound){
possible.transfo <- c("probit")
}
if (0.3955 <= Lower.bound && 0.6755 <= Upper.bound
    && 
    0.6755 >= Lower.bound && 1 >= Upper.bound){
possible.transfo <- c("arcsine")
}
if (0.6755 <= Lower.bound && 1 <= Upper.bound
    && 
    1 >= Lower.bound){
possible.transfo <- c("linear")
}
}

## plot the profile log likelihood
if (object$plotit == "TRUE" ){ 
if (object$link.bis == 1){ # create title depending on type of transformation                
title <- "Aranda-Ordaz symmetric transformation  
Profile log-likelihood plot"
} 
else{
title <- "Aranda-Ordaz asymmetric transformation 
Profile log-likelihood plot"
}
 
if (object$plot.spline == "TRUE")
{
spl <- smooth.spline (object$valid.lambda, object$valid.logLik)
object$valid.lambda <- spl$x
object$valid.logLik <- spl$y
}

plot(
x = object$valid.lambda, y = object$valid.logLik,
xlab = expression (lambda), ylab = "Log-likelihood",
type = "l"
, main = title
)

## vertical line for MLE
lines (x = rep (object$valid.lambda[object$MLE.pos], 2), 
       y = c (min (object$valid.logLik), max (object$valid.logLik)),
       type = "l", lty = "dashed", col = "red")

if ( # check if we have enough data with reasonable precision to calculate a CI
-(qchisq(0.95, 1)/2) + max (object$valid.logLik) > 
min (object$valid.logLik)                  
&& 
min (abs (object$valid.logLik[1:object$MLE.pos] -
          (-(qchisq(0.95, 1)/2) +
          max (object$valid.logLik)))) < 1
&& 
min (abs (
     object$valid.logLik[object$MLE.pos:length(object$valid.logLik)] -
     (-(qchisq(0.95, 1)/2) +
     max (object$valid.logLik)))) < 1){	   
	   
## horizontal line for CI 
lines (y = rep (-(qchisq(0.95, 1)/2) + max (object$valid.logLik), 2), 
       x = c (min (object$valid.lambda), max (object$valid.lambda)), 
       lty = "dashed", col = "blue")

## vertical line for lower bound of CI
lines (x = 
       rep (object$valid.lambda[
            which.min (abs (object$valid.logLik[1:object$MLE.pos] -
                       (-(qchisq(0.95, 1)/2) +
                       max (object$valid.logLik))))], 2),
       y = c (min (object$valid.logLik), 
	          -(qchisq(0.95, 1)/2) + max(object$valid.logLik)),
       lty = "dashed", 
       col = "blue")
	   
## vertical line for upper bound of CI
lines (x = 
       rep (object$valid.lambda[object$MLE.pos + 
	                            which.min (abs 
                               (object$valid.logLik[object$MLE.pos:
                                length(object$valid.logLik)] -
                                (-(qchisq(0.95, 1)/2) +
                                max (object$valid.logLik))))], 2),
       y = c (min (object$valid.logLik), 
              -(qchisq(0.95, 1)/2) + max (object$valid.logLik)),
       lty = "dashed", 
       col = "blue")
} 
else{ 
lines (x = rep (min(object$valid.lambda), 2), 
       y = c (min (object$valid.logLik), max (object$valid.logLik)),
       type = "l", lty = "dashed", col = "blue")

lines (x = rep (max(object$valid.lambda), 2), 
       y = c (min (object$valid.logLik), max (object$valid.logLik)),
       type = "l", lty = "dashed", col = "blue")	   
}

par (xpd = NA, oma = c (4, 0, 0, 0))
legend (par ("usr")[1], par ("usr")[3],
        c ("MLE estimate", "95% Confidence Interval"),
	    col = c ("red","blue"),
        lty = c ("dashed", "dashed"),
        xjust = 0, yjust = 2.75
		)
}

call <- object$call
stand.u.stat <- object$stand.u.stat
MLE <- object$MLE
logLik <- object$logLik
plotit <- object$plotit
fit.MLE.table.coef <- summary.glm(object$fit.MLE)$coefficients
res <- list (CI.table = CI.table, call = call, 
             stand.u.stat = stand.u.stat, plotit = plotit,
             MLE = MLE, possible.transfo = possible.transfo,
             logLik = logLik, 
             fit.MLE.table.coef = fit.MLE.table.coef)
												 
class (res) <- "summary.ao.glm"
res
}

## define the print method for the summary method of the package 
## for GLM
print.summary.ao.glm <- function (x, ...){ 
cat("Call:\n")
print(x$call)
cat ("\n")

cat("Test of asymmetry:\n")
cat (paste ("Test statistic =", round (x$stand.u.stat, 3),"\n"))
cat (paste ("P-value =", round (pnorm (x$stand.u.stat), 3),"\n"))

cat("\nMLE of lambda and log-likelihood:\n")
print(data.frame(MLE = x$MLE, logLik = x$logLik, 
                 row.names = c("")))

if (!is.null(x$CI.table)){
cat("\n95% confidence interval for lambda:\n")
print (x$CI.table)

cat("\nTransformations included in confidence interval:\n")
print(x$possible.transfo)
}

cat("\nMLE coefficients for Aranda-Ordaz regression:\n")
printCoefmat (x$fit.MLE.table.coef)
cat ("\n")
}

## define the dose.p method of the package for GLM 
ao.dose.p <- function (object, cf = 1:2, p = 0.5,...) {
eta <- object$family$linkfun(p) 
b <- object$fit.MLE.coef[cf] 
x.p <- (eta - b[1L])/b[2L]
names(x.p) <- paste("p = ", format(p), ":", sep = "")
pd <- -cbind(1, x.p)/b[2L]
var.covar <- solve (t(object$fit.MLE$weights * object$design.matrix) 
                    %*% object$design.matrix)  
SE <- sqrt(((pd %*% var.covar[cf, cf]) * pd) %*% c(1, 1))
CI.low <- x.p - qnorm(0.975) * SE 
CI.high <- x.p + qnorm(0.975) * SE
res <- list(x.p = x.p, SE = SE, p = p, CI.low = CI.low, 
            CI.high = CI.high)
class(res) <- "ao.glm.dose"
res
}

## define the print method for the dose.p method of the package 
## for GLM 
print.ao.glm.dose <- function (x,...){
M <- data.frame(x$x.p, x$SE, x$CI.low,
                x$CI.high)
colnames(M) <- c("Dose", "SE", "CI.low", "CI.high")
object <- M
print(M)
}

## define the fitted method of the package for GLM 
fitted.ao.glm <- function (object,...){
fitted.values <- object$fitted.values
res <- list (fitted.values = as.numeric(fitted.values))
class (res) <- "fitted.ao.glm"
invisible(res)
}

## define the print method for the fitted method of the package 
## for GLM 
print.fitted.ao.glm <- function (x, ...){
cat (paste ("Estimate(s):\n"))
print(x$fitted.values)
}

## define the predict method of the package for GLM
predict.ao.glm <- function (object, newdata = NULL, ...) {
if(is.null(newdata)){
predicted.values <- object$fitted.values
}
else {
if (!is.null(object$formula)){
tt <- terms(object$formula)
covariates <- delete.response(tt)
m <- model.frame(covariates, newdata, ...)
x <- model.matrix(covariates, m)
} else {
x <- newdata
}

lp <- object$fit.MLE.coef %*% t(x)

predicted.values <- rep (0, length(lp))
for (i in 1:length(lp)){
if (object$link.bis == 1) { 
if (object$MLE != 0){
if (abs(0.5*lp[i]*object$MLE) < 1){
predicted.values[i] <- ((1 + (object$MLE/2)*lp[i])^(1/object$MLE))/
                       ((1 + (object$MLE/2)*lp[i])^(1/object$MLE) + 
                        (1 - (object$MLE/2)*lp[i])^(1/object$MLE))
} else {
if (0.5*lp[i]*object$MLE <= -1){
predicted.values[i] <- 0 
} else {
predicted.values[i] <- 1
} 
}
} 
else {
if (object$MLE == 0){
predicted.values[i] <- exp(lp[i])/(1+exp(lp[i]))
}
}
} else {
if (object$link.bis == 2) {
if (object$MLE != 0){
if (exp(lp[i])*object$MLE > -1){
predicted.values[i] <- 1 - ((1 + exp(lp[i]) * 
                             object$MLE)^(-1/object$MLE))
} else {
predicted.values[i] <- 1
}
}
if (object$MLE == 0){
predicted.values[i] <- 1 - exp(-exp(lp[i]))
} 
}

}
}
}

res <- list (predicted.values = predicted.values)
class (res) <- "predict.ao.glm"
invisible(res)
}

## define the print method for the predict method of the package 
## for GLM
print.predict.ao.glm <- function (x, ...){
cat (paste ("Estimate(s):\n"))
print(x$predicted.values)
}

## define the plot method of the package for GLM 
plot.ao.glm <- function (x, which = 1:4, ...){
## adjust the graphic device
if (length(which) == 1){
par(mfrow = c(1,1))
}
if (length(which) == 2){
par(mfrow = c(1,2))
}
if (length(which) == 3){
par(mfrow = c(2,2))
}
if (length(which) == 4){
par(mfrow = c(2,2))
}
if (length(which) == 5){
par(mfrow = c(3,2))
}

## asign the glm class
class(x$fit.MLE) <- c(x$fit.MLE$class, c("glm", "lm"))
plot(x$fit.MLE,which = which)
}

## define the AIC method of the package for GLM 
AIC.ao.glm <- function (object, ...){
class(object$fit.MLE) <- c(object$fit.MLE$class, c("glm", "lm"))
AIC(object$fit.MLE)
}

## define the logLik method of the package for GLM 
logLik.ao.glm <- function (object, ...){
class(object$fit.MLE) <- c(object$fit.MLE$class, c("glm", "lm"))
logLik(object$fit.MLE)
}


### QR part

## central fitting function of the package for QR
ao.qr.fit <- function (x, y, weights, kappa, 
                       phi, se, estimation, method,
					   epsilon, transfo, R,
					   ...){

## define the ao transformation                              

## theta bounded between 0 and 1
theta <- (y - (min(y)-epsilon))/
         ((max(y)+epsilon) - (min(y)-epsilon)) 						

## fit model repetitively for sequence of transformation parameters 
## and store log-likelihoods
logLikvector <- rep (-1/0, length(phi))

for (i in 1:length(phi)) {

## transformation of theta
if (transfo == "ao.sym"){ # symmetric transfo
if (phi[i] == 0){
ao <- log (theta/(1 - theta)) # define the logit transfo
} else {										
ao <- (2/phi[i]) * (theta^phi[i] - (1 - theta)^phi[i])/
      (theta^phi[i] + (1 - theta)^phi[i])
}
} 
else {
if (transfo == "ao.asym"){ # asymmetric transfo
if (phi[i] == 1){
ao <- log (theta/(1 - theta)) # define the logit transfo
} else {
if (phi[i] == 0){
ao <- log (-log(1 - theta)) # define the cloglog transfo
} else {
ao <- log((1/phi[i])*((1 - theta)^(-phi[i]) - 1))
}
}
}
}

if (estimation == "laplace"){ 
fit.seq <- try (lqm (ao ~ x - 1, 
#tau = kappa, 
weights = weights,
                ...), silent = TRUE)
} else { 
if (estimation == "linprog"){
fit.seq <- try (rq (ao ~ x - 1, tau = kappa, method = method, 
                    weights = weights, # ~ x[,2:ncol(x)]
                    ...), silent = TRUE)
} 
}

## keep -Infinity as log-lik if qr procedure did not converge,
## otherwise extract log-lik from fitted model 
if (class (fit.seq) == "try-error") {
logLikvector[i] <- -1/0
} else {

## extract predictions
if (estimation == "laplace"){
predict.seq <- predict(fit.seq, interval = TRUE, level = 0.95
)
} else { 
if (estimation == "linprog"){
predict.seq <- predict(fit.seq, method = method,
                      type = c("none"), 
                      interval = c("confidence"),
                      level = 0.95, 
newdata = list(x = x), 
                       se = se)
}
} 

min.y <- min(y)
max.y <- max(y)
eps <- epsilon

fitted <- rep (0, length(predict.seq))
for (j in 1:(length(predict.seq))){
if (transfo == "ao.sym") { 
if (phi[i] != 0){
if (abs(0.5*predict.seq[j]*phi[i]) < 1){
fitted[j] <- min.y + (((max.y + eps) - (min.y - eps)) * 
                      ((1 + (phi[i]/2)*predict.seq[j])^(1/phi[i]))/
                      ((1 + (phi[i]/2)*predict.seq[j])^(1/phi[i]) + 
                       (1 - (phi[i]/2)*predict.seq[j])^(1/phi[i])))
} else {
if (0.5*predict.seq[j]*phi[i] <= -1){
fitted[j] <- min.y 
} else {
fitted[j] <- max.y
} 
}
} else {
if (phi[i] == 0){
if (exp(predict.seq[j])> -1){
fitted[j] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - ((1 + exp(predict.seq[j]))^(-1)))		 
} else {
fitted[j] <- max.y
}
}
}
} else {
if (transfo == "ao.asym"){
if (phi[i] != 0){
if (exp(predict.seq[j])*phi[i] > -1){
fitted[j] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - ((1 + exp(predict.seq[j])*phi[i])^(-1/phi[i])))
} else {
fitted[j] <- max.y
}
}
else {
if (phi[i] == 0){
fitted[j] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - exp(-exp(predict.seq[j])))
}
}
}
}
}

fitted.data <- matrix (fitted, ncol = 3)
res <- y - fitted.data[,1]
estimated.sigma <- mean(res * (kappa - ifelse(res <= 0 , 1, 0))) #### !!!!!!!!!
logLikvector[i] <- mean (
                   dal(res, mu = 0, sigma = estimated.sigma,# tau = kappa, 
                       log = T))
}
}

## find the position of the MLE of the transformation parameter
ok.logLik <- logLikvector > -1/0 # record positions of 
                                 # log-lik > -Infinity 
MLE.pos <- which.max (logLikvector [ok.logLik])

## store the range of valid values for lambda and 
## corresponding log-lik
valid.lambda <- phi [which (ok.logLik)]
valid.logLik <- logLikvector [ok.logLik]

## fit model with MLE of transformation parameter

## transformation of theta
if (transfo == "ao.sym"){ # symmetric transfo
if (valid.lambda[MLE.pos] == 0){
ao <- log (theta/(1 - theta)) # define the logit transfo
phi.new <- 1 
transfo.new <- "ao.asym"
} else {										
ao <- (2/valid.lambda[MLE.pos]) * (theta^valid.lambda[MLE.pos] 
       - (1 - theta)^valid.lambda[MLE.pos])/
      (theta^valid.lambda[MLE.pos] + 
	  (1 - theta)^valid.lambda[MLE.pos])
phi.new <- valid.lambda[MLE.pos]
transfo.new <- "ao.sym"
}
} 
else {
if (transfo == "ao.asym"){ # asymmetric transfo
if (valid.lambda[MLE.pos] == 1){
ao <- log (theta/(1 - theta)) # define the logit transfo
} else {
if (valid.lambda[MLE.pos] == 0){
ao <- log (-log(1 - theta)) # define the cloglog transfo
} else {
ao <- log((1/valid.lambda[MLE.pos])*((1 - theta)^(-valid.lambda[MLE.pos]) - 1))
}
}
phi.new <- valid.lambda[MLE.pos]
transfo.new <- "ao.asym"
}
}

if (estimation == "laplace"){ 
fit <- try (lqm (ao ~ x - 1, #tau = kappa, 
weights = weights,
                ...), silent = TRUE)
} else { 
if (estimation == "linprog"){
fit <- try (rq (ao ~ x - 1, tau = kappa, method = method, weights = weights, 
                ...), silent = TRUE)
}
}

## extract coefficients
if (estimation == "laplace"){
fit.s <- summary (fit, R = R)
fit.coef.table <- fit.s$tTable
fit.coef <- fit.coef.table[,1]
} else {
if (estimation == "linprog"){
fit.coef.table <- summary (fit, se = se, R = R)$coef
fit.coef <- fit$coef
}
}

## extract predictions
if (estimation == "laplace"){
predict.all <- predict(fit
, interval = TRUE, level = 0.95
)
} else { 
if (estimation == "linprog"){
predict.all <- predict(fit, 
type = c("none"), 
method = method, 
interval = c("confidence"),
                 level = 0.95,
 newdata = list(x = x), 
                   se = se)
}
} 

## return list as output				   
list (fit.coef.table = fit.coef.table, 
      fit.coef = fit.coef,
	valid.lambda = valid.lambda, 
	valid.logLik = valid.logLik, 
	MLE = valid.lambda[MLE.pos], MLE.pos = MLE.pos,
	logLik.MLE = max (logLikvector [ok.logLik]),
      estimation = estimation, kappa = kappa, 
      resp = y, design.matrix = x,  
	transfo.new = transfo.new, phi.new = phi.new,
      predict.all = predict.all, 
	transfo = transfo,
	epsilon = epsilon, 
      ao = ao,  
	fit = fit
	)
}

## generic function of the package 
ao.qr <- function (x, ...) UseMethod ("ao.qr") 

## default method of the package for QR
ao.qr.default <- function (x, y, weights, method, 
                           kappa, phi, estimation, 
                           epsilon, transfo, se, R, 
                           ...){
x <- as.matrix (x)
y <- as.numeric (y)

fit <- ao.qr.fit (x, y, weights = weights, kappa = kappa, 
                  phi = phi, estimation = estimation, method = method,
	            epsilon = epsilon, se = se, R = R, 
                  transfo = transfo, ...)

class (fit) <- c ("ao.qr")
fit
}

## formula method of the package for QR
ao.qr.formula <- function (formula, data = list(), 
                           weights = rep (1, length (y)), 
                           kappa = 0.5, phi = seq(0,1.5,0.005), 
                           estimation = "laplace", 
                           epsilon = 0.001, transfo = "ao.sym",
                           plotit = "TRUE", method = "br",
				   se = "boot", R = 100, 
...){
						   
## keep the arguments which should go into the model frame
mf <- match.call (expand.dots = TRUE)
m <- match (c ("formula", "data", "weights"), names (mf), 0)
mf <- mf[c (1, m)]
mf$drop.unused.levels <- TRUE
mf[[1]] <- as.name ("model.frame")
mf <- eval.parent (mf)

## allow model.frame to update the terms object before saving it
mt <- attr (mf, "terms")

## define response variable
y <- model.response (mf, "numeric")

## define design matrix
x <- model.matrix (mt, mf, contrasts)

## retrieve the weights (or define them if not provided) 
weights <- model.weights (mf)

fit <- ao.qr.default (x, y, weights = weights, kappa, phi, 
                      estimation = estimation, method = method, 
                      epsilon = epsilon, transfo = transfo,
                      se = se, R = R, 
...)

fit$call <- match.call()
fit$formula <- formula
fit$model <- mf
fit$plotit <- plotit
fit$se <- se
fit$R <- R
fit$method <- method
fit
}	
				   
## print method of the package for QR
print.ao.qr <- function(x, ...){
cat (paste ("Quantile regression\n"))
cat (paste ("Transformation used:", x$transfo,"\n"))
cat (paste ("call:\n")) 
print (x$call)
cat ("\n")

cat (paste ("MLE and log-likelihood:\n"))
print(data.frame(MLE = x$MLE,
                 logLik = x$logLik.MLE, 
                 row.names = c("")))
cat ("\n")

cat (paste ("Coefficients:\n"))
coef.table <- matrix (
x$fit.coef, nrow = 1, ncol = length (x$fit.coef))
colnames(coef.table) <- colnames(x$design.matrix) 
row.names(coef.table) <- c("")
print (coef.table)
}				   

## summary method of the package for QR					  					   
summary.ao.qr <- function (object, ...){
call <- object$call
kappa <- object$kappa
fit.coef.table <- object$fit.coef.table
row.names(fit.coef.table) <- colnames(object$design.matrix)
phi.new <- object$phi.new
valid.lambda <- object$valid.lambda
valid.logLik <- object$valid.logLik
logLik.MLE <- object$logLik.MLE
transfo <- object$transfo
MLE <- object$MLE
plotit <- object$plotit
res <- list (call = call, phi.new = phi.new,  
             kappa = kappa, valid.lambda = valid.lambda,  
             fit.coef.table = fit.coef.table,
			 valid.logLik = valid.logLik,
			 logLik.MLE = logLik.MLE, plotit = plotit,
			 transfo = transfo, MLE = MLE
			 )
			 
class (res) <- "summary.ao.qr"
res
}

## print method for the summary method of the package 
## for QR
print.summary.ao.qr <- function (x, ...){ 
cat("Call:\n")
print(x$call)
cat ("\n")

cat (paste ("Quantile considered (tau):\n"))
print (x$kappa)
cat ("\n")

cat (paste ("MLE and log-likelihood:\n"))
print(data.frame(MLE = x$MLE,
                 logLik = x$logLik.MLE, 
                 row.names = c("")))
cat ("\n")

if (x$plotit == "TRUE") {
plot (x$valid.lambda,x$valid.logLik, type = "l",
xlab = expression (lambda),
ylab = "Log-likelihood")
lines (x = c(x$valid.lambda[which.max(x$valid.logLik)],
             x$valid.lambda[which.max(x$valid.logLik)]), 
       y = c(min(x$valid.logLik),max(x$valid.logLik)),
	   lty = "dashed", col = "red")

par (xpd = NA, oma = c (4, 0, 0, 0))
legend (par ("usr")[1], par ("usr")[3],
c ("MLE estimate"), col = c ("red"),
lty = c ("dashed"), xjust = 0, yjust = 3.5)

}

cat (paste ("Coefficients:\n"))
printCoefmat (x$fit.coef.table)
}

## predict method of the package for QR
predict.ao.qr <- function(object, newdata = NULL, ...){
min.y <- min(object$resp)
max.y <- max(object$resp)
eps <- object$epsilon
phi <- object$phi.new
predict.all <- object$predict.all
se <- object$se

if(is.null(newdata)){
fitted <- rep (0, length(predict.all))
for (i in 1:(length(predict.all))){
if (object$transfo.new == "ao.sym"){ 
if (phi != 0){
if (abs(0.5*predict.all[i]*phi) < 1){
fitted[i] <- min.y + (((max.y + eps) - (min.y - eps)) * 
                      ((1 + (phi/2)*predict.all[i])^(1/phi))/
                      ((1 + (phi/2)*predict.all[i])^(1/phi)+
                       (1 - (phi/2)*predict.all[i])^(1/phi)))
} else {
if (0.5*predict.all[i]*phi <= -1){
fitted[i] <- min.y 
} else {
fitted[i] <- max.y
} 
}
} 
} else {
if (object$transfo.new == "ao.asym"){
if (phi != 0){
if (exp(predict.all[i])*phi > -1){
fitted[i] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - ((1 + exp(predict.all[i])*phi)^(-1/phi)))
} else {
fitted[i] <- max.y
}
}
else {
if (phi == 0){
fitted[i] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - exp(-exp(predict.all[i])))
}
}
} 
}
}
}

else {
if (!is.null(object$formula)){
tt <- terms(object$formula)
covariates <- delete.response(tt)
m <- model.frame(covariates, newdata, ...)
x <- model.matrix(covariates, m)
} else {
x <- newdata
}
if (object$estimation == "laplace"){
predict.all <- predict(object$fit, newdata = x,interval = TRUE, 
                       level = 0.95)
} else {
predict.all <- predict(object$fit, type = c("none"), 
                       interval = c("confidence"), 
                       level = 0.95, 
                       newdata = list(x = x), se = se)
} 
fitted <- rep (0, length(predict.all))
for (i in 1:(length(predict.all))){
if (object$transfo.new == "ao.sym"){ 
if (phi != 0){
if (abs(0.5*predict.all[i]*phi) < 1){
fitted[i] <- min.y + (((max.y + eps) - (min.y - eps)) * 
                      ((1 + (phi/2)*predict.all[i])^(1/phi))/
                      ((1 + (phi/2)*predict.all[i])^(1/phi) +
                       (1 - (phi/2)*predict.all[i])^(1/phi)))
} else {
if (0.5*predict.all[i]*phi <= -1){
fitted[i] <- min.y 
} else {
fitted[i] <- max.y
} 
}
} 

} else {
if (object$transfo.new == "ao.asym"){
if (phi != 0){
if (exp(predict.all[i])*phi > -1){
fitted[i] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - ((1 + exp(predict.all[i])*phi)^(-1/phi)))
} else {
fitted[i] <- max.y
}
}
else {
if (phi == 0){
fitted[i] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - exp(-exp(predict.all[i])))
}
}
}
}
}
}

fitted <- matrix (fitted, ncol = 3)
colnames (fitted) <- colnames (predict.all) 
res <- list (fitted = fitted)
class (res) <- "predict.ao.qr"
invisible(res)
} 

## print method for the predict method of the package 
## for QR
print.predict.ao.qr <- function (x, ...){
cat (paste ("Estimates:\n"))
print (x$fitted)
}

## fitted method of the package 
fitted.ao.qr <- function(object, ...) {
min.y <- min(object$resp)
max.y <- max(object$resp)
eps <- object$epsilon
phi <- object$phi.new
predict.all <- object$predict.all

fitted <- rep (0, length(predict.all))
for (i in 1:(length(predict.all))){
if (object$transfo.new == "ao.sym") { 
if (phi != 0){
if (abs(0.5*predict.all[i]*phi) < 1){
fitted[i] <- min.y + (((max.y + eps) - (min.y - eps)) * 
                      ((1 + (phi/2)*predict.all[i])^(1/phi))/
                      ((1 + (phi/2)*predict.all[i])^(1/phi) + 
                       (1 - (phi/2)*predict.all[i])^(1/phi)))
} else {
if (0.5*predict.all[i]*phi <= -1){
fitted[i] <- min.y 
} else {
fitted[i] <- max.y
} 
}
} 

} else {
if (object$transfo.new == "ao.asym"){
if (phi != 0){
if (exp(predict.all[i])*phi > -1){
fitted[i] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - ((1 + exp(predict.all[i])*phi)^(-1/phi)))
} else {
fitted[i] <- max.y
}
}
else {
if (phi == 0){
fitted[i] <- min.y + ((max.y + eps) - (min.y - eps)) * 
                     (1 - exp(-exp(predict.all[i])))
}
}
}

}
}

fitted <- matrix(fitted, ncol = 3)
colnames (fitted) <- colnames (predict.all)
res <- list(fitted = fitted)
class (res) <- "fitted.ao.qr"
invisible(res)
}

## print method for the fitted method of the package 
## for QR
print.fitted.ao.qr <- function (x, ...){
cat (paste ("Estimate(s):\n"))
print(x$fitted)
}


