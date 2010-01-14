# For rcorr functions
library(Hmisc)
# For pcr functions
library(pls)
# For pcor functions
source("~/research/lib/pcor.R")

dev.out <- function(fname, fdir="../figures/", width=7, height=7, pdf.figures=T) {
	if (pdf.figures) ext = ".pdf" else ext = ".png"
	if (pdf.figures) pdf(paste(fdir,fname,ext,sep=""), width=width, height=height, family="Helvetica")
	else png(paste(fdir,fname,ext,sep=""), width=width*(480/7), height=height*(480/7))
}

ep <- function(x) {
	eval(parse(text=x))
}

catn <- function(s, sep='') {
  cat(paste(s, sep=sep),'\n')
}

barplot.ci <- function(x, param.name, ylim=NULL, ...) {
  params <- sapply(x, function(m) {coef(m)[[param.name]]})
  lower.cis <- sapply(1:length(x), function(m) {confint(x[[m]], param.name)[[1]]})
  upper.cis <- sapply(1:length(x), function(m) {confint(x[[m]], param.name)[[2]]})
  if (is.null(ylim)) {
    ylim=c(min(0,min(lower.cis)), max(upper.cis))
  }
  pts <- barplot(params, ylim=ylim, ...)
  segments(pts, upper.cis, pts, lower.cis)
}

barplot.err <- function(x, x.err, ...) {
	f.args <- list(...)
	if (is.null(f.args$ylim)) {
		f.args$ylim <- c(min(x-x.err,na.rm=T),max(x+x.err,na.rm=T))
		##print(f.args$ylim)
	}
	bp <- barplot(x, ylim=f.args$ylim, ...)
	segments(bp, x-x.err, bp, x+x.err)
}

my.axis <- function(side, at, log.at=F, log=F, expand.range=0.1, ...) {
  if (log) {
  	at <- unlist(at)
  	at <- at[at>0]
    lat <- as.integer(log10(at))
    range <- c(min(lat,na.rm=T),max(lat,na.rm=T))
    range.expansion.factor <- (range[2]-range[1])*expand.range
    range.expanded <- floor(c(range[1]-range.expansion.factor, range[2]+range.expansion.factor))
    latseq <- seq(range.expanded[1], range.expanded[2],1)
    labs <- lapply(latseq, function(m){substitute(10^i ,list(i=m))})
	place.at <- 10^latseq
    if (log.at) {
    	# Put labels at log-transformed locations.
    	place.at <- latseq
    }
    axis(side, at=place.at, labels=as.expression(labs), ...)
  }
  else {
    axis(side, at, ...)
  }
}

charlist <- function(x) {
	unlist(strsplit(x, ""))
}

## Lower triangle
lt <- function(x) {x[lower.tri(x)]}

## Lower triangle
ut <- function(x) {x[upper.tri(x)]}

## Trace
tr <- function(x) {sum(diag(x),na.rm=T)}

id.match <- function(x, y, id.x, id.y) {
	z <- match(x[id.x], y[id.y])
	data.frame(x, y[z,])
}

rcormat <- function(x, method="spearman") {
	rcorr(as.matrix(x), type=method)
}

pcor.covmat <- function(covmat) {
	cormat <- cov2cor(covmat)
	R.inv <- solve(cormat)
	-cov2cor(R.inv)
}

pcor.pc <- function(x, y, z.mat, z.indices=NULL, method='s', ...) {
  z.mat <- as.matrix(scalerankcols(z.mat))
  rc <- rcormat(z.mat, method=method)
  e <- eigen(rc$r)
  # Make scores
  scores <- z.mat %*% e$vectors
  if (is.null(z.indices)) {
    z.indices <- 1:ncol(scores)
  }
  res <- pcor.test(x, y, scores[,z.indices])
  res$pc.scores <- scores
  res$cmat <- rc$r
  res$eig <- e
  class(res) <- "pcor"
  res
}

scale.if.numeric <- function(x, ...) {
	if (class(x) == 'numeric') {
		return(scale(x,...))
	}
	else {
		return(x)
	}
}

rankcols <- function(mat, na.last="keep", ...) {
	m <- apply(mat,2,rank, na.last=na.last, ...)
	if (is.data.frame(mat)) {
		return(as.data.frame(m))
	}
	else {
		return(m)
	}
}

scalecols <- function(mat, ...) {
	return(scale.mat(mat,dim=2))
}

scale.mat <- function(mat, dim=2, ...) {
  # dim: 1=rows, 2=columns
  if (dim==1) {
	m <- t(apply(mat, dim, scale, ...))
    colnames(m) <- colnames(mat)
  }
  else if (dim==2) {
	m <- apply(mat, dim, scale, ...)
  }
  if (is.data.frame(mat)) {
    m <- as.data.frame(m)
    dimnames(m) <- dimnames(mat)
  }
  return(m)
}

scalerank <- function(x, na.last="keep", ...) {
	scale(rank(x,na.last=na.last, ...))
}

scalerankcols <- function(mat, na.last="keep", ...) {
	m <- apply(mat,2,scalerank, na.last=na.last, ...)
	if (is.data.frame(mat)) {
		return(as.data.frame(m))
	}
	else {
		return(m)
	}
}

scale.common <- function(x, ...) {
  # - Rank columns
  # - Find set of items X in common between all columns
  # - From each data subtract mean over X, divide by sd over X to scale
  mat <- as.matrix(x)
  if (dim(mat)[2] == 1) {
    m = mean(x, na.rm=T)
    sd = sd(x, na.rm=T)
    y <- (x-m)/sd
    return(y)
  }
  nobj <- dim(mat)[1]
  y <- na.omit(mat)
  means <- colMeans(y)
  sds <- apply(y, 2, sd)
  mat.scaled <- (mat - rep(means,each=nobj)) * rep(1/sds, each=nobj)
  dimnames(mat.scaled) <- dimnames(x)
  if (is.data.frame(x)) {
    return(as.data.frame(mat.scaled))
  }
  else {
    return(mat.scaled)
  }
}

# Names of, for example, pairwise correlations between variables named in ns
pair.names <- function(ns, sep='--') {
  apply(combn(ns,2), 2, function(m){paste(m,collapse=sep)})
}

fc <- function(x, digits=3) {
	format(x, digits=digits)
}

quick.formula <- function(response, predictors) {
	as.formula(paste(response,"~",paste(predictors,collapse="+")))
}

# Modified summary function for class mvr
summary.pcr <- function( g, print.out=TRUE ){
  if (g$method == 'eigen.pc') {
    nobj <- g$n
  }
  else {
    nobj <- nrow(g$scores)
  }
  factors <- attr(g$terms,"factor")
  npred <- ncol(factors)
  nresp <- nrow(factors) - npred
  yvarnames <- respnames(g)
  respname <- rownames(factors)[attr(g$terms,"response")]
  cum.compnames <- paste(1:g$ncomp, "comps")
  compnames <- names(g$Xvar)
  prednames <- rownames(g$loadings)
  if (print.out) {
    cat("Data: \tX dimension:", median(nobj,na.rm=T), npred)
    cat("\n\tY dimension:", median(nobj,na.rm=T), nresp)
    cat("\nFit method:", g$method)
    cat("\nNumber of components considered:", g$ncomp, "\n")
  }

  ## Get the % variance explained in the data by each component
  ## and the cumulative % variance explained in the response
  xve <- explvar(g)
  if (abs(sum(xve)-100) > sqrt(.Machine$double.eps)) {
    cat("\n**** Error: total X variance =", sum(xve), "> 100% ****\n")
  }
  if (g$method == 'eigen.pc') {
    cum.yve <- cumsum(100 * g$r.squared)
  }
  else {
    cum.yve <- 100 * drop(R2(g, estimate = "train", intercept = FALSE)$val)
  }
  h <- rbind(cumsum(xve), cum.yve)
  dimnames(h) <- list(c("X", respname), compnames)

  ## number of predictors
  npred <- length(prednames)
  ## percentage differences
  h.offset <- h[,c(1,1:(g$ncomp-1))]
  h.offset[,1] <- 0
  h.diffs <- h - h.offset
  yve <- h.diffs[2,]

  ## print summary with cumulative R^2
  if (print.out) {
    cat( "\nCumulative variance explained:\n" )
    print( round( h, 2 ), print.gap=2 )
  }

  ## print summary with percentage differences
  if (print.out) {
    cat( "\nPercentage differences:\n" )
    print( round( h.diffs, 2 ), print.gap=2 )
  }

  ## Print the normalized projection
  ## Order according to the contributions to the first component
  p2 <- g$projection^2
  ord <- order(p2[,1], decreasing=T)
  if (abs(sum(p2)/g$ncomp - 1) > 10*.Machine$double.eps) {
    ## DAD: warning() here instead?
    if (print.out) cat( "\nError: non-normalized projection matrix (diff = ", abs(sum(p2)/g$ncomp-1),")",sep='' )
    p2 <- p2/(sum(p2)/g$ncomp)
  }
  if (print.out) {
    cat( "\nPercentage contributions to components:\n" )
    print( round( p2[ord,]*sign(g$projection), 3 ) )
  }

	## Print total variance contributions
	if (g$method == 'eigen.pc') {
		## = r^2 for each variable...
		covmat <- g$covmat
		cm <- cov2cor(g$covmat)
		r.squared = drop(cm[respname, prednames])^2
	}
	else {
		all.comps <- dim(g$fitted.values)[3]
		# Recreate the predictors and the response
		d <- data.frame(y=g$fitted.values[,,all.comps]+g$residuals[,,all.comps], g$scores %*% t(g$projection))
		colnames(d) <- c(respname,prednames)
		r.squared <- as.vector(corr.list(d, respname, prednames)^2)
	}
	p2 <- g$projection^2
	comp.r2 <- rep(yve, each=nrow(p2))
	dim(comp.r2) <- dim(p2)
	dimnames(comp.r2) <- dimnames(p2)
	pred.ve <- comp.r2 * p2
	## Construct bars from principal components, suitable for barplot(bars)
	pred.bars <- t(pred.ve)
	bars <- pred.ve
	if (print.out) {
		cat( "\nPercentage variance explained by each predictor:\n")
		print( round(colSums(pred.ve),3) )
		print( round(pred.ve, 3) )
	}

	res <- list(bars=bars, bars.pc=bars, bars.pred=pred.bars, proj=g$projection, Yvar=h.diffs[2,], Xvar=h.diffs[1,], n=nobj, ncomp=g$ncomp)
	invisible(res)
}

pcr.summary <- summary.pcr
summary.mvr <- summary.pcr

# Randomly chooses a unit vector on an n-sphere.
unitvec <- function(n) {
	x <- rnorm(n)
	x/sqrt(sum(x^2))
}

my.pcr <- function(form, data, fxn=scalerankcols, ...) {
	vars <- rownames(attr(terms(form),"factors"))
	g <- pcr(form, data=fxn(data[vars]), ...)
	g
}

my.plsr <- function(form, data, fxn=scalerankcols, ...) {
	vars <- rownames(attr(terms(form),"factors"))
	g <- plsr(form, data=fxn(data[vars]), ...)
	g
}

geom.mean <- function(x, na.rm=FALSE) {
	if (na.rm) {
		x <- na.omit(x)
	}
	exp(mean(log(x)))
}

reweighted.mean <- function(x,w,p=FALSE,...) {
  x <- c(x)
  w <- c(w)
  w[is.na(x)] <- NA
  w <- w/sum(w,na.rm=T)
  sum(w*x,na.rm=T)
}

var.wtd.mean.cochran <- function(x,w)
#
#	Computes the variance of a weighted mean following Cochran 1977 definition
#
{
	n = length(w)
	xWbar = wtd.mean(x,w)
	wbar = mean(w)
	out = n/((n-1)*sum(w)^2)*(sum((w*x-wbar*xWbar)^2)-2*xWbar*sum((w-wbar)*(w*x-wbar*xWbar))+xWbar^2*sum((w-wbar)^2))
	return(out)
}

my.wtd.var <- function(x, w, normwt=TRUE, na.rm=TRUE) {
	if (length(x[!is.na(x)])==1) {
		return(0)
	}
	return(wtd.var(x,w,normwt,na.rm))
}

# Based on wtd.mean and wtd.var from Frank Harrell's Hmisc package
# Simultaneous calculation of mean and var
mat.wtd.meanvar <- function (x, wts = NULL, normwt = FALSE, na.rm = TRUE) {
	if (!length(wts)) {
		return(list(mean=apply(x,1,mean, na.rm=na.rm), var=apply(x,1,var, na.rm=na.rm)))
	}
	weight.matrix <- rep(wts, each=nrow(x))
	dim(weight.matrix) <- dim(x)
	if (na.rm) {
	    s <- is.na(x + weight.matrix)
	    x[s] <- NA
	    weight.matrix[s] <- NA
	}
	if (normwt) {
		rs <- rowSums(weight.matrix, na.rm=na.rm)
		lengths <- rowSums(!is.na(x))
		weight.matrix <- weight.matrix*rep(lengths, ncol(x))/rep(rs, ncol(weight.matrix))
		xbar <- rowSums(weight.matrix * x, na.rm=na.rm)
		rs <- rowSums(weight.matrix, na.rm=na.rm)
	}
	else {
		rs <- rowSums(weight.matrix, na.rm=na.rm)
		xbar <- rowSums(weight.matrix * x, na.rm=na.rm)/rs
	}
	xvars <- rowSums(weight.matrix * ((x - rep(xbar, ncol(x)))^2), na.rm=na.rm)/(rs - 1)
	list(mean=xbar, var=xvars)
}

mat.wtd.mean <- function(x, w, na.rm=TRUE) {
	weight.matrix <- rep(w, each=nrow(x))
	dim(weight.matrix) <- dim(x)
	weight.matrix[is.na(x)] <- NA
	rs <- rowSums(weight.matrix, na.rm=TRUE)
	# Normalize the weights, accounting for NAs
	weight.matrix <- weight.matrix/rs
	weighted.x <- weight.matrix*x
	x.mean <- rowSums(weighted.x, na.rm=na.rm)
	x.mean
}

se <- function(x, w) {
	mv <- mat.wtd.meanvar(x, wts=w)
	print(w)
	res = sqrt(sum(mv$var, na.rm=TRUE))
	res
}

get.weights <- function(x) {
	## Find weights that minimize the standard error
	w <- rep(1/ncol(x),ncol(x))
	g <- optim(par=w, fn=se, x=x, method="L-BFGS-B", lower=c(0,0,0,0))
	g$par
}

mean.measurement <- function(x, wts=NULL, na.rm=FALSE) {
	if (is.null(wts)) {
		wts <- seq(1,1,length.out=ncol(x))/ncol(x)
	}
	# For missing values, if na.rm=T, reweight
	#mv <- mat.wtd.meanvar(x, wts, normwt=TRUE, na.rm=na.rm)
	x.mean <- apply(x, 1, mean, na.rm=na.rm) #mv$mean #mat.wtd.mean(x, wts, na.rm=na.rm) #apply(x, 1, wtd.mean, weights=wts, na.rm=na.rm)
	x.var <- apply(x, 1, var, na.rm=na.rm)  #mv$var
	# Make all rows which contain only NA's equal NA, not 0
	if (na.rm) {
		x.mean[apply(x, 1, function(m) {all(is.na(m))})] <- NA
		x.var[apply(x, 1, function(m) {all(is.na(m))})] <- NA
	}
	list(mean=x.mean, var=x.var)
}

combine.model <- function(x, meas.names, wts=NULL, na.rm=FALSE) {
	mean.data <- list()
	var.data <- list()
	for (meas.name in names(meas.names)) {
		meas <- meas.names[[meas.name]]
		data.x <- as.matrix(x[,meas])
        mm <- mean.measurement(data.x, wts=wts[meas], na.rm=na.rm)
        mean.data[[meas.name]] <- mm$mean
        var.data[[meas.name]] <- mm$var
	}
	list(means=as.data.frame(mean.data), vars=as.data.frame(var.data))
}

mean.model <- function(x, meas.names, wts=NULL, na.rm=FALSE, scale=FALSE) {
	mean.data <- NULL
	for (meas.name in names(meas.names)) {
		meas <- meas.names[[meas.name]]
		data.x <- as.matrix(x[,meas])
        if (scale) {
          data.x <- scale.common(data.x)
        }
        mean.data[[meas.name]] <- mean.measurement(data.x, wts=wts[meas], na.rm=na.rm)
	}
	as.data.frame(mean.data)
}

rcor.grid <- function(x, y, na.rm=FALSE, method='s') {
	res <- list()
	res$r <- matrix(0,nrow=ncol(x),ncol=ncol(y))
	res$n <- matrix(0,nrow=ncol(x),ncol=ncol(y))
	# Ignore NAs from each pair individually
	name.pairs <- as.matrix(expand.grid(1:ncol(x), 1:ncol(y)))
	for (i in 1:nrow(name.pairs)) {
		m <- name.pairs[i,]
		d <- na.omit(data.frame(x=x[,m[[1]]], y=y[,m[[2]]]))
		res$r[m[[1]],m[[2]]] <- cor(d$x, d$y, method=method)
		res$n[m[[1]],m[[2]]] <- nrow(d)
	}
	dimnames(res$r) <- list(colnames(x), colnames(y))
	dimnames(res$n) <- dimnames(res$r)
	class(res) <- "rcov"
	res
}

# DAD: currently does not properly handle rcov(x,y) with x and y matrices.
# DAD: fix
rcov <- function(x, y=NULL, na.rm=FALSE) {
	res <- list()
	if (is.null(y)) {
		res$r <- matrix(0,nrow=ncol(x),ncol=ncol(x))
		res$n <- matrix(0,nrow=ncol(x),ncol=ncol(x))
		if (na.rm) {
			x <- na.omit(x)
			res$r <- cov(x)
			res$n <- res$n + nrow(as.matrix(x))
		}
		else {
			# Ignore NAs from each pair individually
			name.pairs <- as.matrix(expand.grid(1:ncol(x), 1:ncol(x)))
			for (i in 1:nrow(name.pairs)) {
				m <- name.pairs[i,]
				y <- na.omit(x[,m])
				res$r[m[[1]],m[[2]]] <- cov(y[,1], y[,2])
				res$n[m[[1]],m[[2]]] <- nrow(y)
			}
		}
		dimnames(res$r) <- list(colnames(x), colnames(x))
		dimnames(res$n) <- dimnames(res$r)
	}
	else {
      ## Loop over the pairs
      ## If x and y are vectors, return a scalar r and n.
		d <- data.frame(x=x,y=y)
		#colnames(d) <- c(deparse(substitute(x)), deparse(substitute(y)))
		#if (na.rm) {
			d <- na.omit(d)
		#}
		res$r <- cov(d$x, d$y)
		res$n <- nrow(d)
	}
    class(res) <- "rcov"
	res
}

rcov2cor <- function(x) {
  rc <- x
  rc$r <- cov2cor(x$r)
  rc
}

print.rcov <- function(x) {
  cat("r\n")
  print(round(x$r,2))
  cat("\nn\n")
  print(x$n)
}

rcor <- function(x, y=NULL, ...) {
	if (is.null(y)) {
		res <- rcormat(x, ...)
	}
	else {
		rc <- rcorr(x,y,...)
		res <- NULL
		res$r <- rc$r['x','y']
		res$n <- rc$n['x','y']
	}
	return(res)
}



cov.estimate <- function(x, meas.names, wts=NULL, na.rm=TRUE, use="pairwise.complete.obs", cov.fxn=rcov, trans.fxn=NULL, regularize=TRUE) {
	# Create
	# meas.names is a list of measurement sets, each named as in: list(abund=c("abund.1","abund.2"), expr=c("expr.2","expr.4"),...)
	# data is a data.frame with columns corresponding to the names in meas.names
	# wts [optional] is a list containing the error variances for each measurement named in meas.names
	x <- as.matrix(x)
	all.f <- names(meas.names)

	use <- match.arg(use, c("everything","all.obs","complete.obs","pairwise.complete.obs"))
	if (use == 'complete.obs') {
		# Eliminate cases where all variables do not have at least two valid measurements
		meas.na <- function(y) {
			s <- sapply(meas.names, function(m) {sum(!is.na(y[m])) >= 2})
			all(s)
		}
		at.least.two <- apply(x,1,meas.na)
		x <- x[at.least.two,]
	}

	if (is.null(wts)) {
		# Weight each measurement equally
		for (nm in all.f) {
			for (m in meas.names[[nm]]) {
				wts[[m]] <- 1/length(meas.names[[nm]])
			}
		}
	}
	# Compute the means of the data
	comb <- combine.model(x, meas.names, wts, na.rm=na.rm)
	mean.data <- comb$means
	mean.data.xform <- mean.data
	var.data <- comb$vars


	# Transform data after means have been computed
	if (!is.null(trans.fxn)) {
		x <- apply(x, 2, trans.fxn)
		mean.data.xform <- as.data.frame(apply(mean.data, 2, trans.fxn))
	}

	# Compute the covariance matrix
	# s_ij = {
	# 	i!=j : S_{\bar{X_i},\bar{X_j}} with \bar{X_i} = \sum_k \alpha_k X_{ik} and \sum_k \alpha_k = 1
	# 	i==j : \sum_{m<n} \beta_{mn} S_{X_{im}, X_{in}} with \sum_{m<n} \beta_{mn} = 1
	# }
	#
	# Use \alpha_k = 1/n and \beta_{kl} = 2/(n(n-1)), unless weights are provided
	n.vars <- length(all.f)
    cov.Z <- list()
	# Upper triangle first
	cov.Z$r <- matrix(0,nrow=n.vars, ncol=n.vars,dimnames=list(all.f, all.f))
	cov.Z$n <- matrix(0,nrow=n.vars, ncol=n.vars,dimnames=list(all.f, all.f))
	# Diagonal, i==j
	for (var.i in 1:length(all.f)) {
		meas = all.f[var.i]
        vars <- meas.names[[meas]]
        meas.x <- x[,vars]
        if (length(vars) < 2) {
			## Handles case of singleton measurements -- note that these will not be unbiased estimates!
			rc <- cov.fxn(meas.x, meas.x)
			cov.Z$r[var.i, var.i] <- rc$r
			cov.Z$n[var.i, var.i] <- rc$n
        }
        else {
        	## More than one measurement.
			ps.names <- combn(meas.names[[meas]],2)
			n.combs <- ncol(ps.names)
			all.covs <- apply(ps.names, 2, function(m) {
				meas.x <- x[,m[[1]]]
				meas.y <- x[,m[[2]]]
				rc <- cov.fxn(meas.x, meas.y)
				# Weights
				w <- wts[[m[[1]]]]*wts[[m[[2]]]]
				c(rc$r, rc$n, w) #, m[[1]], m[[2]])
			})
			dimnames(all.covs) <- list(c("r","n","w"), NULL)
			all.covs <- as.data.frame(t(all.covs))
			cov.Z$r[var.i, var.i] <- wtd.mean(all.covs$r, weights=all.covs$w)
			cov.Z$n[var.i, var.i] <- wtd.mean(all.covs$n, weights=all.covs$w)
		}
	}
	# Off-diagonal, i<j
    # Here, the covariance between two noisy variables is an unbiased estimate of the true covariance, assuming additive uncorrelated noise
    # Our best estimate of the true value of the variables is given by the means.
	ps <- combn(1:n.vars,2)
	ps.names <- combn(colnames(mean.data),2)
	meas <- meas.names[all.f]
	for (i in 1:ncol(ps)) {
		row <- ps[1,i]
		col <- ps[2,i]
		x.mean <- mean.data.xform[,ps.names[1,i]]
		y.mean <- mean.data.xform[,ps.names[2,i]]
        rc <- cov.fxn(x.mean, y.mean)
 		cov.Z$r[row,col] <- rc$r
        cov.Z$n[row,col] <- rc$n
	}
    # We've computed the diagonal and the upper triangle;
    # now make symmetric.
	vars <- diag(cov.Z$r)
	cov.Z$r <- cov.Z$r + t(cov.Z$r)
	diag(cov.Z$r) <- vars

	# Now we're done computing the raw covariance matrix.
	# We may wish to ensure that it's positive semi-definite.
	cov.Z.unc <- cov.Z$r
	if (regularize) {
		regularize.bock <- T
		min.eig.ratio <- 1e4
		if (regularize.bock) {
			# Bock & Petersen Biometrika 1975
			orig.vars <- apply(mean.data.xform, 2, var, na.rm=T)
			# St is the unbiased estimate of the covariance matrix
			St <- cov.Z$r
			d.Se <- pmax(orig.vars - diag(cov.Z$r),0) #max(orig.vars)/min.eig.ratio)
			# Se is the error estimate, such that St = Sy - Se
			Se <- diag(d.Se)
			Sy = St + Se
			# DAD: more principled way to choose large.small.ev.ratio?
			S.corr <- posdef.bock(Sy, Se, large.small.ev.ratio=min.eig.ratio)
			cov.Z$r <- S.corr
		}
		else {
			cov.Z$r <- posdefify(cov.Z$r, eps.ev=1/min.eig.ratio)

		}
	}

    ns <- diag(cov.Z$n)
    cov.Z$n <- cov.Z$n + t(cov.Z$n)
    diag(cov.Z$n) <- ns
    cov.Z$n <- round(cov.Z$n)
	dimnames(cov.Z$r) <- list(all.f, all.f)
	dimnames(cov.Z$n) <- list(all.f, all.f)
    class(cov.Z) <- "rcov"

    res <- list(covmat=cov.Z, means=mean.data, vars=var.data, means.xform=mean.data.xform, covmat.uncorr=cov.Z.unc)
	res
}

cov.test.generator <- function(n.obs, n.vars, n.meas, sds=NULL) {
	if (is.null(sds)) {
		sds <- runif(n.vars, 1, 2)
	}
	d <- NULL
	u <- rnorm(n.obs, mean=0, sd=mean(sds))
	xu <- as.matrix(sapply(sds, function(sd) {u + rnorm(n.obs, sd=mean(sds))}))
	colnames(xu) <- paste("v",1:n.vars,sep='.')
	meas.names <- vector("list",n.vars)
	names(meas.names) <- paste("v",1:n.vars,sep='.')
	for (i in 1:n.meas) {
		xi <- as.matrix(sapply(1:n.vars, function(j) {xu[,j] + rnorm(n.obs,sd=sds[j])}))
		colnames(xi) <- paste("v",1:n.vars,i,sep='.')
		d <- cbind(d,xi)
	}
	for (n in names(meas.names)) {
		meas.names[[n]] <- paste(n,1:n.meas,sep='.')
	}
	res <- list(meas.data=d, meas.names=meas.names, true.data=xu, u=u)
	res
}

equalize.noise <- function(meas.data, meas.names, true.data) {
	# Create new meas.data such that r(meas.data, true.data) is equal for all variables
	#r <- list()
	true.vars <- apply(true.data, 2, var)
	noise.equalized.data <- NULL
	min.corr <- 1
	for (m in names(meas.names)) {
		mn <- unclass(meas.names[[m]])
		d <- data.frame(meas.data[,mn], y=true.data[,m])
		colnames(d) <- c(mn, "y")
		cl <- corr.list(d, "y", mn)
		if (min(cl) < min.corr) {
			min.corr <- min(cl)
		}
	}
	#print(min.corr)
	for (m in names(meas.names)) {
		mn <- unclass(meas.names[[m]])
		d <- data.frame(meas.data[,mn], y=true.data[,m])
		colnames(d) <- c(mn, "y")
		# add gaussian noise with sd = sqrt(1/(var(true)*q^2) - var(true)) to equalize noise
		xi <- as.matrix(sapply(mn, function(q) {
			x <- meas.data[,q];
			norm.sd <- sqrt(true.vars[m]*(1/(min.corr^2) - 1))
			true.data[,m] + rnorm(length(x), sd=norm.sd)
		}))
		noise.equalized.data <- cbind(noise.equalized.data, xi)
	}
	noise.equalized.data[,colnames(meas.data)]
}

cov.estimate.reps <- function(data.reps, scale=TRUE, na.rm=TRUE, boot.covars=FALSE, boot.reps=100, trans.fxn=noop) {
	# data.reps is an n x m x r matrix -- n observations, m variables, r replicate measurements of each variable's observations
	d <- dim(data.reps)
	dn <- dimnames(data.reps)
	flat.data <- NULL
	meas.names <- vector("list", d[3]);
	names(meas.names) <- dn[2]
	for (n in dn[2]) {
		meas.names[n] <- paste(dn[2],1:d[3],sep='')
	}

	for (i in 1:d[3]) {
		m <- data.reps[,,i]
		colnames(m) <- paste(dn[2],i,sep='')
		flat.data <- cbind(flat.data,m)
	}
	cov.estimate(flat.data, meas.names, scale, na.rm, boot.covars, boot.reps, trans.fxn)
}


pcr.covmat.full <- function(form, covmat, n=NA, data=NULL, stripped=FALSE, regularize=FALSE, ...) {
  ## Principal component regression from a pre-computed covariance matrix
  ## form = a formula describing the regression equation.  Presently only understands a single response and purely additive predictors, e.g. y~x1+x2+x3
  ## covmat = a covariance matrix describing the data
  ## n = an integer or matrix describing the number of observations used to generate the observations in covmat
  ## data = the input data; required if stripped==TRUE
  ## stripped = don't use data, don't generate scores, fitted values, or residuals.
  if (is.null(data)) {
  	# We can't produce scores, fitted values, or residuals without the data.
  	stripped = T
  }
  if (!stripped) {
	data <- as.matrix(data)
  }
  flds <- rownames(attr(terms(form),"factors"))
  resp.f <- flds[attr(terms(form),"response")]
  inds <- 1:length(flds)
  ## Predictors are all the terms except the response
  pred.f <- flds[inds[flds!=resp.f]]
  all.f <- c(resp.f, pred.f)

  ## Dimensions of the data
  ncomp <- length(pred.f)
  if (!stripped) {
	nobj <- nrow(data)
  }
  else {
    nobj <- mean(n, na.rm=T)
  }
  nresp <- length(resp.f)
  npred <- length(pred.f)

  ## Prepare the data
  if (!stripped) {
	all.x <- as.matrix(data)
	## Means
	means=colMeans(all.x, na.rm=T)
    Xmeans <- means[pred.f]
	Ymean <- means[resp.f]
	resp <- all.x[,resp.f] - rep(Ymean, each=nobj)
	pred <- all.x[,pred.f] - rep(Xmeans, each=nobj)
  }

  ## Extract the relevant covariance matrix
  cov.Z <- covmat[all.f, all.f]
  cor.Z <- cov2cor(cov.Z)

  ## Now extract eigen decomposition of correlation matrix to perform PCA
  ## on the predictors.
  cov.X <- cov.Z[pred.f, pred.f]
  cor.X <- cor.Z[pred.f, pred.f]
  eig <- eigen(cor.X)
  evec <- eig$vector

  compnames <- paste("Comp.", 1:ncomp,sep='')
  cum.compnames <- paste(1:ncomp, "comps")
  dimnames(evec) <- list(pred.f, compnames)

  ## Set up projection matrix
  P <- matrix(nrow=nrow(evec)+1, ncol=ncol(evec)+1, 0)
  P[1,1] <- 1
  P[2:nrow(P),2:ncol(P)] <- evec

  ## Covariance: C(ZP) = P^T C(Z) P
  ## Should always have  cor(ZP)[2:n,2:n] = I after rotation
  cov.ZP <- t(P) %*% cov.Z %*% P
  # Must be true that cov(X,Y) <= max{ cov(X,X), cov(Y,Y) }
  # Enforce maximum if it is exceeded.
  vars <- diag(cov.ZP)
  max.cov.ZP <- sqrt(vars %*% t(vars))
  cor.ZP <- cov2cor(cov.ZP)

  dimnames(cor.ZP) <- list(c(resp.f, compnames),c(resp.f, compnames))

  ## Variance explained in the predictors
  Xvar <- eig$values
  names(Xvar) <- compnames

  ## with X = predictors; uX = predictor means
  ##      D = eigenvectors of cor(X)
  ##      y = response; uY = response mean
  ##
  ## alpha ~ D^T ([X-uX]^T [X-uX])^{-1} (X-uX)^T (y-uY)
  ##       = D^T cov(X,X)^{-1} cov(X,y)
  alpha.hat <- t(evec) %*% solve(cov.X) %*% cov.Z[pred.f, resp.f]
  dimnames(alpha.hat) <- list(compnames, resp.f)

  ##beta.hat <- evec %*% alpha.hat
  ##dimnames(beta.hat) <- list(pred.f, resp.f)

  r.squared <- cor.ZP[resp.f, compnames]^2
  if (regularize) {
  	# For eigenvalues more than eig.thresh times less than
  	# the largest eigenvalue, set R^2 to NA.
  	eig.thresh = 1e-6
  	r.squared[eig$values/eig$values[1]<eig.thresh] <- NA
  }
  names(r.squared) <- compnames

  ## Scores
  comp.scores <- NULL
  y.hat <- NULL
  if (!stripped) {
      comp.scores <- pred %*% evec
      colnames(comp.scores) <- compnames
    }

  ## Loadings
  loadings <- evec
  dimnames(loadings) <- list(pred.f, compnames)
  class(loadings) <- "loadings"
  Yloadings <- t(alpha.hat)
  dimnames(Yloadings) <- list(resp.f, compnames)
  class(Yloadings) <- "loadings"

  ## Coefficients
  coefs <- array(0, dim = c(npred, nresp, ncomp))
  for (i in 1:ncomp) {
    coefs[,,i] <- evec[,1:i,drop=F] %*% alpha.hat[1:i,,drop=F]
  }
  dimnames(coefs) <- list(pred.f, resp.f, cum.compnames)

  ## Residuals and fitted values
  if (!stripped) {
    ## Fitted values
    y.hat <- array(0, dim = c(nobj, nresp, ncomp))
    for (i in 1:ncomp) {
      ## Reconstruction of correlation matrix using first i eigenvectors
      v <- evec * 0
      v[,1:i] <- evec[,1:i]
      a.i <- alpha.hat * 0
      a.i[1:i,] <- alpha.hat[1:i,]
      y.i <- pred %*% v %*% a.i
      y.hat[,,i] <- y.i + Ymean
    }
    dimnames(y.hat) <- list(NULL,resp.f,cum.compnames)
    ## Residuals
    residuals <- resp - y.hat
  }


  z <- list(call=match.call(), method='eigen.pc', coefficients=coefs, loadings=loadings, Yloadings=Yloadings, projection=unclass(loadings),
            Xvar=Xvar, Xtotvar=sum(unlist(Xvar)), ncomp=ncomp, terms=terms(form), n=n, covmat=covmat)
  if (!stripped) {
    z$Xmeans <- Xmeans
    z$Ymeans <- Ymean
    z$scores <- comp.scores
    z$fitted.values <- y.hat
    z$residuals <- residuals
  }
  class(z) <- "mvr"
  z$cov.ZP <- cov.ZP
  z$cor.ZP <- cor.ZP
  z$r.squared <- r.squared
  z$alpha.hat <- alpha.hat
  z$eig <- eig
  z
}

pcr.covmat <- function(form, covmat, n=NA, ncomp=NULL) {
	## Principal component regression from a pre-computed covariance matrix.
	## form = a formula describing the regression equation.  Presently only understands a single response and purely additive predictors, e.g. y~x1+x2+x3.
	## covmat = a covariance matrix describing the data.
	## n = an integer or matrix describing the number of observations used to generate each of the covariances in covmat; n is not used but is passed into the results.
	cormat <- cov2cor(covmat)
	flds <- rownames(attr(terms(form),"factors"))
	resp.f <- flds[attr(terms(form),"response")]
	inds <- 1:length(flds)
	## Predictors are all the terms except the response
	pred.f <- flds[inds[flds!=resp.f]]
	all.f <- c(resp.f, pred.f)
	if (is.null(ncomp)) {
		ncomp <- length(pred.f)
	}

	## Dimensions of the data
	nresp <- length(resp.f)
	npred <- length(pred.f)
	nobj <- n

	## Names
	compnames <- paste("Comp.", 1:ncomp,sep='')
	cum.compnames <- paste(1:ncomp, "comps")

	## Extraction of correlations between predictors and response
	cor.X <- cormat[pred.f, pred.f]
	cor.Xy <- cormat[pred.f, resp.f]
	cor.X.inv <- solve(cor.X)

	## Eigen decomposition of the predictor correlation matrix
	eig <- eigen(cor.X)
	U <- eig$vectors[,1:ncomp]
	dimnames(U) <- list(pred.f, compnames)

	## Variance explained in the predictors
	Xvar <- eig$values[1:ncomp]
	names(Xvar) <- compnames

	## Loadings
	loadings <- U[, 1:ncomp]
	dimnames(loadings) <- list(pred.f, compnames)
	class(loadings) <- "loadings"
	## Regression equation is y = X U alpha + e
	## Regression equation is y = X beta + e
	beta.hat <- cor.X.inv %*% cor.Xy
	alpha.hat <-  t(U) %*% beta.hat
	Yloadings <- t(alpha.hat)
	dimnames(Yloadings) <- list(resp.f, compnames)
	class(Yloadings) <- "loadings"

	## Coefficients
	## These are the beta in y = X beta + e
	## which, with the projected regression y = X U alpha + e,
	## gives beta = U alpha.  Following the pls/mvr convention, we
	## generate cumulative betas, including components up to i.
	coefs <- array(0, dim = c(npred, nresp, ncomp))
	for (i in 1:ncomp) {
		coefs[,,i] <- U[,1:i,drop=F] %*% alpha.hat[1:i,,drop=F]
	}
	dimnames(coefs) <- list(pred.f, resp.f, cum.compnames)

	## Coefficients of determination, R^2
	## Each component (XU)i, represented by the ith eigenvector U[,i], has variance
	## equal to the corresponding ith eigenvalue.
	## The component R^2 are given by cor(X U_i, y)^2 = cov(X U_i, y)^2/ eigenvalue_i = (t(U_i) * cov(X,y))^2/eigenvalue_i
	r.squared <- sapply(1:ncomp, function(i) { (t(U[,i]) %*% cor.Xy)^2/eig$values[i] })
	names(r.squared) <- compnames

	z <- list(call=match.call(), method='eigen.pc', coefficients=coefs, loadings=loadings, Yloadings=Yloadings, projection=U,
			Xvar=Xvar, Xtotvar=sum(unlist(Xvar)), ncomp=ncomp, terms=terms(form), r.squared=r.squared, n=n, covmat=covmat, eig=eig)
	class(z) <- "mvr"
	z
}

test.pcr.covmat <- function(method=c("pearson","spearman")) {
  ## Test the pcr.covmat function
  ## Ensure that it produces identical values to pcr()
  n <- 10000
  f <- function(n, method) {
    y <- rnorm(n)
    x <- sapply(1:5, function(m){y+rnorm(y,sd=m)})
    d <- data.frame(y,x);
    pred.f <- colnames(d)[2:6]
    colnames(x) <- pred.f
    resp.f <- "y"
    form <- as.formula(paste(resp.f,"~", paste(pred.f,collapse="+")))

    xform <- switch(method, pearson=scalecols, spearman=rankcols)
    g.svd <- pcr(form, data=xform(d))
    rc <- rcormat(d, meth=method)
    g.cov <- pcr.covmat(form, covmat=rc$r, n=rc$n)
    cor.X <- rc$r[pred.f,pred.f]
    cor.Xy <- rc$r[pred.f,resp.f]
    e <- eigen(cor.X)
    u <- e$vectors
    b.svd <- coef(g.svd) #t(u) %*% solve(cor.X) %*% cor.Xy
    b.cov <- coef(g.cov)
    list(svd=b.svd, cov=b.cov)
  }
  sq.diff <- function(n, method) {
    x <- f(n, method)
    (x$svd-x$cov)^2
  }
  method <- match.arg(method, c("pearson","spearman"))
  diffs <- replicate(10, sq.diff(1000, method))
  cat("Cumulative RMS deviation between SVD and Cov methods (", method, ") = ", sqrt(mean(diffs)), "\n", sep='')
}

print.mvr <- function (x, ...) {
    switch(x$method, kernelpls = {
        regr = "Partial least squares"
        alg = "kernel"
    }, widekernelpls = {
        regr = "Partial least squares"
        alg = "wide kernel"
    }, simpls = {
        regr = "Partial least squares"
        alg = "simpls"
    }, oscorespls = {
        regr = "Partial least squares"
        alg = "orthogonal scores"
    }, svdpc = {
        regr = "Principal component"
        alg = "singular value decomposition"
    }, eigen.pc = {
        regr = "Principal component"
        alg = "eigenvalue decomposition"
    }, stop("Unknown fit method."))
    cat(regr, "regression, fitted with the", alg, "algorithm.")
    if (!is.null(x$validation))
        cat("\nCross-validated using", length(x$validation$segments),
            attr(x$validation$segments, "type"), "segments.")
    cat("\nCall:\n", deparse(x$call), "\n", sep = "")
    invisible(x)
}


# From the package sfsmisc
posdefify <- function (m, method = c("someEVadd", "allEVadd"), symmetric, eps.ev = 1e-07) {
    stopifnot(is.numeric(m) && is.matrix(m))
    method <- match.arg(method)
    ev.m <- eigen(m, symmetric = symmetric)
    n <- length(lam <- ev.m$values)
    Eps <- eps.ev * abs(lam[1])
    if (lam[n] < Eps) {
        switch(method, someEVadd = lam[lam < Eps] <- Eps, allEVadd = lam <- lam + Eps - lam[n])
        Q <- ev.m$vectors
        o.diag <- diag(m)
        m <- Q %*% (lam * t(Q))
        D <- sqrt(pmax(Eps, o.diag)/diag(m))
        m[] <- D * m * rep(D, each = n)
    }
    m
}

posdef.bock <- function(Sy, Se, large.small.ev.ratio=NULL) {
	# Create positive-definite covariance matrix according to the method
	# of Bock and Petersen, Biometrika 62(3):673-678 (1975).
	#
	# Assume that St = Sy - Se is an unbiased estimate for the covariance
	# matrix S, but is not necessarily positive-definite.  Goal is to
	# compute a corrected St.corr such that it is the maximum-likelihood
	# estimate under the restriction to positive-definite matrices and
	# assumed multivariate normality of the covariates.
	#
	# Solve generalized eigenvalue problem (Sy - lambda_i Se) x_i = 0
	# Approach:  Factorize Se = t(U) %*% U
	# Then with x = U^-1 z, we have
	# (U^-T Sy U^-1) z_i = lambda_i z_i or C Z = Z L
	# so by finding the eigendecomposition of C =(U^-T Sy U^-1) = Z, L we
	# can identify L as the eigenvalues and X = U^-1 Z as the eigenvectors
	# of the original problem.
	#
	# Then create L.star such that L.star_i = max(L_i, 1), and with B = X^-1,
	# return corrected matrix S = t(B) %*% (L.star - I) %*% B.

	# First, decompose Se into t(U) %*% U
	eSe <- eigen(Se)
	U <- diag(sqrt(eSe$values)) %*% t(eSe$vectors)
	#print(Se)
	#print(eSe)
	if (sum(abs(Se-t(U)%*%U))>sqrt(.Machine$double.eps)) {
		print(sum(abs(Se-t(U)%*%U)))
		stop("Se not symmetric.")
	}

	# Now compute C
	U.inv <- solve(U)
	C = t(U.inv) %*% Sy %*% U.inv
	# Solve eigenvalue problem
	eC = eigen(C)
	L = diag(eC$values)
	X = U.inv %*% eC$vectors

	# Now perform Bock and Petersen correction
	B = solve(X)
    if (is.null(large.small.ev.ratio)){
      min.ev <- 1 #+ .Machine$double.eps
    }
    else{
      min.ev <- 1 + eC$values[1]/large.small.ev.ratio
    }
	L.star = diag(pmax(eC$values,min.ev))
	I = diag(rep(1,ncol(L.star)))
	St.corr = t(B) %*% (L.star - I) %*% B
	#S.unc = t(B) %*% (L - I) %*% B  # the uncorrected result
	St.corr
}

# Returns the probability of a vector having a dot-product
# larger than dotprod.  (Scaled to unit vector!)
pvec <- function(vraw,dotprod,nreps=1000) {
	v <- as.vector(vraw)
	normv <- v/sqrt(sum(v^2)) # Scale
	l <- length(v)
	q <- sapply(1:nreps, function(n,y) {unitvec(l) %*% y}, y=normv)
	length(q[q>dotprod])/(1.0*nreps)
}

pcor <- function(c, prefix="") {
	cat(prefix, c$data.name, "r =", format(c$estimate,digits=3), "P =", format(c$p.value,digits=3), "\n", sep=" ")
}

pstat <- function(c, prefix="") {
	cat(prefix, c$data.name, "r =", format(c$statistic,digits=3), "P =", format(c$p.value,digits=3), "\n", sep=" ")
}

cortest <- function(x,y, meth="spearman", exact=FALSE, ...) {
	c <- cor.test(x, y, method=meth, exact=exact, ...)
	c$data.name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
	c
}

corr.list <- function(x, v1, vs, fxn=cortest, ret="estimate", ...) {
	y <- sapply(vs, function(m){
		fxn(x[,v1], x[,m],...)[[ret]]
	})
	y <- t(as.matrix(y))
	dimnames(y) <- list(v1,vs)
	return(y)
}

corr.my <- function(v1, vs, fxn=cortest, ret="estimate", ...) {
	y <- sapply(names(vs), function(m){fxn(v1,vs[[m]],...)[[ret]]})
	return(as.vector(y))
}

comb.comp <- function(g, indices) {
	x <- as.matrix((g$loadings[,indices] %*% as.matrix(g$eig$values[indices]))/sum(g$eig$values[indices]))
	rownames(x) <- rownames(g$loadings)
	x
}

entropy.dist <- function(p, base=2) {
	w <- getOption("warn")
	options(warn=-1)
	e <- sum(-p*log(p,base), na.rm=T)
	options(warn=w)
	e
}

conditional.entropy.dist <- function(pxy, base=2, ...) {
	# Compute conditional entropy H(X|Y) = H(X,Y) - H(Y)
	py <- rowSums(pxy)
	return(joint.entropy.dist(pxy,base,...)-entropy.dist(py,base,...))
}

joint.entropy.dist <- function(pxy,base=2) {
	# Computes H(X,Y) = sum_{x in X} sum_{y in Y} -p(x,y) log(p(y|x))
	entropy.sum <- 0
	w <- getOption("warn")
	options(warn=-1)
	entropy.sum <- sum(-pxy*log(pxy,base), na.rm=T)
	options(warn=w)
	entropy.sum
}

congruence <- function(f1raw, f2raw) {
	f1 <- as.vector(f1raw)
	f2 <- as.vector(f2raw)
	(f1 %*% f2 / sqrt(sum(f1^2)*sum(f2^2)))[1]
}

congruence.test <- function(f1raw, f2raw, n.trials=1000) {
	f1 <- as.vector(f1raw)
	f2 <- as.vector(f2raw)

	est <- congruence(f1, f2)
	pv <- pvec(f1, est, n.trials)
	cong <- list(est, estimate=est, n.obs=length(f1), p.value=pv)
	cong
}

down.congruence <- function(x1, x2) {
	flds <- intersect(names(x1), names(x2))
	congruence(x1[flds], x2[flds])
}

plot.cong <- function (x, x.axisnames=NULL, y.axisnames=NULL, colors=colorRampPalette(c("blue","black","yellow"),space="rgb")(50), zlim=c(-1,1), ...) {
	x <- as.matrix(x)
    nc <- ncol(x)
	nr <- nrow(x)
    corr <- x
    image(1:nc, 1:nr, t(corr)[1:nc,nr:1], col = colors,
        axes = FALSE, xlab='', ylab='', zlim=zlim, ...)

    if (!is.null(x.axisnames)) {
      if (any(is.na(x.axisnames))) {
        xaxn <- colnames(x)
      }
      else {
        xaxn <- x.axisnames
      }
      axis(3, at = 1:nc, labels = xaxn, las=2, ...)
    }
    if (!is.null(y.axisnames)) {
      if (any(is.na(y.axisnames))) {
        yaxn <- rev(rownames(x))
      }
      else {
        yaxn <- rev(y.axisnames)
      }
      axis(2, at = 1:nr, labels = yaxn, las=1, ...)
    }
    box()
}


# Congruence test for matrices, ala dist
cong.dist <- function(x, cong.fxn=congruence) {
	#res <- matrix(ncol=nrow(x), nrow=nrow(x), dimnames=list(rownames(x),rownames(y)))
	res <- vector(mode="numeric", length=(nrow(x)*(nrow(x)-1))/2)
	n = 1
	for (i in 1:(nrow(x)-1)) {
		for (j in (i+1):nrow(x)) {
			res[n] <- 1-cong.fxn(x[i,], x[j,])^2
			#cat("e: ",j + i*nrow(x),"\n")
			n <- n+1
		}
	}
	#cat("l:",n,"\n")
	attr(res,"Labels") <- rownames(x)
	attr(res,"Size") <- nrow(x)
	class(res) <- "dist"
	return(res)
}

# Congruence test for matrices, ala rcorr
rcong <- function(x,y, cong.fxn=congruence) {
	res <- matrix(ncol=ncol(y), nrow=ncol(x), dimnames=list(colnames(x),colnames(y)))
	for (i in 1:ncol(x)) {
		for (j in 1:ncol(y)) {
			res[i,j] <- cong.fxn(x[,i], y[,j])
		}
	}
	res
}

# Compute the root-mean-square distance
rmsd <- function(x,y) {
	sqrt(mean((x-y)^2))
}

rmsd.dist <- function(x) {
	#sapply(
}

# Reconstruct a matrix from the eigenvectors and eigenvalues
# indicated by indices.
eigen.reconstruct <- function(eig, indices) {
	vecs <- eig$vectors[,indices]
	vals <- diag(eig$values[indices],length(indices))
	res <- vecs %*% vals %*% t(vecs)
	res
}

# This function allows a correlation/covariance matrix to be submitted
# directly for principal components analysis.
matrix.prcomp <- function(raw.cmat, select.flds=NULL, raw.data=NULL, scores=FALSE, meth="spearman") {
	#print(deparse(substitute(raw.cmat)))
	if (is.null(select.flds)) {
		if (is.null(colnames(raw.cmat))) {
			select.flds <- 1:ncol(raw.cmat)
		}
		else {
			select.flds <- colnames(raw.cmat)
		}
	}
	cmat <- raw.cmat[select.flds,select.flds]
	#princomp(x=raw.data, covmat=cmat, scores=scores)
	e <- eigen(cmat)
	loadings <- e$vectors
	# Sign convention: fewest negative signs.
	loading.signs <- sapply(1:ncol(loadings), function(m) {
		if (sign(loadings[1,m])<0) {
			return(-1)
		}
		return(1)
	})
	loadings <- loadings %*% diag(unlist(loading.signs))
	rownames(loadings) <- rownames(cmat)
	colnames(loadings) <- paste("Comp.",1:length(e$values),sep='')

	prop.var <- e$values/sum(e$values)
	names(prop.var) <- colnames(loadings)
	res <- list(eig=e, loadings=loadings, prop.var=prop.var)
	if (scores) {
		res$scores <- as.data.frame(as.matrix(raw.data[,select.flds]) %*% as.matrix(e$vectors))
		names(res$scores) <- colnames(loadings)
		# This is inefficient, as we're computing many more correlations than are needed.
		df <- data.frame(raw.data[,select.flds], res$scores)
		res$r <- rcormat(df)$r[select.flds, colnames(loadings)]
		class(res$r) <- "matrix.prcomp.r.squared"
		res$r.squared <- res$r^2
		class(res$r.squared) <- "matrix.prcomp.r.squared"
	}
	class(res) <- "matrix.prcomp"
	res
}

print.matrix.prcomp.r.squared <- function(mpcr2) {
	print(round(unclass(mpcr2),2))
}

print.matrix.prcomp <- function(mpc) {
	cat("Loadings:\n")
	print(round(mpc$loadings,2))
	cat("\nProportion of variance explained:\n")
	print(round(mpc$prop.var, 3))
	mpc
}

# Reconstruct a data matrix from the data and the principal components
reconstruct.data <- function(pc, components, raw.data) {
	mat <- pc$loadings[,components]
	recon.data <- as.matrix(raw.data) %*% (mat %*% t(mat))
	recon.data
}

# Compute the squared deviation of the reconstructed data from the original data.
reconstruction.error <- function(raw.data, recon.data, vars=NULL) {
	if (is.null(vars)) {
		res <- (raw.data - recon.data)^2
	}
	else {
		res <- (raw.data[,vars]-recon.data[,vars])^2
	}
	res
}

# Compute the mean squared error (MSE) in reconstructing the i'th variable
# using the j'th eigenvalue/eigenvector.  Use analytical formula.
mse.prcomp <- function(pc,data=NULL) {
	res <- pc$loadings
	if (is.null(data)) {
		for (i in 1:nrow(res)) {
			for (j in 1:ncol(res)) {
				res[i,j] <- (1-pc$eig$values[j]*pc$loadings[i,j]^2)
			}
		}
	}
	else {
		flds <- rownames(pc$loadings)
		d <- as.matrix(data[,flds])
		n <- mean(apply(!is.na(d),2,sum))
		for (i in 1:nrow(res)) {
			for (j in 1:ncol(res)) {
				res[i,j] <- ((n-1)/n)*(1-pc$eig$values[j]*pc$loadings[i,j]^2)
			}
		}
	}
	res
}

# Compute the mean squared error (MSE) in reconstructing the i'th variable
# using the j'th eigenvalue/eigenvector.
mse.prcomp.exact <- function(pc,raw.data) {
	res <- pc$loadings
	flds <- rownames(pc$loadings)
	scaled.raw.data <- scalerankcols(raw.data[,flds])

	for (j in 1:ncol(res)) {
		d <- reconstruct.data(pc,j,scaled.raw.data)
		for (i in 1:nrow(res)) {
			res[i,j] <- mean((d[,i]-scaled.raw.data[,i])^2,na.rm=T)
		}
	}
	res
}


combine.prcomp <- function(pc, inds) {
	a <- pc$eig$values
	v <- pc$eig$vectors
	res <- v[,inds] %*% matrix(a[inds]/sum(a[inds]),ncol=1,nrow=length(inds))
	rownames(res) <- rownames(v)
	res
}

# Plot means of equi-distantly spaced bins
ez.plotmeans <- function(x,y, n.bins=20, xlimit=NULL, log.x=FALSE, ...) {
    d <- data.frame(x,y)
    if (is.null(xlimit)) {
        xlimit <- c(min(x), max(x))
		#print(xlimit)
    }
    breaks <- seq(xlimit[1], xlimit[2], (xlimit[2]-xlimit[1])/n.bins)
    if (log.x) {
        breaks <- 10^seq(log10(xlimit[1]),log10(xlimit[2]),log10(xlimit[2]/xlimit[1])/n.bins)
		#breaks <- axTicks(1, axp=c(xlimit[1], xlimit[2], 1), log=T)
		#print(breaks)
    }
    d$cat <- cut(d$x, breaks=breaks)
    plotmeans(y~cat, data=d, n.label=F, xaxt='n', ...)
    axis(1, at=1:length(breaks), labels=breaks)

}

# Get function of data in equi-distantly spaced bins
ez.bin <- function(x,y, fxn=mean, n.bins=20, xlimit=NULL, log.x=FALSE, log.y=FALSE, ...) {
    d <- na.omit(data.frame(x,y))
    if (is.null(xlimit)) {
        xlimit <- c(min(d$x), max(d$x))
		#print(xlimit)
    }
    breaks <- seq(xlimit[1], xlimit[2], (xlimit[2]-xlimit[1])/(n.bins))
    #print(length(breaks))
    #print(breaks)
    if (log.x) {
        breaks <- 10^seq(log10(xlimit[1]),log10(xlimit[2]),log10(xlimit[2]/xlimit[1])/(n.bins))
		#breaks <- axTicks(1, axp=c(xlimit[1], xlimit[2], 1), log=T)
		#print(breaks)
    }
    d$cat <- cut(d$x, breaks=breaks, include.lowest=TRUE)
    if (log.y) {
	    fxn.out <- as.vector(by(d$y, d$cat, fxn, na.rm=TRUE))
	}
	else {
	    fxn.out <- 10^(as.vector(by(log10(d$y), d$cat, fxn, na.rm=TRUE)))
	}
    #breaks.plus.one <- append(c(0),breaks[1:length(breaks)-1])
    #print(length(breaks.plus.one))
    if (log.x) {
	    breaks.mids <- 10^((log10(breaks[2:length(breaks)])+log10(breaks[1:length(breaks)-1]))/2)
	}
	else {
	    breaks.mids <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
	}
	#print(breaks.mids)
    d <- data.frame(breaks.mids, fxn.out)
    colnames(d) <- c("bin.mid",deparse(substitute(fxn)))
    d
}

# Get function of data in equi-distantly spaced bins
ez.equal.bin <- function(x,y, fxn=mean, n.bins=20, log.x=FALSE, log.y=FALSE, ...) {
    d <- na.omit(data.frame(x,y))
    breaks <- quantile(d$x, probs=seq(0,1,1/(n.bins))) #[1:n.bins]
    u.breaks <- unique(breaks)
    #print(length(breaks))
    #print(breaks)
    #print(unique(breaks))
    if (length(u.breaks) < length(breaks)) {
    	d$cat <- cut(rank(d$x, ties.method='random'), breaks=seq(0,length(d$x),length(d$x)/n.bins), include.lowest=FALSE)
    }
    else {
	    d$cat <- cut(d$x, breaks=breaks, include.lowest=TRUE)
	}
    if (log.y) {
	    fxn.out <- as.vector(by(d$y, d$cat, fxn, na.rm=TRUE))
	}
	else {
	    fxn.out <- 10^(as.vector(by(log10(d$y), d$cat, fxn, na.rm=TRUE)))
	}
    #breaks.plus.one <- append(c(0),breaks[1:length(breaks)-1])
    #print(length(breaks.plus.one))
    if (log.x) {
	    breaks.mids <- 10^((log10(breaks[2:length(breaks)])+log10(breaks[1:length(breaks)-1]))/2)
	}
	else {
	    breaks.mids <- (breaks[2:length(breaks)]+breaks[1:length(breaks)-1])/2
	}
	#print(breaks.mids)
    d <- data.frame(breaks.mids, fxn.out)
    colnames(d) <- c("bin.mid",deparse(substitute(fxn)))
    d
}


# Get function of data in windows
ez.window <- function(x,y, n.per.window, step.size, x.fxn=median, y.fxn=mean, log.x=FALSE, log.y=FALSE, ...) {
    d <- na.omit(data.frame(x,y))
    d <- d[order(d$x),]
    m <- sapply(seq(n.per.window, nrow(d), step.size), function(window.max) {
    	n <- nrow(d)
    	if (log.x) {
	    	t.x <- log10(d$x[max(1,(window.max-n.per.window)):min(n,window.max)])
	    	x.point <- 10^(x.fxn(t.x))
	    } else {
	    	t.x <- d$x[max(1,(window.max-n.per.window)):min(n,window.max)]
	    	x.point <- x.fxn(t.x)
	    }
	    if (log.y) {
	    	t.y <- log10(d$y[max(1,(window.max-n.per.window)):min(n,window.max)])
	    	y.point <- 10^(y.fxn(t.y))
	    } else {
	    	t.y <- d$y[max(1,(window.max-n.per.window)):min(n,window.max)]
	    	y.point <- y.fxn(t.y)
	    }
    	c(x.point, y.point)
    	})
    #print(m)
    d <- data.frame(t(m))
    colnames(d) <- c(paste("bin.",deparse(substitute(x.fxn)),sep=''),deparse(substitute(y.fxn)))
    d
}

noop <- function(x) {
	x
}

multi.ecdf <- function(x, log=F, col=NULL, lty="solid", lwd=1, legend.at=NULL, xlim=NULL, ylim=NULL,
	equal.height=F, relative.heights=NULL, xlab="x", ylab="Empirical CDF", weight.list=NULL, ...) {
	extra.args <- list(...)
	if (is.data.frame(x) || is.matrix(x)) {
		col.names <- colnames(x)
		x <- lapply(1:ncol(x),function(m){x[,m]})
		names(x) <- col.names
	}
	if (!is.list(x)) {
		x <- list(x)
	}
	if (is.null(col)) {col <- rainbow(length(x))}
	## Extend properties into vectors, if necessary
	cols <- as.vector(replicate(length(x)/length(col) + 1,col))
	ltys <- as.vector(replicate(length(x)/length(c(lty))+1,lty))
	lwds <- as.vector(replicate(length(x)/length(c(lwd))+1,lwd))
	#cat("h2\n")
	if (log) {
		trans <- log10
		inv.trans <- function(y) {10^y}
	}
	else {
		trans <- noop
		inv.trans <- noop
	}

	densities <- lapply(1:length(x), function(n) {
      y <- x[[n]]
      ny <- na.omit(y)
      d <- ecdf(trans(ny))
      d
    })
	#cat("h4\n")

	## X limits
	if (is.null(xlim)) {
		## Make xlims
		valid.x <- x
		if (log) {
			valid.x <- lapply(x, function(y) {y[y>0]})
		}
		xmin <- min(sapply(valid.x,min,na.rm=T),na.rm=T)
		xlim <- c(xmin, max(sapply(x,max,na.rm=T),na.rm=T))
	}

    length.out <- 1000
    at <- seq(trans(xlim[1]), trans(xlim[2]), length.out=length.out)

	if (log) {log.str <- "x"} else {log.str <- ""}
	plot(inv.trans(at), densities[[1]](at), type='n', col=col[1], xlim=xlim, ylim=ylim, lty=lty, lwd=lwd, log=log.str, xlab=xlab, ylab=ylab, ...)
	for (i in 1:length(x)) {
		d <- densities[[i]]
		lines(inv.trans(at), d(at), col=cols[i], lty=ltys[i], lwd=lwds[i], ylim=ylim, ...)
	}

	if (!is.null(legend.at)) {
		legend.names = names(x)
		if (is.null(legend.names)) {
		  legend.names = as.character(1:length(x))
		}
		legend.cols <- col[1:min(length(x), length(col))]
		legend(legend.at[1], legend.at[2], col=legend.cols, legend=legend.names, lty=ltys)
	}
}

## Takes a list of variables, plots kernel densities
multi.density <- function(x, log=F, kernel="r", col=NULL, lty="solid", fill=FALSE, lwd=1, legend.at=NULL, xlim=NULL, ylim=NULL,
	equal.height=F, relative.heights=NULL, xlab="x", ylab="Density", weight.list=NULL, cex.legend=1, ...) {
	extra.args <- list(...)
	if (is.data.frame(x) || is.matrix(x)) {
		col.names <- colnames(x)
		x <- lapply(1:ncol(x),function(m){x[,m]})
		names(x) <- col.names
	}
	if (!is.list(x)) {
		x <- list(x)
	}
	#cat("h1\n")
	## Colors
	if (is.null(col)) {col <- rainbow(length(x))}
	## Extend properties into vectors, if necessary
	cols <- as.vector(replicate(length(x)/length(col) + 1,col))
	ltys <- as.vector(replicate(length(x)/length(c(lty))+1,lty))
	lwds <- as.vector(replicate(length(x)/length(c(lwd))+1,lwd))
	#cat("h2\n")
	if (log) {
		trans <- log10
		inv.trans <- function(y) {10^y}
	}
	else {
		trans <- noop
		inv.trans <- noop
	}
	#cat("h3\n")
	## Compute the densities
	densities <- lapply(1:length(x), function(n) {
		if (is.null(weight.list)) {
			y <- x[[n]]
			ny <- na.omit(y)
			d <- density(trans(ny), na.rm=T, kern=kernel)
		}
		else {
			y <- data.frame(x=x[[n]], w=weight.list[[n]])
			ny <- na.omit(y)
			d <- density(trans(ny$x), na.rm=T, kern=kernel, weights=ny$w/sum(ny$w,na.rm=T))
		}
		d
		})
	#cat("h4\n")
	max.heights <- lapply(densities, function(d) {max(d$y, na.rm=T)})
	max.height <- max(unlist(max.heights), na.rm=T)
	if (!is.null(relative.heights)) {
		relative.heights <- relative.heights/max(relative.heights, na.rm=T)
		equal.height=T
	}
	else {
		relative.heights <- seq(1,1,length.out=length(x))
	}
	if (is.null(ylim)) {
		if (equal.height) {
			ylim=c(0,1.05)
		}
		else {
			ylim=c(0,1.05*max.height)
		}
	}

	## X limits
	if (is.null(xlim)) {
		## Make xlims
		valid.x <- x
		if (log) {
			valid.x <- lapply(x, function(y) {y[y>0]})
		}
		xmin <- min(sapply(valid.x,min,na.rm=T),na.rm=T)
		xlim <- c(xmin, max(sapply(x,max,na.rm=T),na.rm=T))
	}

	## Plot the first dataset
	if (log) {log.str <- "x"} else {log.str <- ""}
	if (equal.height) {height.div <- max.height} else {height.div <- 1.0}
	## Infer X and Y labels, if not passed in
	if (missing(xlab) | is.null(xlab)) {
		if (!is.null(names(x)[1])) {
			xlab <- names(x)[1]
		}
		else {
			xlab <- "x"
		}
	}
	if (missing(ylab) | is.null(ylab)) {
		if (equal.height) {
			ylab <- "Peak-normalized density"
		}
		else {
			ylab <- "Density"
		}
	}

	## Actually plot the data
	d <- densities[[1]]
	plot(inv.trans(d$x), d$y/height.div, type='n', col=col[1], xlim=xlim, ylim=ylim, lty=lty, lwd=lwd, log=log.str, xlab=xlab, ylab=ylab, ...)
	for (i in 1:length(x)) {
		d <- densities[[i]]
		height.div <- 1.0
		if (equal.height) {
			height.div <- max.heights[[i]]
		}
		if (fill) {
			polygon(inv.trans(d$x), relative.heights[[i]]*d$y/height.div, col=cols[i], lty=ltys[i], lwd=lwds[i], ylim=ylim, ...)
		}
		else {
			lines(inv.trans(d$x), relative.heights[[i]]*d$y/height.div, col=cols[i], lty=ltys[i], lwd=lwds[i], ylim=ylim, ...)
		}
	}

	## cat("h6\n")
	## Legend
	if (!is.null(legend.at)) {
		legend.names = names(x)
		if (is.null(legend.names)) {
		  legend.names = as.character(1:length(x))
		}
		legend.cols <- col[1:min(length(x), length(col))]
		legend(legend.at[1], legend.at[2], col=legend.cols, legend=legend.names, lty=ltys, cex=cex.legend)
	}
}

multidens <- multi.density

multi.lm <- function(response, predictors, data, rank=F, na.last=T){
	if (rank) {
		data <- rankcols(data[append(response,predictors)], na.last=na.last)
	}
	adds <- paste(predictors,collapse="+",sep="")
	f = paste(response,"~",adds,sep="")
	lm(formula(f),data=data)
}
