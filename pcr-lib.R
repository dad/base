library(pls)

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
    r <- cor(d, meth=method)
    g.cov <- pcr.covmat(form, covmat=r, n=nrow(d))
    cor.X <- r[pred.f,pred.f]
    cor.Xy <- r[pred.f,resp.f]
    e <- eigen(cor.X)
    u <- e$vectors
    b.svd <- coef(g.svd) #t(u) %*% solve(cor.X) %*% cor.Xy
    b.cov <- coef(g.cov)
    list(svd=b.svd, cov=b.cov)
  }
  sq.diff <- function(n, method) {
    x <- f(n, method)
    ((x$svd-x$cov)/x$svd)^2
  }
  method <- match.arg(method, c("pearson","spearman"))
  diffs <- replicate(10, sq.diff(1000, method))
  cat("Cumulative RMS % deviation between SVD and Cov methods (", method, ") = ", sqrt(mean(diffs)), "\n", sep='')
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

