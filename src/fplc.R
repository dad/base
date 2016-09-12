fplc.sampling <- function(x, from.ml=NULL, to.ml=NULL, ml.fld='ml', na.rm=TRUE, ...) {
	# Determine sampling width for FPLC data
	xs <- x
	if (!is.null(from.ml)) {
		xs <- tsubset(xs, paste(ml.fld, ">=", from.ml))
	}
	if (!is.null(to.ml)) {
		xs <- tsubset(xs, paste(ml.fld, "<=", to.ml))
	}
	mld <- xs[,ml.fld]
	n <- length(mld)
	widths <- mld[2:n]-mld[1:(n-1)]
	mean(widths, na.rm=na.rm, ...)
}

find.break <- function(y, val, from.below=T) {
	# Find mL points where quant is above and below val 
	# Assume y is sorted by mL
	# First check to see whether val is inside available values -- if not, there is no defined break
	if (val >= min(y$mAU, na.rm=T) & val <= max(y$mAU, na.rm=T)) {
		if (from.below) {
			below <- subset(y, mAU<val)
			lower.mL <- below[nrow(below),'ml']
			above <- subset(y, mAU>=val)
			upper.mL <- above[1,'ml']
			#res <- approx(c(below[nrow(below),'mL'], above[1,'mL']), c(below[nrow(below),'mAU'], above[1,'mAU']), xout=lower.mL + ())
		} else {		
			below <- subset(y, mAU>=val)
			lower.mL <- below[nrow(below),'ml']
			above <- subset(y, mAU<val)
			upper.mL <- above[1,'ml']
		}
		break.mL <- mean(lower.mL, upper.mL)
		res <- list(at=break.mL, quant=fplc.sample.ml(y, break.mL))

	} else {
		# val not in available values. 
		if (from.below)
			res <- list(at=y[1,'ml'], quant=y[1,'mAU'])
		else
			res <- list(at=y[nrow(y),'ml'], quant=y[nrow(y),'mAU'])
	}
	res
}

peakdim.full <- function(x, from.ml, to.ml, prop=0.5, quant='mAU') {
	# Find full width at half-maximum (or prop-maximum)
	res <- fplc.max(x, from.ml, to.ml, quant=quant)
	prop.max = res$max*prop
	# Find mL points where quant is above and below half-max 
	lower.sub <- subset(x, ml>from.ml & ml<=res$at)
	if (nrow(lower.sub)>0){
		lower.break <- find.break(lower.sub, prop.max, from.below=TRUE)
	} else {
		lower.break <- list(at=from.ml)
	}
	upper.sub <- subset(x, ml>res$at & ml<=to.ml)
	if (nrow(upper.sub)>0){
		upper.break <- find.break(upper.sub, prop.max, from.below=FALSE)
	} else {
		upper.break <- list(at=to.ml)
	}

	list(max=res$max, lower=lower.break$at, upper=upper.break$at, width=upper.break$at-lower.break$at)
}
peakdim <- function(peak, x, prop=0.5, quant='mAU') {
	peakdim.full(x, peak$lower, peak$upper, prop=prop, quant=quant)
}

fullwidth.halfmax.full <- function(x, from.ml, to.ml, prop=0.5, quant='mAU') {
	peakdim.full(x, from.ml, to.ml, prop=prop, quant=quant)
}

area.under.curve.full <- function(x, from.ml, to.ml, baseline=0, quant='mAU', na.rm=TRUE) {
	# Quantify FPLC data
	# Determine mean sampling bin size in mL
	binwidth <- fplc.sampling(x, na.rm=na.rm)
	xs <- subset(x, ml>from.ml & ml<=to.ml)
	nbins <- nrow(xs)
	vol <- sum(xs[,quant]-baseline, na.rm=na.rm)*binwidth
	vol
}
area.under.curve <- function(peak, x, ...) {
	area.under.curve.full(x, peak$lower, peak$upper, ...)
}

fwhm.gaussian <- function(x, from.ml, to.ml, quant='mAU') {
	# Find a normal distribution with the same FWHM
	gmax <- fplc.max(x, from.ml, to.ml, quant)
	fwhm <- fullwidth.halfmax.full(x, from.ml, to.ml, prop=0.5, quant=quant)
	sigma <- fwhm$width/(2*sqrt(2*log(2)))
	mult <- gmax$max/dnorm(gmax$at, mean=gmax$at, sd=sigma)
	list(mean=gmax$at, sd=sigma, mult=mult)
}

auc.fwhm.gaussian <- function(x, from.ml, to.ml, quant='mAU') {
	# Area under the curve, assuming Gaussian
}

fplc.stat <- function(x, from.ml, to.ml, quant='mAU', stat=median, ...) {
	# Quantify FPLC data
	# Determine mean sampling bin size in mL
	binwidth <- fplc.sampling(x)
	xs <- subset(x, ml>=from.ml & ml<=to.ml)
	stat(xs[,quant], ...)
}

fplc.at.max <- function(x, at=NULL, from.ml, to.ml, quant='mAU') {
	# Quantify FPLC data
	# Determine maximum value
	y <- subset(x, ml>from.ml & ml<=to.ml)
	y$score <- y[,quant]
	mean(subset(y, score==max(y$score))[,quant])
}

fplc.max <- function(x, from.ml, to.ml, quant='mAU') {
	# Quantify FPLC data
	# Determine maximum value and position
	y <- subset(x, ml>from.ml & ml<=to.ml)
	y$score <- y[,quant]
	max.score <- max(y$score, na.rm=T)
	list(max=max.score, at=mean(subset(y, score==max.score)[,'ml']))
}

fplc.sample.ml <- function(x,sample.ml) {
	if (is.data.frame(x)) {x <- list(x)}
	sapply(x, function(xs){
		approx(xs$ml, xs$mAU, xout=sample.ml)$y
		})
}

mlrange <- function(lower,upper) {
	list(lower=lower, upper=upper)
}


# DAD: bandwidth bw based on indices makes little sense much of the time.
# Move to support bw based on data
find.extrema <- function(x, y=NULL, bw=20, method=c('max','min'), limit.value=0.0, ...) {
	# Ensure right lengths
	if (is.list(x)) {
		y <- x[[2]]
		x <- x[[1]]
	} else {
		if (is.null(y)) {
			y <- x # Treat vector as array of heights
			x <- 1:length(y) # Make list of integers for positions
		}
		stopifnot(length(y) == length(x))
	}
	n <- length(y)

	# Is extremum? (min or max)
	# DAD: could be made much more efficient, but for FPLC data, no need.
	method <- match.arg(method)
	if (method == 'max') {
		fn <- max
	}
	if (method == 'min') {
		fn <- min
	}
	is.extremum <- function(i, x, bw=bw, fn=fn, ...) {
		!is.na(x[i]) & (x[i] == fn(x[(i-bw):(i+bw)], ...))
	}

	# Move across data, identifying local extrema within window of +/- bw
	maxness <- c(rep(NA,bw),sapply((bw+1):(n-bw), is.extremum, x=y, bw=bw, fn=fn, ...),rep(NA,bw))
	maxi <- na.omit((1:length(x))[maxness])

	res <- data.frame(x=x[maxi], y=y[maxi], i=maxi)
	res
}

# Subset of data x that falls within peak [lower,upper) interval
peak.subset <- function(x, peaks) {
	if (class(peaks)=='fplc.peak') {
		peaks <- list(peaks)
	}
	res <- lapply(peaks, function(p){subset(x, ml>=p$lower & ml<p$upper)})
	if (length(res)==1) {
		return(res[[1]])
	}
	res
}

find.closest <- function(x, start.x, target.y, vx='ml', vy='mAU') {
	# Given a range (start.x, target.x), (start.y, target.y)
	# Find first value in y that exceeds limit.y

	res <- list(lower=NA, upper=NA)
	x$.ind <- 1:nrow(x)
	x$.x <- x[,vx]
	x$.target <- x[,vy]
	x$.dist <- (start.x - x$.x)^2
	# Subsets outside the acceptable range: we want the closest (minimum-distance)
	# values from these subsets
	x.sub.lower <- subset(x, .target<target.y & .x<start.x)
	x.sub.upper <- subset(x, .target<target.y & .x>start.x)

	# If there are no fitting results
	if (nrow(x.sub.lower)==0) {
		res$lower <- x[1,vx]
	} else {
		close <- min(x.sub.lower$.dist, na.rm=T)
		res$lower <- x[subset(x.sub.lower, .dist==close)$.ind,vx]
	}
	if (nrow(x.sub.upper)==0) {
		res$upper <- x[nrow(x),vx]
	} else {
		close <- min(x.sub.upper$.dist, na.rm=T)
		res$upper <- x[subset(x.sub.upper, .dist==close)$.ind,vx]
	}

	res
}

filter.peak.list <- function(plist, x, min.peak.dist=NULL, min.base=NULL) {
	# Filter NULL values
	pl <- plist[!sapply(plist, is.null)]
	# Filter based on peak distance
	filt.pl <- pl
	if (!is.null(min.peak.dist) & length(pl)>1) {
		if (min.peak.dist>0) {
			bad.peaks <- NULL # indices of bad peaks
			for (xi in 1:(length(pl)-1)) {
				px <- pl[[xi]]
				for (yi in (xi+1):length(pl)) {
					if (!(yi %in% bad.peaks)) {
						py <- pl[[yi]]
						# If peaks are too close together...
						if (abs(px$x-py$x)<min.peak.dist) {
							# Keep only the tallest peak
							if (px$y>py$y) bad.peaks <- c(bad.peaks,yi) else bad.peaks <- c(bad.peaks,xi)
						}
					}
				}
			}
			filt.pl <- pl[!(1:length(pl) %in% bad.peaks)]
		}
	}
	if (!is.null(min.base)) {
		new.pl <- lapply(filt.pl, function(p){
			psub <- peak.subset(x,p)
			cl <- find.closest(psub, p$x, target.y=min.base)
			p$lower <- cl$lower
			p$upper <- cl$upper
			p
		})
		filt.pl <- new.pl
	}
	filt.pl
}

identify.peaks <- function(x, y=NULL, bw=20, min.peak=0, min.peak.dist=NULL, min.base=0, ...) {
	# Ensure right lengths
	x.orig <- x
	if (is.list(x)) {
		y <- x[[2]]
		x <- x[[1]]
	} else {
		if (is.null(y)) {
			y <- x # Treat vector as array of heights
			x <- 1:length(y) # Make list of integers for positions
		}
		stopifnot(length(y) == length(x))
	}
	maxs <- find.extrema(x, y, bw, limit.value=min.peak, method='max', ...)
	mins <- find.extrema(x, y, bw, limit.value=min.base, method='min', ...)
	# Put mins on bottom and top of range as sentinels
	padded.minx <- c(min(x,na.rm=T),mins$x,max(x,na.rm=T))
	#print(padded.minx)
	# For each maximum, find two nearest minima
	inds <- 1:length(padded.minx)
	peaklist <- lapply(1:nrow(maxs), function(m){
		maxpos <- maxs$x[m]
		maxy <- maxs$y[m]
		peak <- NULL
		# If this is a peak we may want to keep...
		if (maxy>min.peak) {
			peak <- list(x=maxpos, y=maxy, min=min.base)
			# Find closest mins on either side
			# First find closest min value
			q <- (padded.minx - maxpos)^2

			mini <- inds[q==min(q,na.rm=T)][[1]]
			minval <- padded.minx[mini]
			# Then find the min on the other side
			if (maxpos>minval){
				# Also grab
				peak$lower <- minval
				peak$upper <- padded.minx[mini+1]
			} else {
				peak$upper <- minval
				peak$lower <- padded.minx[mini-1]
			}
			class(peak) <- 'fplc.peak'
		}
		peak
	})
	filter.peak.list(peaklist, x.orig, min.peak.dist=min.peak.dist, min.base=min.base)
}

print.fplc.peak <- function(x) {
	cat('mL = ', x$x, "\nmax mAU = ", x$y, '\nlower = ', x$lower, '\nupper = ', x$upper, '\nmin = ', x$min, '\n', sep='')
}

quantify.peaks <- function(peaks, x, y=NULL, apply=FALSE, method=c('minmin', 'gaussian')) {
	n <- length(peaks)
	fwhm <- lapply(peaks, function(p){peakdim(p, x, prop=0.5)})
	#print(fwhm)
	restable <- data.frame(
		# Sequential identifier
		id = 1:n,
		# Elution volume
		ml = sapply(peaks, function(p){p$x}),
		# Peak maximum
		max = sapply(peaks, function(p){p$y}),
		# Area under the curve
		auc = sapply(peaks, function(p){area.under.curve(p, x)}),
		# Full width at half maximum
		fwhm = sapply(fwhm, function(f){f$width})
		)
	restable
}

plot.peaks <- function(x, peaks=NULL, annotate=TRUE, min.peak=0, min.base=0, xlab='Volume (mL)', ylab='UV absorbance (mAU)', las=1, ...) {
	plot(x[[1]], x[[2]], xlab=xlab, ylab=ylab, type='l', las=las, ...)
	plist <- peaks
	if (annotate) {
		if (is.null(plist)) {
			plist <- identify.peaks(x, min.peak=min.peak, min.base=min.base)
		}	
		abline(v=sapply(plist, function(x){x$x}));
	}
	invisible(plist)
}

extract.peak <- function(x, peak) {
	subset(x, ml>=peak$lower & ml<peak$upper)
}

peak.polygon <- function(peaks, x, y=NULL, col=NULL, min=0, ...) {
	if (!is.list(peaks)) {
		peaks <- list(peaks)
	}
	n <- length(peaks)
	if (is.null(col)){
		col <- myrainbow(n)
	} else {
		if (length(col)<n) col <- rep(col,n)
	}
	for (xi in enum(peaks)){
		peak <- xi$x
		peaksub <- extract.peak(x, peak)
		min.base <- if (!is.null(peak$min)) peak$min else 0
		flatpolygon(peaksub[[1]], peaksub[[2]], min=min.base, col=col[xi$i], ...)
	}
}

peakpos <- function(x) {
	if (class(x)=='fplc.peak') {
		x <- list(x)
	}
	sapply(x, function(y){y$x})
}

peakheight <- function(x) {
	if (class(x)=='fplc.peak') {
		x <- list(x)
	}
	sapply(x, function(y){y$y})
}