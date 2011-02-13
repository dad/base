
label.gene.ratios <- function(x, target.genes, r='ratio.hl.normalized.mean', a=abund.f, cex.pt=1, cex.gene.lab=0.6, ...) {
  y <- subset(x, gene %in% target.genes)
  if (nrow(y)>0) {
  	points(y[,abund.f], y[,r], cex=cex.pt, ...)
  	text(y[,abund.f], y[,r], y$gene, pos=4, cex=cex.gene.lab, ...)
  }
}

plot.gene.ratios <- function(x, target.genes, ratio.f='ratio.hl.normalized', abund.f='ratio.hl.count', err=TRUE, labels=target.genes, ...) {
  y <- subset(x, gene %in% target.genes)
  if (nrow(y)>0) {
  	if (err) {
		points.err(y[,abund.f], y[,ratio.f], x.lower=y[,abund.f], x.upper=y[,abund.f], y.lower=y[,p.0(ratio.f,"lower.95")], y.upper=y[,p.0(ratio.f,"upper.95")], ...)
	} else {
		points(y[,abund.f], y[,ratio.f], ...)
	}
	text(y[,abund.f], y[,ratio.f], labels, pos=2)
  }
}
p.gene <- plot.gene.ratios

load.ms.data <- function(fname, invert=FALSE) {
	# Light = KG071, YFP-WT; Heavy = KG079, YFP-C4
	d.raw <- read.table(fname, header=T,stringsAsFactors=F)
	#desc <- read.delim("~/research/data/scerevisiae/scer-descs.txt", header=T)
	all.orfs <- unique(as.vector(c(d.raw$orf, yeast$bg$orf)))
	zd <- match(all.orfs, d.raw$orf)
	zy <- match(all.orfs, yeast$avg$orf)
	#z.desc <- match(d.raw$orf, desc$orf)
	d <- data.frame(d.raw[zd,], yeast$avg[zy,], yeast$raw[zy,yeast$fields$abund], yeast$bg[zy,c('gene','desc')])
	y <- yeast$raw
	intensity.fld <- 'intensity'
	abund.fld <- 'abund.1'
	d$ratio.significance.corrected <- p.adjust(d$ratio.significance, "fdr")
	d$intensity.proportion <- d[,intensity.fld]/sum(as.numeric(d[,intensity.fld]), na.rm=T)
	tot.abund <- colSums(na.omit(d[,c(abund.fld,intensity.fld)]))
	d$estimated.abundance <- d[,intensity.fld]*tot.abund[1]/tot.abund[2]
	d$est.abund.proportion <- d$estimated.abundance/sum(yeast$raw[,abund.fld],na.rm=T)
	d[is.na(d$gene),'gene'] <- d[is.na(d$gene),'orf']
	d <- d[order(d$est.abund.proportion, decreasing=T),]
	d$order <- 1:nrow(d)
	if (invert) {
		invert.flds <- c('ratio.hl','ratio.hl.normalized','ratio.hl.mean','ratio.hl.normalized.mean')
		d[,invert.flds] <- 1/d[,invert.flds]
	}
	d
}

sampling.error <- function(n, log.sd.cutoff) {
	log.sd.cutoff/sqrt(n)
}

label.silac <- function(x, abund.f='ratio.hl.count', ratio.f="ratio.hl.normalized.mean", lab.f="gene", ...) {
	text(x=x[,abund.f], y=x[,ratio.f], x[,lab.f], ...)
}

plot.silac.updown <- function(x, count.cutoff, up.orfs, down.orfs, lab.fld="gene", abund.f="ratio.hl.count", sd.envelope=0.995, cex.pt=1, cex.gene.lab=0.6, xlim=NULL, ...) {
	ms.sub <- subset(x, ratio.hl.count>=count.cutoff)
	sig.sub <- ms.sub[match(c(up.orfs,down.orfs),ms.sub$orf),]
	up.sig.sub <- ms.sub[match(up.orfs,ms.sub$orf),]
	down.sig.sub <- ms.sub[match(down.orfs,ms.sub$orf),]
	#sig.sub <- subset(ms.sub, orf %in% c(up.orfs,down.orfs))
	#up.sig.sub <- subset(ms.sub, orf %in% up.orfs)
	#down.sig.sub <- subset(ms.sub, orf %in% down.orfs)
	if (is.null(xlim)) {
		xlim <- c(min(x[,abund.f], na.rm=T)+1,max(x[,abund.f], na.rm=T)*1.1)
		print(xlim)
	}
	plot(ms.sub[,abund.f], ms.sub$ratio.hl.normalized.mean, type='n', log='xy', xaxt='n', pch=16, las=1, xlim, ...)
	abline(h=1, col='lightgray')
	# Error limits
	f <- function(n, log.sd.cutoff) {
		log.sd.cutoff/sqrt(n)
	}
	sd.logs <- subset(x, ratio.hl.count>2)$ratio.hl.normalized.sd
	sd.log.prop <- sd.envelope
	sd.log.cutoff <- sd.logs[order(sd.logs)][ceiling(length(sd.logs)*sd.log.prop)]
	env.lims <- exp(seq(log(xlim[1]/2),log(xlim[2]*2), length.out=100))
	lines(env.lims, exp(sampling.error(env.lims, sd.log.cutoff)), col='gray75', ...)
	lines(env.lims, exp(-sampling.error(env.lims, sd.log.cutoff)), col='gray75', ...)
	# Points
	points(ms.sub[,abund.f], ms.sub$ratio.hl.normalized.mean, pch=18, col="gray70", cex=0.5*cex.pt, ...)
	points(sig.sub[,abund.f], sig.sub$ratio.hl.normalized.mean, col='gray10', pch=16, cex=cex.pt, ...)
	segments(sig.sub[,abund.f], sig.sub$ratio.hl.normalized.lower.95, sig.sub[,abund.f], sig.sub$ratio.hl.normalized.upper.95, col='black', ...)
	d <- list(...)
	xlim <- d$xlim
	if (is.null(xlim)) {
		my.axis(1, ms.sub[,abund.f], log=T, las=1, ...)
	} else {
		my.axis(1, xlim, log=T, las=1, ...)
	}
	if (nrow(up.sig.sub) > 0 & !is.null(lab.fld)) {
		label.silac(up.sig.sub, pos=4, cex=cex.gene.lab, col="gray5")
		#text(x=up.sig.sub[,abund.f], y=up.sig.sub$ratio.hl.normalized.mean, up.sig.sub[,lab.fld], pos=4, cex=cex.gene.lab, col="gray5")
	}
	if (nrow(down.sig.sub) > 0 & !is.null(lab.fld)) {
		label.silac(down.sig.sub, pos=4, cex=cex.gene.lab, col="gray5")
		#text(x=down.sig.sub[,abund.f], y=down.sig.sub$ratio.hl.normalized.mean, labels=down.sig.sub[,lab.fld], pos=4, cex=cex.gene.lab, col="gray5")
	}
	list(all=ms.sub, up=up.sig.sub, down=down.sig.sub)
}

get.silac.updown <- function(x, count.cutoff=2, sd.envelope=0.95) {
	ms.sub <- subset(x, ratio.hl.count>=count.cutoff)
	sd.logs <- subset(x, ratio.hl.count>=max(2,count.cutoff))$ratio.hl.normalized.sd
	sd.log.prop <- sd.envelope
	sd.log.cutoff <- sd.logs[order(sd.logs)][ceiling(length(sd.logs)*sd.log.prop)]
	
	ms.up <- subset(ms.sub, ratio.hl.normalized.mean>1 & ratio.hl.normalized.lower.95>exp(sd.log.cutoff/sqrt(ratio.hl.count)))
	ms.down <- subset(ms.sub, ratio.hl.normalized.mean<1 & ratio.hl.normalized.upper.95<exp(-sd.log.cutoff/sqrt(ratio.hl.count)))
	list(up=ms.up[order(ms.up$ratio.hl.normalized.mean, decreasing=TRUE),], down=ms.down[order(ms.down$ratio.hl.normalized.mean, decreasing=FALSE),], sd.log.cutoff=sd.log.cutoff)
}

plot.silac <- function(x, diff.list=NULL, lab.fld="gene", abund.f="ratio.hl.count", cex.pt=1, cex.gene.lab=0.6, xlim=NULL, xlab="Ratio count", ylab='Protein ratio', las=1, ...) {
	ms.sub <- subset(x, ratio.hl.count>0)
	sig.sub <- rbind(diff.list$up, diff.list$down)
	# Limits
	d <- list(...)
	xlim <- d$xlim
	if (is.null(xlim)) {
		xlim <- c(1,max(ms.sub[,abund.f], na.rm=T)*1.1)
	}
	if (is.null(d$ylim)) {
		# Set y limits so that we can always see the full envelope, and also encompass all significantly changed proteins and their error bars.
		ylim <- c(min(exp(-diff.list$sd.log.cutoff),min(sig.sub$ratio.hl.normalized.lower.95,na.rm=T)),
			max(exp(diff.list$sd.log.cutoff),max(sig.sub$ratio.hl.normalized.upper.95,na.rm=T)))
	}
	# Initial plot for labels and limits
	plot(ms.sub[,abund.f], ms.sub$ratio.hl.normalized.mean, type='n', log='xy', xaxt='n', xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, las=las, ...)
	abline(h=1, col='lightgray')
	# Error limits
	env.lims <- exp(seq(log(xlim[1]/2),log(xlim[2]*2), length.out=100))
	lines(env.lims, exp(sampling.error(env.lims, diff.list$sd.log.cutoff)), col='gray75', ...)
	lines(env.lims, exp(-sampling.error(env.lims, diff.list$sd.log.cutoff)), col='gray75', ...)
	# Points
	points(ms.sub[,abund.f], ms.sub$ratio.hl.normalized.mean, pch=18, col="gray70", cex=0.5*cex.pt, ...)
	points(sig.sub[,abund.f], sig.sub$ratio.hl.normalized.mean, col='gray10', pch=16, cex=cex.pt, ...)
	segments(sig.sub[,abund.f], sig.sub$ratio.hl.normalized.lower.95, sig.sub[,abund.f], sig.sub$ratio.hl.normalized.upper.95, col='black', ...)
	if (is.null(xlim)) {
		my.axis(1, ms.sub[,abund.f], log=T, las=las, ...)
	} else {
		my.axis(1, xlim, log=T, las=las, ...)
	}
	# Labels for proteins
	if (nrow(diff.list$up) > 0 & !is.null(lab.fld)) {
		label.silac(diff.list$up, pos=4, cex=cex.gene.lab, col="gray5")
	}
	if (nrow(diff.list$down) > 0 & !is.null(lab.fld)) {
		label.silac(diff.list$down, pos=4, cex=cex.gene.lab, col="gray5")
	}
}

plot.silac.all <- function(x, count.cutoff=1, lab.fld="gene", abund.f="ratio.hl.count", sd.envelope=0.95, cex.pt=1, cex.gene.lab=0.6, xlim=NULL, xlab="Ratio count", ylab='Protein ratio', ...) {
	ms.sub <- subset(x, ratio.hl.count>0)
	diff.genes <- get.silac.updown(ms.sub, count.cutoff, sd.envelope)
	sig.sub <- rbind(diff.genes$up, diff.genes$down)
	# Limits
	d <- list(...)
	xlim <- d$xlim
	if (is.null(xlim)) {
		xlim <- c(1,max(ms.sub[,abund.f], na.rm=T)*1.1)
	}
	if (is.null(d$ylim)) {
		# Set y limits so that we can always see the full envelope, and also encompass all significantly changed proteins and their error bars.
		ylim <- c(min(exp(-diff.genes$sd.log.cutoff),min(sig.sub$ratio.hl.normalized.lower.95,na.rm=T)),
			max(exp(diff.genes$sd.log.cutoff),max(sig.sub$ratio.hl.normalized.upper.95,na.rm=T)))
	}
	# Initial plot for labels and limits
	plot(ms.sub[,abund.f], ms.sub$ratio.hl.normalized.mean, type='n', log='xy', xaxt='n', pch=18, las=1, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	abline(h=1, col='lightgray')
	# Error limits
	env.lims <- exp(seq(log(xlim[1]/2),log(xlim[2]*2), length.out=100))
	lines(env.lims, exp(sampling.error(env.lims, diff.genes$sd.log.cutoff)), col='gray75', ...)
	lines(env.lims, exp(-sampling.error(env.lims, diff.genes$sd.log.cutoff)), col='gray75', ...)
	# Points
	points(ms.sub[,abund.f], ms.sub$ratio.hl.normalized.mean, pch=18, col="gray70", cex=0.5*cex.pt, ...)
	points(sig.sub[,abund.f], sig.sub$ratio.hl.normalized.mean, col='gray10', pch=16, cex=cex.pt, ...)
	segments(sig.sub[,abund.f], sig.sub$ratio.hl.normalized.lower.95, sig.sub[,abund.f], sig.sub$ratio.hl.normalized.upper.95, col='black', ...)
	# Axes
	if (is.null(xlim)) {
		my.axis(1, ms.sub[,abund.f], log=T, las=1, ...)
	} else {
		my.axis(1, xlim, log=T, las=1, ...)
	}
	# Labels for proteins
	if (nrow(diff.genes$up) > 0 & !is.null(lab.fld)) {
		label.silac(diff.genes$up, pos=4, cex=cex.gene.lab, col="gray5")
	}
	if (nrow(diff.genes$down) > 0 & !is.null(lab.fld)) {
		label.silac(diff.genes$down, pos=4, cex=cex.gene.lab, col="gray5")
	}
	# Return differentially affected proteins
	diff.genes
}

