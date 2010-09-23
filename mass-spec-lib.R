sampling.error <- function(n, log.sd.cutoff) {
	log.sd.cutoff/sqrt(n)
}

plot.silac.updown <- function(x, count.cutoff, up.orfs, down.orfs, lab.fld, abund.f="intensity", sd.envelope=0.995, cex.pt=1, cex.gene.lab=0.6, ...) {
	ms.sub <- subset(x, ratio.hl.count>=count.cutoff)
	sig.sub <- subset(ms.sub, orf %in% c(up.orfs,down.orfs))
	up.sig.sub <- subset(ms.sub, orf %in% up.orfs)
	down.sig.sub <- subset(ms.sub, orf %in% down.orfs)
	plot(ms.sub[,abund.f], ms.sub$ratio.hl.normalized.mean, type='n', log='xy', xaxt='n', pch=16, las=1, ...)
	abline(h=1, col='lightgray')
	# Error limits
	f <- function(n, log.sd.cutoff) {
		log.sd.cutoff/sqrt(n)
	}
	sd.logs <- subset(x, ratio.hl.count>2)$ratio.hl.normalized.sd
	sd.log.prop <- sd.envelope
	sd.log.cutoff <- sd.logs[order(sd.logs)][ceiling(length(sd.logs)*sd.log.prop)]
	env.lims <- exp(seq(log(abund.xlim[1]/2),log(abund.xlim[2]*2), length.out=100))
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
		text(x=up.sig.sub[,abund.f], y=up.sig.sub$ratio.hl.normalized.mean, up.sig.sub[,lab.fld], pos=4, cex=cex.gene.lab, col="gray5")
	}
	if (nrow(down.sig.sub) > 0 & !is.null(lab.fld)) {
		text(x=down.sig.sub[,abund.f], y=down.sig.sub$ratio.hl.normalized.mean, labels=down.sig.sub[,lab.fld], pos=4, cex=cex.gene.lab, col="gray5")
	}
	list(all=ms.sub, up=up.sig.sub, down=down.sig.sub)
}
