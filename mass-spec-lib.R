
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
		text(x=up.sig.sub[,abund.f], y=up.sig.sub$ratio.hl.normalized.mean, up.sig.sub[,lab.fld], pos=4, cex=cex.gene.lab, col="gray5")
	}
	if (nrow(down.sig.sub) > 0 & !is.null(lab.fld)) {
		text(x=down.sig.sub[,abund.f], y=down.sig.sub$ratio.hl.normalized.mean, labels=down.sig.sub[,lab.fld], pos=4, cex=cex.gene.lab, col="gray5")
	}
	list(all=ms.sub, up=up.sig.sub, down=down.sig.sub)
}

get.silac.updown <- function(x, count.cutoff=2, sd.envelope=0.995) {
	ms.sub <- subset(x, ratio.hl.count>=count.cutoff)
	sd.logs <- subset(x, ratio.hl.count>max(2,count.cutoff))$ratio.hl.normalized.sd
	sd.log.prop <- sd.envelope
	sd.log.cutoff <- sd.logs[order(sd.logs)][ceiling(length(sd.logs)*sd.log.prop)]
	
	ms.up <- subset(ms.sub, ratio.hl.normalized>1 & ratio.hl.normalized.lower.95>exp(sd.log.cutoff/sqrt(ratio.hl.count)))
	ms.down <- subset(ms.sub, ratio.hl.normalized<1 & ratio.hl.normalized.upper.95<exp(-sd.log.cutoff/sqrt(ratio.hl.count)))
	list(up=ms.up[order(ms.up$ratio.hl.normalized, decreasing=TRUE),], down=ms.down[order(ms.down$ratio.hl.normalized, decreasing=FALSE),], sd.log.cutoff=sd.log.cutoff)
}

