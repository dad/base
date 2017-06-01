source("fplc.R")

directory <- "../../../Dropbox (Drummond Lab)/pab1-chaperones/data/Fig1"
#directory <- "c:/Dropbox (Drummond Lab)/heat-agg/data/fplc/PAB1 plus chaperones"

files <- list.files(path=directory, pattern='*.asc$')
#print(files)
load.f <- function(f){read.table(paste(directory,f,sep='/'), skip=2, sep='\t', header=T)}

quant.fname <- paste(directory, 'quant.txt', sep='/')
n.written <- 0

file.out <- T
if (file.out) dev.out(paste(directory, 'quant-figures',sep='/'),fdir='', width=5, height=5, output.type='pdf')
par(mar=c(4,4,1,1))
for (fi in files) {
	#print(fi)
	x <- load.f(fi)
	print(fi)
	#print(nrow(x))
	peaks <- identify.peaks(x, min.peak=0.15, min.peak.dist=0.05, min.base=0.0)
	quant <- quantify.peaks(peaks, x)
	quant <- data.frame(fname=rep(fi, nrow(quant)), quant)
	writecol <- fi == files[1]
	write.table(format(quant, digits=3, nsmall=3), file=quant.fname, append=!writecol, quote=F, row.names=F, col.names=writecol, sep='\t')
	n.written <- n.written + nrow(quant)
	plot.peaks(x, peaks, annotate=F, main=fi, xlim=c(min(x$ml,na.rm=T),max(x$ml,na.rm=T)))
	peak.polygon(peaks, x)
	text(peakpos(peaks), peakheight(peaks), 1:length(peaks), pos=3, offset=0.1, col=myrainbow(length(peaks)))
}
if (file.out) dev.off()

cat("# Wrote", n.written, "lines to", quant.fname,'\n')
