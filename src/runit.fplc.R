library(RUnit)
source("fplc.R")

test.max <- function() {
	mult <- 1.0
	res <- 1/100; 
	s <- seq(1,30,res); 
	d <- data.frame(ml=s, mAU=mult*(dnorm(s, mean=2, sd=0.1)+dnorm(s, mean=4, sd=0.1)+dnorm(s, mean=7, sd=1)))
	ex <- find.extrema(d, meth='max')
	checkEqualsNumeric(ex[1,'x'],2)
	checkEqualsNumeric(ex[2,'x'],4)
	checkEqualsNumeric(ex[3,'x'],7)
}

test.min <- function() {
	mult <- 1.0
	res <- 1/100; 
	s <- seq(1,5,res); 
	d <- data.frame(ml=s, mAU=mult*(dnorm(s, mean=2, sd=0.5)+dnorm(s, mean=4, sd=0.5)))
	#plot(d)
	ex <- find.extrema(d, meth='min', bw=10)
	#print(ex)
	checkEqualsNumeric(ex[1,'x'],3)
}

test.find.nearest <- function() {
	mult <- 5.0
	res <- 1/100; 
	s <- seq(0,3,res); 
	d <- data.frame(ml=s, mAU=mult*(dnorm(s, mean=2, sd=0.5)))
	d2 <- data.frame(ml=s, mAU=mult*(dnorm(s, mean=2, sd=0.5))-1)
	#plot(d)
	plist1 <- identify.peaks(d)
	plist2 <- identify.peaks(d2)
	checkEquals(length(plist1), 1)
	checkEquals(length(plist2), 1)
	cl1 <- find.closest(d, plist1[[1]]$x, target.y=0.0)
	cl2 <- find.closest(d2, plist2[[1]]$x, target.y=0.0)
	#print(cl)
	#abline(h=0.0)
	#abline(v=c(cl$lower, cl$upper))
	checkTrue(cl2$lower>cl1$lower)
	checkTrue(cl2$upper<cl1$upper)
}

test.requant.peaks <- function() {
	mult <- 1.0
	res <- 1/100; 
	s <- seq(1,30,res); 
	d <- data.frame(ml=s, mAU=mult*(dnorm(s, mean=2, sd=0.1)+dnorm(s, mean=4, sd=0.1)+dnorm(s, mean=7, sd=1)))
	plist <- identify.peaks(d)
	q1 <- quantify.peaks(plist, d)
	d2 <- d
	mult.diff <- 0.1
	d2$mAU <- mult.diff*d2$mAU
	q2 <- quantify.peaks(plist, d2)
	for (xi in 1:nrow(q1)) {
		checkEqualsNumeric(q2[xi,'auc'], mult.diff*q1[xi,'auc'])
	}
}

test.peak.center <- function() {
	mult <- 1.0
	res <- 1/10; s <- seq(1,20,res); d <- data.frame(ml=s, mAU=mult*dnorm(s, mean=10, sd=1))
	p <- identify.peaks(d)
	p1 <- p[[1]]
	checkEqualsNumeric(10, p1$x)
}

test.area.under.curve <- function() {
	mult <- 1.0
	res <- 1/10; 
	s <- seq(1,20,res); 
	d <- data.frame(ml=s, mAU=mult*dnorm(s, mean=10, sd=1))
	p <- identify.peaks(d, min.base=NULL)
	p1 <- p[[1]]
	auc <- area.under.curve(p1, d)
	checkEqualsNumeric(mult, auc)
}

test.filter.peaks <- function() {
	p1 <- list(x=10, y=11, lower=9, upper=11)
	p2 <- list(x=15, y=10, lower=14, upper=16)
	plist <- list(p1,p2)
	pl <- filter.peak.list(plist)
	checkEquals(length(pl),2)
	pl <- filter.peak.list(plist, min.peak.dist=6)
	checkEquals(length(pl),1) # filter out
	checkEqualsNumeric(pl[[1]]$y, 11) # check to see the largest peak is retained
}

test.peaks.not.too.close <- function() {
	mult <- 1.0
	res <- 1/10; 
	s <- seq(1,30,res); 
	d <- data.frame(ml=s, mAU=mult*(dnorm(s, mean=10, sd=1)+dnorm(s, mean=20, sd=1)))
	p <- identify.peaks(d, min.peak.dist=0.1)
	checkEquals(length(p),2)
}

test.peaks.too.close <- function() {
	mult <- 1.0
	res <- 1/10; 
	s <- seq(1,30,res); 
	d <- data.frame(ml=s, mAU=mult*(dnorm(s, mean=10, sd=1)+dnorm(s, mean=13, sd=1)))
	p <- identify.peaks(d, min.peak.dist=5)
	checkEquals(length(p),1)
}
