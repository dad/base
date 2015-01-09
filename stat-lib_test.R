source("stat-lib.R")

c1 <- circle(c(1,1),r=0.5)
c2 <- circle(c(2,2),r=0.5)
c3 <- circle(c(2,2),r=1)

stopifnot(collides(c1,c3))
stopifnot(collides(c2,c3))
stopifnot(!collides(c1,c2))
stopifnot(collides(c1,list(c2,c3)))
