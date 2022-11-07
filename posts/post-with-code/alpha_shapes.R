## alpha_shapes.R
## last updated: 2011-04-02

matfile <- "P4_14Jun11_spont_ctl__hotMaps.mat"
alpha <- 6.4

## Type "colors()" at the R prompt to get a list of colournames that you
## can use here.
active.col <- "green"
inside.col <- "black"
outside.col <- "lightgrey"

## Set the following to TRUE/FALSE to state
## if you need PDF or Postscript output.
## 
need.pdf <- TRUE
need.postscript <- FALSE

## Which waves do you want to see plotted?
waves.to.plot <- "all"                  #if you want to see them all.
##waves.to.plot <- c(71)             # if you just want to see a few.

######################################################################
## Nothing else below to edit.
######################################################################

require(R.matlab)
require(alphahull)
require(lattice)
##trellis.device(color = FALSE)

get.coords <- function(X, noise=0.0001) {
  ## X is a 2d matrix.
  ## Return the 2d location of the coordinates, with tiny bit of noise
  ## to stop geometry problems.
  pts <- which(X>0, arr.ind=TRUE)
  pts + matrix(rnorm(nrow(pts)*2)*noise, ncol=2)
}

get.id <- function(X) {
  ## X is a 2d matrix.
  ## Return the ID of each active electrode.
  which(X>0, arr.ind=FALSE)
}


inahullM <- function(hull, X) {
  ## Each row of X is a point P to pass to inahill(hull, P)
  apply(X, 1, function(p) {inahull3(hull, p)})
}

plot.one.wave <- function(n) {
  hull <- hulls[[n]]
  inside <- inahullM(hull, E)
  title <- sprintf("w %d: d = %d/%d = %.2f; a = %.1f\n", n,
                   n.active[n], n.insiders[n], density[n],
                   areas[n])
  ## origin in top-left corner
  plot(hull, asp=1, main=title,
       xaxt="n", yaxt="n",
       xlab="", ylab="",
       xlim=c(0,64), ylim=c(64,0), wpoints=FALSE)
  active <- R[[n]]
  stopifnot(length(active) == n.active[n])
  pt.col <- ifelse(inside, inside.col, outside.col)
  pt.col[active] <- active.col
  points(E, pch=20, col=pt.col, cex=0.4)
  if (FALSE) {
    plot(shapes[[n]], add=T,
         col='blue',wpoints=F,lwd=0.4)  #add the alphashape
  }
  rect(0,0, 64,64, lty=2, lwd=0.5)

  ## Add a scalebar
  segments(0, 66, 10, 66)
}

id.to.row <- function(ids, E.id) {
  ## IDS is set of electrode ids; look them up and
  ## return the corresponding locations.
  rows <- match(ids, E.id)
  if (any(is.na(rows)))
    browser()
  rows
}

row.to.coords <- function(rows, E) {
  E[rows,]
}

inahull3<-function (ahull.obj, p, tol=10e-5) 
{
  ## This function emailed to me from Beatriz, author of alphahull.

compl <- ahull.obj$complement
wh<-which(compl[,3]==ahull.obj$alpha)
compl<-compl[-wh,]
compl<-rbind(ahull.obj$arcs[,1:3],compl[,1:3])


    halfpl <- which(compl[, "r"] < 0)
    n.halfpl <- length(halfpl)
    ball <- which(compl[, "r"] > 0)
    n.ball <- length(ball)
    in.compl <- FALSE
    if (n.halfpl >= 1) {
        h <- 1
        while ((h <= n.halfpl) & in.compl == FALSE) {
            sig = compl[halfpl[h], 3]
            a = compl[halfpl[h], 1]
            b = compl[halfpl[h], 2]
            if (sig <= -3) {
                if (p[1] > (a+tol)) {
                  if (sig == -3) {
                    in.compl <- TRUE
                  }
                }
                else if (p[1] < (a-tol)) {
                  if (sig == -4) {
                    in.compl <- TRUE
                  }
                }
            }
            else {
                if (p[2] > a + b * p[1]+tol) {
                  if (sig == -1) {
                    in.compl <- TRUE
                  }
                }
                else if (p[2] < a + b * p[1]-tol) {
                  if (sig == -2) {
                    in.compl <- TRUE
                  }
                }
            }
            h <- h + 1
        }
    }
    if (in.compl == FALSE) {
        k <- 1
        while ((k <= n.ball) & in.compl == FALSE) {
            r = compl[ball[k], 3]
            c1 = compl[ball[k], 1]
            c2 = compl[ball[k], 2]
            d <- sqrt((p[1] - c1)^2 + (p[2] - c2)^2)
            if (d < r-tol) {
                in.compl <- TRUE
            }
            k <- k + 1
        }
    }
    return(in.ahull = !in.compl)
}

######################################################################
## End of functions
######################################################################


## This is the set of all active electrodes.

csvfile <- gsub(".mat$", "_alpha.csv", matfile)
pdffile <- gsub(".mat$", "_alpha.pdf", matfile)
psfile <- gsub(".mat$", "_alpha.ps", matfile)
pdfolapfile <- gsub(".mat$", "_alpha_olap.pdf", matfile)

matdata <- readMat(matfile)

E <- get.coords(matdata$overallHotMap, noise=0.001)
E.id <- get.id(matdata$overallHotMap)
plot(E, pch=20, main="electrode positions in recording")


hotmaps <- matdata$hotMaps
names(hotmaps) <- NULL  ## important - else levelplots won't display! with NAs.
##length(hotmaps) <- 10
n.waves <- length(hotmaps)
cat(sprintf("About to analyse %d waves\n", n.waves))
##I <- lapply(hotmaps, get.id)
get.id2 <- function(hotmap) {
  m = hotmap[[1]]
  a = m[[1]]
  dim(a)
  get.id(m)
}
get.id2(hotmaps[[2]])

I <- lapply(hotmaps, get.id2)

## Derive the position of each wave from the E map, which has
## been jittered.
R <- lapply(I, id.to.row, E.id)
C <- lapply(R, row.to.coords, E)


## C <- lapply(hotmaps, get.coords)

n.active <- sapply(I, length)           #number of active electrodes

hulls  <- lapply(C, ahull, alpha=alpha)
shapes <- lapply(C, ashape, alpha=alpha)
areas  <- sapply(hulls, areaahull)
insiders <-sapply(hulls, function(hull) { which(inahullM(hull, E))})
n.insiders <- sapply(insiders, length)
density <- round(n.active/n.insiders,3)

cat("About to calculate overlaps\n")
shared.electrodes <- shared.insiders <- matrix(NA, nrow=n.waves, ncol=n.waves)
for (i in 1:(n.waves-1)) {
  for (j in (i+1):n.waves) {
    shared.electrodes[i,j] <- length(intersect(I[[i]], I[[j]]))
    shared.insiders[i,j] <- length(intersect(insiders[[i]], insiders[[j]]))
  }
}


min.active <- outer(n.active, n.active, pmin)
min.insiders <- outer(n.insiders, n.insiders, pmin)
shared.electrodes.frac <-  shared.electrodes/min.active
shared.insiders.frac <-    shared.insiders/min.insiders

## adapted from ?image, with volcano.
##f <- function(m) {t(m)[ncol(m):1,]}
f <- function(m) {m}

p1 <- levelplot(f(shared.electrodes), main="shared.active (n)")
p2 <- levelplot(f(shared.insiders), main="shared.insiders (n)")
p3 <- levelplot(f(shared.electrodes.frac), main='shared.active (fraction)')
p4 <- levelplot(f(shared.insiders.frac), main='shared.insiders (fraction)')

pdf(file=pdfolapfile, width=8, height=8)
print(p1, position=c(0.0, 0.0, 0.5, 0.5), more=TRUE)
print(p2, position=c(0.5, 0.0, 1.0, 0.5), more=TRUE)
print(p3, position=c(0.0, 0.5, 0.5, 1.0), more=TRUE)
print(p4, position=c(0.5, 0.5, 1.0, 1.0), more=FALSE)
dev.off()


#######################################################################
## Which waves do we want to plot?
if (is.character(waves.to.plot) && (waves.to.plot == "all")) {
  waves.to.plot = 1:n.waves
}

## qpdf can compress:
## qpdf --stream-data=compress P4__hotMaps_alpha.pdf a.pdf

if (need.pdf) {
  pdf(file=pdffile,width=8, height=11)
  par(mfrow=c(3,2), xaxt="n", yaxt="n")
  for (i in waves.to.plot) {
    plot.one.wave(i)
  }
  dev.off()
}

if (need.postscript) {
  postscript(file=psfile,width=8, height=11, onefile=FALSE, horizontal=FALSE)
  par(mfrow=c(3,2), xaxt="n", yaxt="n")
  for (i in waves.to.plot) {
    plot.one.wave(i)
  }
  dev.off()
}





df <- data.frame(wave=1:n.waves,
                 area=round(areas,3),
                 n.active=n.active,
                 n.inside=n.insiders,
                 density=density)
write.csv(df, file=csvfile, row.names=FALSE)

## area and n.inside are well-correlated.
plot(df$area, df$n.inside, xlab='alphahull area', ylab='#electrodes inside')


