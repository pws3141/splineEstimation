# input: vector of means for each spline; size of each  spline
# out: quadratic functions for each spline

splineOptimFunction <- function(par, means, breaks) {
        # par = list(a, b, c)
        # partition is the right hand x value
        len <- length(breaks)
        a <- par$a
        b <- par$b
        c <- par$c
        forwardSeq <- c(2:len, 1)
        aF <- a[forwardSeq]
        bF <- b[forwardSeq]
        cF <- c[forwardSeq]
        x1 <- breaks
        x0 <- breaks[c(len, 1:(len - 1))]
        a <- (3 / 2) * 
                ((bF - 2 * aF) * (x1^2 - x0^2) + (means + aF - cF) * (x1 - x0)) /
                ((x1^3 - x0^3) - 3 * (x1^2 - x0^2) + 3 * (x1 - x0))
        b <- bF - 2 * (a - aF)
        c <- a - aF + cF
        res <- list(a = a, b = b, c = c)
}

splineOptimCriteria <- function(par, means, breaks) {
        # par = list(a, b, c)
        # partition is the right hand x value
        len <- length(breaks)
        a <- par$a
        b <- par$b
        c <- par$c
        forwardSeq <- c(2:len, 1)
        aF <- a[forwardSeq]
        bF <- b[forwardSeq]
        cF <- c[forwardSeq]
        x1 <- breaks
        x0 <- breaks[c(len, 1:(len - 1))]
        minOne <- (1 / (x1 - x0)) * ((1 / 3) * a * (x1^3 - x0^3) +
                                     (1 / 2) * b * (x1^2 - x0^2) +
                                     c * (x1 - x0)) - means
        minTwo <- (a - aF) * x1^2 + (b - bF) * x1 + (c - cF)
        minThree <- 2 * (a - aF) * x1 + (b - bF)
        min <- minOne^2 + minTwo^2 + minThree^2
        sum(min)
}

systemMatrixA <- function(breaks, i) {
        n <- length(breaks) - 1
        pBackTmp <- breaks[i]
        pTmp <- breaks[i + 1]
        deltaP <- pTmp - pBackTmp
        if (i < n) {
                rowOne <- c(pTmp^2, pTmp, 1, -pTmp^2, -pTmp, -1)
                rowTwo <- c(2 * pTmp, 1, 0, -2 * pTmp, -1, 0)
                rowOne <- c(rep(0, length = 3 * (i - 1)), rowOne, 
                            rep(0, length = 3 * (n - (i + 1))))
                rowTwo <- c(rep(0, length = 3 * (i - 1)), rowTwo, 
                            rep(0, length = 3 * (n - (i + 1))))
        }
        if (i == n) {
                pOne <- breaks[1]
                rowOne1 <- c(-pOne^2, -pOne, -1)
                #rowOne1 <- c(1, 0, 0)
                rowOne2 <- c(pTmp^2, pTmp, 1)
                #rowOne2 <- c(-1, 0, 0)
                rowTwo1 <- c(-2 * pOne, -1, 0)
                #rowTwo1 <- c(0, 1, 0)
                rowTwo2 <- c(2 * pTmp, 1, 0)
                #rowTwo2 <- c(0, -1, 0)
                #rowThree1 <- c(0, 0, 1)
                #rowThree2 <- c(0, 0, -1)
                #rowOne <- c(rowOne1, rep(0, length = 3 * (n - 2)),
                            #rowOne2)
                rowOne <- c(rowOne1, rep(0, length = 3 * (n - 2)),
                            rowOne2)
                #rowTwo <- c(rowTwo1, rep(0, length = 3 * (n - 2)),
                            #rowTwo2)
                rowTwo <- c(rowTwo1, rep(0, length = 3 * (n - 2)),
                            rowTwo2)
        }
        rowThree <- (1 / (pTmp - pBackTmp)) * 
                        c((1 / 3) * (pTmp^3 - pBackTmp^3),
                          (1 / 2) * (pTmp^2 - pBackTmp^2), 1)
        rowThree <- c(rep(0, length = 3 * (i - 1)), rowThree, 
                    rep(0, length = 3 * (n - i)))
        matTmp <- cbind(rowOne, rowTwo, rowThree)
        matTmp
}

ouSplines <- function(means, breaks, method = c("system", "optim")) {
        # partition is the right hand x value
        # approximate n quadratic functions 
        # using mean of each function as constraint
        if (all(method == c("system", "optim"))) method = "system"
        n <- length(means)
        if (missing(breaks)) breaks <- (0:n) 
        if (method == "optim") {
                init <- numeric(3 * n)
                opt <- optim(par = init, function(abc) {
                                     abc <- list(a = abc[1:n], b = abc[(n+1):(2*n)],
                                                 c = abc[(2*n + 1):(3*n)])
                                     abcNew <- splineOptimFunction(par = abc, means = means, 
                                                                 breaks = breaks[-1])
                                     min <- splineOptimCriteria(par = abcNew, means = means,
                                                               breaks = breaks[-1])
                                     min
                                             }, method = "BFGS", control = list(maxit = 5000))
                optPar <- opt$par
                optPar <- list(a = optPar[1:n], b = optPar[(n+1):(2*n)],
                         c = optPar[(2*n + 1):(3*n)])
                res <- splineOptimFunction(par = optPar, means = means, 
                                   breaks = breaks[-1])
                res <- sapply(seq_len(n), function(i) {
                                      tmpSeq <- seq(from = i, by = n, length = 3)
                                      res[tmpSeq]
                                             })
                res <- do.call(cbind, res)
        }
        if (method == "system") {
                A <- lapply(seq_len(n), function(i) {
                                    systemMatrixA(breaks = breaks, 
                                                  i = i)
                                             })
                A <- t(do.call(cbind, A))
                rownames(A) <- NULL
                b <- unlist(lapply(seq_len(n), function(i) {
                                           mm <- means[i]
                                           c(rep(0, 2), mm)
                                             }))
                #b <- c(b, 0, 0, 0)
                res <- solve(a = A, b = b)
                res <- matrix(res, ncol = 3, byrow=TRUE)
                colnames(res) <- c("a", "b", "c")
                #res <- res[-(n+1),]
                opt <- list()
        }
        rownames(res) <- paste0("spline", seq_len(n))
        opt$par <- res
        opt$breaks <- breaks
        opt
}

splineCheck <- function(optPar, means, breaks) {
        # check that each RH boadary is consistent
        # par = matrix(ncol = 3, nrow = n)
        # breaks is the right hand x value
        len <- length(breaks) - 1
        a <- optPar[, 1]
        b <- optPar[, 2]
        c <- optPar[, 3]
        forwardSeq <- c(2:len, 1)
        aF <- a[forwardSeq]
        bF <- b[forwardSeq]
        cF <- c[forwardSeq]
        x1 <- breaks[2:(len + 1)]
        x0 <- breaks[1:len]
        yValRHBoundary <- a * x1^2 + b * x1 + c
        yValLHBoundary <- aF * x1^2 + bF * x1 + cF
        yValDiff <- yValRHBoundary - yValLHBoundary
        dxValRHBoundary <- 2 * a * x1 + b 
        dxValLHBoundary <- 2 * aF * x1 + bF 
        dxValDiff <- dxValRHBoundary - dxValLHBoundary
        meanVal <- (1 / (x1 - x0)) * ((a / 3) * (x1^3 - x0^3) + 
                                      (b / 2) * (x1^2 - x0^2) +
                                      c * (x1 - x0))
        meanValDiff <- meanVal - means
        yVal <- cbind(yValRHBoundary, yValLHBoundary, yValDiff)
        colnames(yVal) <- c("RH", "LH", "Diff")
        dxVal <- cbind(dxValRHBoundary, dxValLHBoundary, dxValDiff)
        colnames(dxVal) <- c("RH", "LH", "Diff")
        meanVal <- cbind(meanVal, means, meanValDiff)
        colnames(meanVal) <- c("intMean", "mean", "Diff")
        res <- list(meanVal = meanVal, yVal = yVal, dxVal = dxVal)
        res
}

# splinePlot aim: plot round the circle?? Is this possible (for me to do)?
# would certainly look nice
# Q: what happens at partitionEnd boundary????

splinePlot <- function(spline, r = 5, circle = TRUE, savePlot = TRUE) {
        # r = diameter of circle in plot
        degree <- c(2, 1, 0)
        par <- spline$par
        breaks <- spline$breaks
        # assume equal partition spacing
        deltaPartition <- breaks[2] - breaks[1]
        n <- length(breaks) - 1
        lenSeq <- 100
        xy <- lapply(seq_len(n), function(i) {
                             pTmp <- breaks[i+1]
                             pBackTmp <- breaks[i]
                        tmpSeq <- seq(from = pBackTmp, to = pTmp, length = lenSeq)
                        tmpSeq <- tmpSeq[-lenSeq]
                        xQuadratic <- outer(tmpSeq, degree, "^")
                        splineTmp <- par[i,]
                        yTmp <- as.vector(tcrossprod(splineTmp, xQuadratic))
                        #yTmp
                        list(x = tmpSeq, y = yTmp)
                                             })
        x <- unlist(lapply(seq_len(n), function(i) xy[[i]]$x))
        y <- unlist(lapply(seq_len(n), function(i) xy[[i]]$y))
        if (circle == TRUE) {
        theta <- 2 * pi * x / max(breaks)
        xCircle <- cos(theta) * r
        yCircle <- sin(theta) * r
        plotScale <- 1/15
        R <- r + plotScale * y
        pname <- "splineCirclePlot.pdf"
        if (savePlot == TRUE) pdf(pname)
                par(mai = c(0, 0, 0, 0))
                plot(xCircle, yCircle, t = "l", asp = 1, #col = "light gray",
                                xlim = c(-10, 10), ylim = c(-10, 10),
                                col = rgb(0, 0, 0, .5),
                                axes = FALSE, xlab = "", ylab = "")
                lines(cos(theta) * R, sin(theta) * R, lwd = 2)
                fanLinesTheta <- 2 * pi * breaks / max(breaks)
                xFanLines <- cos(fanLinesTheta) * r
                yFanLines <- sin(fanLinesTheta) * r
                lenFan <- 2
                for (i in 1:n) {
                        lines(x = c(0, lenFan * xFanLines[i]),
                              y = c(0, lenFan * yFanLines[i]), 
                                col = rgb(0, 0, 0, .5))
        }
        }
        if (circle == FALSE) {
                pname <- "splinePlot.pdf"
                if (savePlot == TRUE) pdf(pname)
                #par(mai = c(1, 1, 1, 1))
                plot(x, y, t = "l", lwd = 2,
                                col = rgb(0, 0, 0, 1),
                                axes = TRUE, xlab = "", ylab = "")
                abline(v = breaks, col = rgb(0, 0, 0, .5))
        }
        if (savePlot == TRUE) invisible(dev.off())
}
