# input: vector of means for each spline; size of each  spline
# out: quadratic functions for each spline

splineOptimFunction <- function(par, means, partition) {
        # par = list(a, b, c)
        # partition is the right hand x value
        len <- length(partition)
        a <- par$a
        b <- par$b
        c <- par$c
        forwardSeq <- c(2:len, 1)
        aF <- a[forwardSeq]
        bF <- b[forwardSeq]
        cF <- c[forwardSeq]
        x1 <- partition
        x0 <- partition[c(len, 1:(len - 1))]
        a <- (3 / 2) * 
                ((bF - 2 * aF) * (x1^2 - x0^2) + (means + aF - cF) * (x1 - x0)) /
                ((x1^3 - x0^3) - 3 * (x1^2 - x0^2) + 3 * (x1 - x0))
        b <- bF - 2 * (a - aF)
        c <- a - aF + cF
        res <- list(a = a, b = b, c = c)
}

splineOptimCriteria <- function(par, means, partition) {
        # par = list(a, b, c)
        # partition is the right hand x value
        len <- length(partition)
        a <- par$a
        b <- par$b
        c <- par$c
        forwardSeq <- c(2:len, 1)
        aF <- a[forwardSeq]
        bF <- b[forwardSeq]
        cF <- c[forwardSeq]
        x1 <- partition
        x0 <- partition[c(len, 1:(len - 1))]
        minOne <- (1 / (x1 - x0)) * ((1 / 3) * a * (x1^3 - x0^3) +
                                     (1 / 2) * b * (x1^2 - x0^2) +
                                     c * (x1 - x0)) - means
        minTwo <- (a - aF) * x1^2 + (b - bF) * x1 + (c - cF)
        minThree <- 2 * (a - aF) * x1 + (b - bF)
        min <- minOne^2 + minTwo^2 + minThree^2
        sum(min)
}

systemMatrixA <- function(partition, i) {
        n <- length(partition) - 1
        pBackTmp <- partition[i]
        pTmp <- partition[i + 1]
        if (i != n) {
                rowOne <- c(pTmp^2, pTmp, 1, -pTmp^2, -pTmp, -1)
                rowTwo <- c(2 * pTmp, 1, 0, -2 * pTmp, -1, 0)
                rowOne <- c(rep(0, length = 3 * (i - 1)), rowOne, 
                            rep(0, length = 3 * (n - (i + 1))))
                rowTwo <- c(rep(0, length = 3 * (i - 1)), rowTwo, 
                            rep(0, length = 3 * (n - (i + 1))))
        }
        if (i == n) {
                rowOne1 <- c(-pTmp^2, -pTmp, -1)
                rowOne2 <- c(pTmp^2, pTmp, 1)
                rowTwo1 <- c(-2 * pTmp, -1, 0)
                rowTwo2 <- c(2 * pTmp, 1, 0)
                rowOne <- c(rowOne1, rep(0, length = 3 * (n - 2)),
                            rowOne2)
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

ouSplines <- function(means, partition, method = c("system", "optim")) {
        # partition is the right hand x value
        # approximate n quadratic functions 
        # using mean of each function as constraint
        if (all(method == c("system", "optim"))) method = "system"
        n <- length(means)
        if (missing(partition)) partition <- (0:n) 
        init <- numeric(3 * n)
        if (method == "optim") {
                opt <- optim(par = init, function(abc) {
                                     abc <- list(a = abc[1:n], b = abc[(n+1):(2*n)],
                                                 c = abc[(2*n + 1):(3*n)])
                                     abcNew <- splineOptimFunction(par = abc, means = means, 
                                                                 partition = partition[-1])
                                     min <- splineOptimCriteria(par = abcNew, means = means,
                                                               partition = partition[-1])
                                     min
                                             }, method = "BFGS", control = list(maxit = 5000))
                optPar <- opt$par
                optPar <- list(a = optPar[1:n], b = optPar[(n+1):(2*n)],
                         c = optPar[(2*n + 1):(3*n)])
                res <- splineOptimFunction(par = optPar, means = means, 
                                   partition = partition[-1])
                res <- sapply(seq_len(n), function(i) {
                                      tmpSeq <- seq(from = i, by = n, length = 3)
                                      res[tmpSeq]
                                             })
                res <- do.call(cbind, res)
        }
        if (method == "system") {
                A <- lapply(seq_len(n), function(i) {
                                    systemMatrixA(partition = partition, 
                                                  i = i)
                                             })
                A <- t(do.call(cbind, A))
                rownames(A) <- NULL
                b <- unlist(lapply(seq_len(n), function(i) {
                                           mm <- means[i]
                                           c(rep(0, 2), mm)
                                             }))
                res <- solve(a = A, b = b)
        }
        rownames(res) <- paste0("spline", seq_len(n))
        opt$par <- res
        opt$partition <- partition
        opt
}

splineCheck <- function(optPar, means, partition) {
        # check that each RH boadary is consistent
        # par = matrix(ncol = 3, nrow = n)
        # partition is the right hand x value
        len <- length(partition)
        a <- optPar[, 1]
        b <- optPar[, 2]
        c <- optPar[, 3]
        forwardSeq <- c(2:len, 1)
        aF <- a[forwardSeq]
        bF <- b[forwardSeq]
        cF <- c[forwardSeq]
        x1 <- partition[2:len]
        x0 <- partition[c(len, 1:(len - 1))]
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

splinePlot <- function(spline, r = 5, savePlot = TRUE) {
        # r = diameter of circle in plot
        degree <- c(2, 1, 0)
        par <- spline$par
        partition <- spline$partition
        # assume equal partition spacing
        deltaPartition <- partition[2] - partition[1]
        n <- length(partition)
        lenSeq <- 100
        xy <- lapply(seq_len(n), function(i) {
                             pTmp <- partition[i+1]
                             pBackTmp <- partition[i]
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
        theta <- 2 * pi * x / 12
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
        fanLines <- t(sapply(seq_len(n), function(i) {
                                   whichPoint <- (i - 1) * lenSeq + 1
                                   c(xCircle[whichPoint], yCircle[whichPoint])
                        }))
        lenFan <- 2
        for (i in 1:n) {
                lines(x = c(0, lenFan * fanLines[i, 1]),
                      y = c(0, lenFan * fanLines[i, 2]), 
                        col = rgb(0, 0, 0, .5))
        }
        if (savePlot == TRUE) invisible(dev.off())
}
