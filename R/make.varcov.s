"make.varcov"<-
    function(regressor.matrix, y, subset = is.finite(y), casewt = NULL)
{
    df <-
        if (is.logical(subset))
            sum(subset) - 1
        else
            length(y[subset])-1
    if(is.null(casewt)) {
        var.x <- var(regressor.matrix[subset,  ]) * df
        cov.xy <- var(regressor.matrix[subset,  ], y[subset]) * df
        var.y <- var(y[subset]) * df
    }
    else {
        casewt <- as.vector(casewt[subset])
        mean.x <- matrix(casewt/sum(casewt), nrow = 1) %*%
            as.matrix(regressor.matrix[subset,  ])
        mean.y <- sum(y[subset] * casewt)/sum(casewt)
        center.x <-
            as.matrix(regressor.matrix[subset,  ]) -
                rep(mean.x, rep(sum(subset), length(mean.x)))
        center.y <- y[subset] - mean.y
        var.x <- t(center.x) %*% (casewt * center.x)
        cov.xy <- t(center.x) %*% (casewt * center.y)
        var.y <- sum(casewt * center.y^2)
    }
    res <- list(var.x = var.x, var.y = var.y, cov.xy = cov.xy, df = df)
    attr(res,"call") <- match.call()
    class(res) <- "varcov"
    res
}
