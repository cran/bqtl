"make.regressor.matrix"<-
    function(state.matrix, mode.mat = NULL)
{
    dsm <- dim(state.matrix)
    dsn <- dimnames(state.matrix)[[2]]
    dimnames(state.matrix) <- NULL
    if(is.null(mode.mat))
        if(dsm[3] == 3) {
            ## assume F2 with both additive and dominance effects
            mode.mat <- cbind(add = c(-1, 0, 1), dom = c(-1, 1, -1))
        }
        else {
            ## assume BC1
            mode.mat <- matrix(c(-1, 1), nrow = 2)
        }
    if(nrow(mode.mat) == 1) {
        dim(state.matrix) <- c(dsm[1] * dsm[2], dsm[3])
        regressor.matrix <- state.matrix %*% mode.mat
        dim(regressor.matrix) <- dsm[1:2]
        dimnames(regressor.matrix) <- list(NULL, dsn)
    }
    else {
        dim(state.matrix) <- c(dsm[1] * dsm[2], dsm[3])
        regressor.matrix <- state.matrix %*% mode.mat
        dimnames(regressor.matrix) <- NULL
        dim(regressor.matrix) <- c(dsm[1:2], ncol(mode.mat))
        regressor.matrix <- aperm(regressor.matrix, c(1, 3, 2))
        mode.names <-  dimnames(mode.mat)[[2]]
        dim(regressor.matrix) <- c(dsm[1], dsm[2] * ncol(mode.mat))
        dimnames(regressor.matrix) <-
            list(NULL,
                 if (is.null(mode.names)) dsn else
                 paste(rep(mode.names, 
                           dsm[2]), rep(dsn, rep(ncol(mode.mat), dsm[2])),
                       sep = "."))
    }
    regressor.matrix
}

