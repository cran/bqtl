"make.loc.right"<-
    function(marker.frame, marker.distances)
{
    ncm <- ncol(marker.frame)
    nrm <- nrow(marker.frame)
    res <- matrix(0, nrow = nrm, ncol = ncm)
    dimnames(res) <- list(NULL, names(marker.frame))
    init <- rep(ncm, nrm)
    for(i in rev(1:ncm)) {
        is.terminal <-
            if(i < ncm && marker.distances[i] == 0)
                rep(TRUE, nrm)
            else
                !is.na(marker.frame[, i])
        init[is.terminal] <- i
        res[, i] <- init
    }
    res
}
