"%equiv%" <-
    function(x,y)
{
    if (!is.matrix(x)) dim(x) <- c(length(x),1)
    if (!is.matrix(y)) dim(y) <- c(length(y),1)
    dm <- c(ncol(x),ncol(y));
    nr <- nrow(x);
    if (nr!=nrow(y)) stop("nrow(x) != nrow(y)")
    t.x <- x[,rep(seq(dm[1]),dm[2]),drop=F]
    t.y <- y[,rep(seq(dm[2]),rep(dm[1],dm[2])),drop=F]
    matrix(apply(t.x==t.y,2,sum) == nr,nc=dm[2],dimnames=list(dimnames(x)[[2]],dimnames(y)[[2]]))
}
