"summary.map.frame" <-
    function(object,...)
{
    chrs<-table(object$chr.num,
                factor(object$is.marker,c(TRUE,FALSE),c("marker","pseudo-marker")))
    lens <- do.call("rbind",tapply(object$cM,object$chr.num,range))
    res <- cbind(chrs,lens)
    dimnames(res)[[2]][3:4]<-c("cM.low","cM.high")
    res
}

