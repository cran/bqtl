"summary.map.frame" <-
    function(x)
{
    chrs<-table(x$chr.num,
                factor(x$is.marker,c(TRUE,FALSE),c("marker","pseudo-marker")))
    lens <- do.call("rbind",tapply(x$cM,x$chr.num,range))
    res <- cbind(chrs,lens)
    dimnames(res)[[2]][3:4]<-c("cM.low","cM.high")
    res
}

