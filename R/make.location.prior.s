"make.location.prior"<-
    function(x, add.2.end = 0, normalize = TRUE)
{
###x is a lambda parametter vector
    chromo.break <- x == 0
    map.dist <- map.dx(x)
    len.map <- length(map.dist)
    res <-
        ifelse(chromo.break,
               add.2.end/2 + c(0, map.dist[ - len.map]), 
               map.dist/2) + ifelse(c(TRUE, chromo.break[ - len.map]), 
                                    map.dist,
                                    add.2.end/2 + c(0, map.dist[ - len.map])/2)
    if(normalize)
        res/sum(res)
    else res
}
