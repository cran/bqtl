"marker.fill"<-
    function(map.frame, reso , return.nint=FALSE)
{
### given some markers and a desired resolution
### add pseudo-markers at 'good' distances
### and create names that show their relative
### positions
### map.frame$mgd is the distance in M or cM from the left telomere
### reso is in same metric
    int.suffix <- function(nint)
    {
        res <- sapply(nint, function(x)
                  {
                      strng <- if(x > 1) {
                          c("", format((1:(x - 1))/10^(trunc(log10(x - 1) +
                                                             1))))
                      }
                      else {
                          ""
                      }
                      substring(strng, 2)
                  }
                      )
        unlist(res)
    }
    
### assume that map.frame is a map.frame.object
### start of each chromosome
    if (any(is.na(pmatch(c("cM","marker.name","chr.num"),names(map.frame)))))
        stop("need $cM, $marker.name, and $chr.num components in first arg")
    dx <- map.frame$cM
    names(dx) <-   map.frame$marker.name
    mgd <- tapply(dx, map.frame$chr.num,c)
    names(mgd) <- rep("",length(mgd))
    
    res.main <- lapply(mgd, function(x, reso, int.suffix, return.nint)
                   {
                       genes <- grep("gene", names(x))
                       if(length(genes) > 0)
                           mark.loc <- x[ - genes]
                       else mark.loc <- x
###     nint <- pmax(1, (diff(mark.loc)+reso/2) %/% reso)
###                       nint <- pmax(1, ceiling(diff(mark.loc) / reso))
                       nint <- if (length(mark.loc)==1) integer(0) else pmax(1, ceiling(diff(mark.loc)/reso))
                       if(return.nint) {
                           nint <- c(nint, 1)
                           names(nint) <- names(mark.loc)
                           return(nint)
                       }
                       else
                       {
                           dxes <- c(rep(diff(mark.loc)/nint, nint), Inf)
                           names(dxes) <-
                               paste(rep(names(mark.loc), c(nint, 1)), 
                                     int.suffix(c(nint, 1)), sep = "")
                           return(dxes)
                       }
                   }
                       , reso, int.suffix,return.nint)
    unlist(res.main)
}
