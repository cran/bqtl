"make.map.frame"<-
    function(dx,chr.num=NULL,
             prior=make.location.prior(lambda),morgan=100,nint=NULL,reso=NULL)
{
    if (is.null(nint)&&!is.null(reso))
        nint <- marker.fill(make.map.frame(dx),reso, TRUE )
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
    if (is.data.frame(dx)) {
        ## allow a data.frame with the first column containing the distances
        if (ncol(dx)>1) {
            dx.synonyms <- c("cM","cm","dx","M","m")
            dx.name <- dx.synonyms[is.element(dx.synonyms,names(dx))][1]
            if ( length(dx.name)!=0 && dx.name!="NA" && !is.na(dx.name) ) {# found explicitly labelled distance
                marker.names <-
                    if (any(pmatch(names(dx),"marker.names",0)!=0))
                        dx$marker.name
                    else {
                        warning("could not find marker.names using row.names(dx) instead")
                        row.names(dx)
                    }
                if (is.element(dx.name,c("m","M"))){
                    morgan <- 100
                    mgd <- dx[,dx.name]/100
                }
                else
                    mgd <- dx[,dx.name]
                if ( missing(chr.num) && ("chr.num"%in%names(dx)) )
                    chr.num <- dx$chr.num
                
            }
            else {
                warning(paste("using first column",names(dx)[[1]], "as map distance"))
                marker.names <- row.names(dx)
                mgd <- dx[[1]]
            }
        }
        else { # 1 column data.frame
            marker.names <- row.names(dx)
            mgd <- dx[[1]]
        }
    }
    else { #not a  data.frame
        mgd <- dx 
        if (is.null(names(dx)))
            names(dx) <- paste("marker",seq(along=chr.num),sep=".")
        marker.names <-  names(dx)
    }
    chr.num <-
        if (is.null(chr.num))
            cumsum( diff(c(Inf,mgd))<0 )
        else
            eval(chr.num)
### check ordering within chromos
    mk.order <- order(chr.num,mgd)
    if (any(out.order<-
            tapply(mgd,chr.num,function(x) length(x)>1 && any(diff(x)<0))))
        stop(paste("markers not in increasing order on chr:",
                   paste(names(out.order)[out.order],collapse=",")))
    
    ##  check arg lengths
    if (!all.equal( length(nint), length(mgd), length(chr.num) ) )
        stop("lengths of nint, chr.num, dx not all equal")
    ## default is 'right' if only one marker on a chromo
    pos.type <- as.numeric(c(0,chr.num[-length(chr.num)])==chr.num) +
        as.numeric(c(chr.num[-1],0)==chr.num)*2
    pos.type <- factor(c(1,1,2,3)[pos.type+1],1:3,c("right","left","center"))
    is.marker <- rep(TRUE,length(chr.num))
    
    increments <- diff(mgd/morgan)

    not.r <- ifelse(pos.type=="right",0.0,
                    c(exp(-increments),0.0))
    if (!is.null(nint)) {
        orig <- cumsum(c(1,nint[-length(nint)]))
        not.r <- rep(not.r^(1/nint),nint)
        marker.names <- paste(rep(marker.names,nint),int.suffix(nint),sep="")
        is.marker <- rep(is.marker,nint)
        is.marker[-orig] <- FALSE
        chr.num <- rep(chr.num,nint)
        mgd <-
            rep(mgd,nint) +
                morgan*rep(c(increments,0)/nint,nint)*
                    (unlist(sapply(nint,seq))-1)
        pos.type <- rep(pos.type,nint)
        pos.type[-orig] <- "center"
    }
    lambda <- ifelse(not.r==0,0,2*not.r-1)
    prior <- eval(prior)
    pos.plot <- mgd
    locus<-paste("C",chr.num,as.numeric(format(100*mgd/morgan,digits=4)),
                 sep=".")
    res <- data.frame(
                      marker.name=
                      I(as.character(make.names(marker.names,TRUE))),
                      cM=as.numeric(mgd),
                      prior=as.numeric(prior),
                      pos.type=I(pos.type),
                      is.marker=I(is.marker),
                      pos.plot=as.numeric(pos.plot),
                      lambda=as.numeric(lambda),
                      locus=I(locus),
                      chr.num=I(chr.num))
    attr(res,"morgan") <- morgan
    class(res)<- c("map.frame",class(res))
    res
}
