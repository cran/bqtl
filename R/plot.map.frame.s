"plot.map.frame"<-
    function(x,y, fun = if (y.type=="matrix") matlines else lines,
             type="l", include.rug=TRUE, rug.lwd=0.1,title.string=NULL,
             y.range=NULL,
             ylab=deparse(substitute(y)),xlab="Location", ...)
{
    mf.obj <- x
    x.as.text <- deparse(substitute(x))
    if (any(is.na(pmatch(c("chr.num","pos.plot"),names(mf.obj)))))
        stop("need $chr.num and $pos.plot components in first arg")
### try to be smart: guess if this is a vector or matrix arg.
    if (missing(y)){
        y.range <- range(mf.obj$chr.num)
        x.range <- range(unlist(mf.obj$cM))
        plot(x.range,y.range,type="n",xlab=xlab,
             ylab= if (missing(ylab)) "Chromosome" else ylab)
        for (i in unique(mf.obj$chr.num)) {
            i.index <- i==mf.obj$chr.num
            x <- mf.obj$pos.plot[i.index]
            x.range <- range(x)
            lines(x.range,c(i,i))
            seg.col <- c(1,2)
            seg.lty <- c(1,3)
            if (any( which.seg <- mf.obj$is.marker[i.index] )) {
                n.seg <- sum(which.seg)
                segments(x[which.seg],rep(i,n.seg)-0.1,
                         x[which.seg],rep(i,n.seg)+0.1,lty=seg.lty[1],
                         col=seg.col[1])}
            if (any( which.seg <- !mf.obj$is.marker[i.index] )) {
                n.seg <- sum(which.seg)
                segments(x[which.seg],rep(i,n.seg)-0.1,
                         x[which.seg],rep(i,n.seg)+0.1,lty=seg.lty[2],
                         col=seg.col[2])}
        }
        title(paste(title.string,x.as.text))
        return(invisible())
    }
    
    y.type <-
        if ( length(dmm <- dim(y)) && max(dmm)<length(y))
            "matrix"
        else
            "vector"
    if (is.null(y)) stop("NULL value for y not accepted in plot.map.frame")
    if (is.null(y.range))
        y.range <- range(y,na.rm=TRUE)

    if (exists("is.R") && is.R() && y.type=="matrix") { # S compat workaround
        matlines <-
            function(x, y, type = 'l', lty=1:5, lwd = 1, pch=NULL, col=1:6, ...)
                matplot(x=x, y=y, type=type, lty=lty, lwd=lwd, pch=pch, col=col,
                        add=TRUE, ...)
        matpoints <-
            function(x, y, type = 'p', lty=1:5, lwd = 1, pch=NULL, col=1:6, ...)
                matplot(x=x, y=y, type=type, lty=lty, lwd=lwd, pch=pch, col=col,
                        add=TRUE, ...)
    }

    for (i in unique(mf.obj$chr.num)) {
        i.index <- i==mf.obj$chr.num
        yp <- if (y.type=="matrix") y[i.index,,drop=FALSE] else y[i.index]
        if (length(yp)==0) 
            print(paste("no data for chr.num = ",i))
        else { 
            x <- mf.obj$pos.plot[i.index]
            x.range <- range(x)
            plot(x.range,y.range,type="n",ylab=ylab,xlab=xlab )
            fun( x, yp, type=type, ...)
            if (include.rug) rug(x[mf.obj$is.marker[i.index]],lwd=rug.lwd)
            title(paste(title.string,"Chromosome", i))
        }
    }
    NULL
    invisible()
}


