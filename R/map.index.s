"map.index"<-
    function(x,...)
    UseMethod("map.index")
"map.index.default"<-
    function(x,chromo,cM)
{
    if (length(chromo)==1){
        subset <-    x$chr.num==chromo
        if (!missing(cM)){
            this.cM <- x[ subset ,"cM" ]
            loc <-
                seq(along=this.cM)[ abs(this.cM-cM) == min(abs(this.cM-cM)) ][1]
            seq(along=subset)[subset][loc]
        }
        else
        {
            seq(along=subset)[subset]
        }
    }
    else
    {# two locations find a loci in span
        if (length(chromo)>2) stop("> 2 values for chromo not allowed")
        if (chromo[1]>chromo[2]) stop("chromo[1]>chromo[2]")
        subset <-
            ( (x$chr.num <= chromo[2]) & (x$chr.num >= chromo[1] ) )
        if (missing(cM)){
            return( seq(along=subset)[subset] )
        }
        else {
            if (length(cM) != 2) stop("cM must have 2 elements")
            if (chromo[2]==chromo[1] && cM[1]>cM[2]) stop("cM[1]>cM[2]")
            seq(from=map.index(x,chromo[1],cM[1]),
                to=map.index(x,chromo[2],cM[2]))
        }
    }
}
"map.index.analysis.object"<-
    function(x,...)
    map.index(x$map.frame,...)
