"make.analysis.obj"<-
    function(data, map.frame, marker.frame, marker.levels=NULL,method="F2",
             casewt=NULL,varcov=FALSE,mode.mat=NULL)
{
    
### map.frame - a map.frame.object OR a
### lambda - vector of pseudo-marker distance parms - has names attribute
### marker.frame has dimnames attribute to match 
###
### methods are: F2, BC1, RI.self, and RI.sib
    "check.markers" <-
        function(x,y)
        {
            tab <- table(unlist(sapply(x,as.character)))
            bad.x <- !(names(tab) %in% c(y,"NA"))
            bad.y<- !(y %in%names(tab))
            list(bad.markers=tab[bad.x],unused.levels=y[bad.y])
        }
    
    this.call <- match.call()
    
    if( inherits(data,"map.frame") )
        warning("first arg should not be a map.frame! ?? ")
    
    if( inherits(marker.frame,"map.frame") )
        warning("third arg should not be a map.frame! ?? ")
    
    method.types <- c("F2","BC1","RI.self","RI.sib")
    method.code <- pmatch(method,method.types)
    if (is.na(method.code))
        stop(paste("method =",method, "not found"))
    method <- method.types[method.code]
    
    if (length(dm.data <- dim(data))==0) dm.data<- c(length(data),1)
    if (dm.data[1]!=nrow(marker.frame)) stop(
               "nrow(data)!=nrow(marker.frame)")
    
    
### RI strains use coding of 1,2,6
    
    if (length(marker.levels)==0) marker.levels <-
        switch(method,
               F2=f2.levels(),
               BC1=bc1.levels(),
               RI.self=,
               RI.sib=ri.levels() 
               )
### set mode.mat
    if (is.null(mode.mat)){
        mode.mat <-
            if (method=="F2")
                matrix(c(-1,0,1,-1,1,-1),nc=2,
                       dimnames=list(marker.levels[1:3],c("add","dom")))
            else
                matrix(c(-1,1),nc=1,dimnames=list(marker.levels[1:2],NULL))
    }
    else {   # check the supplied version
        dmmt <- dim(mode.mat)
        if (method=="F2") {
            if ( length(dmmt)==0 || dmmt[1]!=3 || dmmt[2] > 2)
                stop("F2 mode.mat must be 3 x 2 matrix")
            if (dmmt[2]==1)
                warning("F2 mode.mat has just one column")
            else
                if (length(dimnames(mode.mat)[[2]])==0)
                    stop("need column names for F2 mode.mat")
        }
        else {
            if ( length(dmmt)==0 || dmmt[1]!=2 || dmmt[2] > 1)
                stop("non-F2 mode.mat must be 2 x 1 matrix")
        }
    }
    
    if(inherits(map.frame,"map.frame")) {
        lambda <- map.frame$lambda
        marker.names <-
            as.character(map.frame$marker.name[is.element(map.frame$is.marker,
                                                          c("TRUE","T"))])
        names.lambda <- as.character(map.frame$marker.name)
        m.chromo <- map.frame$chr.num
    }
    else {
        ## not actually a map.frame
        names.lambda <- names(map.frame)
        if (length(map.frame) != length(names.lambda))
            stop("problem with map.frame arg - must be data frame or vector with names()")
        lambda <- map.frame
        marker.names <- dimnames(marker.frame)[[2]]
        m.chromo <- cumsum(c(1,lambda==0))[-length(lambda)-1]
    }
    
### check the setup
### markers found to match marker.names in map.frame
    if (any(lost.markers <-
            !is.element(marker.names,dimnames(marker.frame)[[2]])
            )){
        cat(deparse(this.call$map.frame),"defines markers NOT found in",
            deparse(this.call$marker.frame),":\n")
        print(marker.names[lost.markers])
        stop("bailing out")
    }
    
### marker values match marker.levels??

    ck.mark <- check.markers(marker.frame[,marker.names],marker.levels)
    if (length(ck.mark$bad) != 0)
    {
        cat("Marker values found, but not defined in marker.levels:\n")
        print(ck.mark$bad)
        stop("bailing out")
    }
    if (any(marker.levels[1:2]%in%ck.mark$unused.levels))
        warning(paste("unused marker.levels:",
                      paste(ck.mark$unused.levels,collapse=" ")))
    if (method=="F2" && all(marker.levels[3:5]%in%ck.mark$unused.levels) )
        warning(paste("method=F2, but did not find",
                      paste(marker.levels[3:5],collapse=", "),
                      "in",deparse(this.call$marker.frame),"\n"))
    m.state.matrix <-
        array(0,c(dm.data[1],length(lambda),if (method=="F2") 3 else 2),dimnames=list(NULL,names.lambda,NULL))
    m.loc.right <- array(0,c(dm.data[1],length(lambda)),dimnames=list(NULL,names.lambda))
    
    
### the principal reason for the explicit loop is to avoid large objects in R
###  --- the space required sometimes kills the whole function
    for (i in unique(m.chromo)) {
        which.chromo <- m.chromo==i
        m.frame <- as.data.frame(matrix(NA, nr = dm.data[1], nc = sum(which.chromo)))
        names(m.frame) <- names.lambda[m.chromo==i]
        names.to.use <-
            names(m.frame)[is.element(names(m.frame),marker.names)]
        
        if (length(names.to.use) == 0) stop(
                  "names(lambda) did not match names(marker.frame)")
        for (i.name in names.to.use) m.frame[,i.name] <- marker.frame[,i.name]
        
        lambda.i <- lambda[which.chromo][ - sum(which.chromo)]
        m.state.matrix[,which.chromo,] <-
            make.state.matrix( make.marker.numeric(m.frame,marker.levels),
                              lambda.i , method)
        left.loc <- min(seq(which.chromo)[which.chromo])-1
        m.loc.right[,which.chromo] <- left.loc + make.loc.right(m.frame, lambda.i)
    }
    
    m.regressor.matrix <- make.regressor.matrix(m.state.matrix,mode.mat)
    
    m.varcov <-
        if (varcov!=FALSE)
            make.varcov(m.regressor.matrix, data[,varcov],casewt=casewt)
        else
            NULL
    reg.names <- dimnames(m.regressor.matrix)[[2]]
### sanity check on m.regressor.matrix
    if (any(null.var <-
            apply(apply(m.regressor.matrix,2,range),2,diff) <
            0.5/nrow(m.regressor.matrix)))
        warning(paste("null or tiny variance in",
                      paste(reg.names[null.var],collapse=" "),"\n"))
    
    dim(reg.names) <- c(dim(m.state.matrix)[[3]]-1,dim(m.state.matrix)[[2]])
    dimnames(reg.names) <-
        list(
             if (method!="F2") NULL else dimnames(mode.mat)[[2]],
             names.lambda
             )
    version.stamp <- if (exists("is.R")&&is.R()) version.bqtl else "splus"
    structure(
              list(data=data.frame(data,m.regressor.matrix), varcov = m.varcov, 
                   reg.names=reg.names,method=method,
                   state.matrix = m.state.matrix, loc.right = m.loc.right,
                   map.frame=map.frame,casewt=casewt,mode.mat=mode.mat,
                   call = this.call,version=version.stamp),
              class="analysis.object"
              )
}
