"update.bqtl"<-
    function(setup,loc.mat,ana.obj)
{
    dmloc <- dim(loc.mat)
    if (length(dmloc) == 0){
        nloc <-   length(loc.mat)
        n.alt <- 1
    }
    else {
        nloc <- dmloc[1]
        n.alt <- dmloc[2]
    }
    if ( nloc != setup$nparm[4] )
        stop("nrow(loc.mat) not = number of loci in setup")
    if (any( loc.mat > ncol(ana.obj$loc.right))||any(loc.mat<1) )
        stop("invalid locus number")
    subset <- attr(setup,"subset")
    lambda <- switch(ana.obj$method,
                   RI.self = {ana.obj$map.frame$lambda/(2 - ana.obj$map.frame$lambda)},
                   RI.sib = {ana.obj$map.frame$lambda/(4 - 3 * ana.obj$map.frame$lambda)},
                   ana.obj$map.frame$lambda)
    setup <-
        c(list("upbqtl"),
          setup[-c(1,22)],
          list( loc.mat = as.integer(loc.mat-1)) ,
          list( n.alt = as.integer(n.alt)) ,
          list( loc.order = as.integer(apply(as.matrix(loc.mat),2,order)-1)) ,
          list( res = double(n.alt)),
          list( orig.x = as.double(setup$x)) ,
          list( loc.right = as.integer(ana.obj$loc.right[subset,]-1)) ,
          list( map.lambda = as.double(lambda)) ,
          list( state.matrix = as.double(ana.obj$state.matrix[subset,,])) ,
          list( n.state.loc = as.integer(dim(ana.obj$state.matrix)[2])))
    
    res <- do.call(".C",setup)
    
    res$res
}

