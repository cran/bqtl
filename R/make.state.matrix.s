"make.state.matrix"<-
    function(marker.frame, marker.distances, method = "F2")
{
###it is assumed
### F2 method has as.numeric(marker.states)%in%c( 1,2,3,4,5,6 )
### for BC1 method as.numeric(marker.states)%in%c( 1,2,5,6 )
### or NA's
###  if NA's are present, they are recoded to 6
  
                                        # no this is wrong NA cause error

###  the assumed setup is as follows:
###    strains are A and a
###
###  type   F2.code   BC.code RI.code
###  AA      1         1       1
###  Aa      2         2 
###  aa      3                 2 
###  A-      4           
###  a-      5           
###  --      6         6       6
###
### marker. frame is assumed to be all numeric
###
    sum.or.one <- function(x)  {y <- sum(x);if (y>0) y else 1.0}
    md.unique <- unique(marker.distances)
    md.ind <- factor(marker.distances, md.unique)
    n.ind <- length(md.unique)
###
### RI strains use an approximation - they aren't quite a Markov Chain
    
    if (method=="RI.self") md.unique <- md.unique/(2-md.unique)
    if (method=="RI.sib") md.unique <- md.unique/(4-3*md.unique)
    
    
    mf.num <- as.matrix(marker.frame)
    if (any(is.na(mf.num)))
        stop("marker.frame must be all numeric")

    if (n.ind==0){            # only one marker
      if (method=="F2")
        return(
               array(rbind(diag(3),c(1,2,0)/3,c(0,2,1)/3,c(1,2,1)/4)[mf.num,],
                     c(dim(mf.num),3)))
      else
        return(array(rbind(diag(2),NA,NA,NA,c(1,1)/2)[mf.num,],
                     c(dim(mf.num),2)))
    }

    switch(method,
           F2 = {
               tmat <- matrix(c(1, 1, 1, -1, 0, 1, 1, -1, 1), 3, byrow
                              = F)
               tinv <- matrix(c(1, -2, 1, 2, 0, -2, 1, 2, 1)/4, 3, 
                              byrow = F)
               t1 <- tmat[, 1, drop = F] %*% tinv[1,  , drop = F]
               t2 <- tmat[, 2, drop = F] %*% tinv[2,  , drop = F]
               t3 <- tmat[, 3, drop = F] %*% tinv[3,  , drop = F]
               tx.mat <- array(rep(t1, n.ind), c(3, 3, n.ind)) +
                   outer(t2, md.unique) + outer(t3, md.unique^2)
               unc.pr <- t1[1,  ]
               dom.A.pr <- c(1,2,0)/3
               dom.a.pr <- c(0,2,1)/3
               lp <- 3
           }
           ,
           RI.self=,
           RI.sib=,
           BC1 = {
               t1 <- matrix(1/2, nc = 2, nr = 2)
               t2 <- matrix(c(1, -1, -1, 1)/2, nc = 2)
               tx.mat <- array(rep(t1, n.ind), c(2, 2, n.ind)) +
                   outer(t2, md.unique)
               unc.pr <- t1[1,  ]                        
               lp <- 2
           }
           ,
           stop(paste(method, " method unknown")))
    left.marker <-
        matrix(c(diag(lp),rep(NA, if (lp==2) 8 else 9)),nr=lp)[,mf.num]
    dmm <- c(lp, dim(marker.frame))
    dim(left.marker) <- dmm
    center.marker <- right.marker <- left.marker
    dimnames(center.marker) <- list(NULL, NULL, dimnames(marker.frame)[[2]])
    zero.dx <- marker.distances == 0
    chr.number <- sum(zero.dx)
    is.miss <- mf.num[,c(TRUE, zero.dx)]==6
    is.miss <- rep(is.miss,rep(lp,length(is.miss)))
    left.marker[,  , c(TRUE, zero.dx)][is.miss] <- unc.pr
    is.miss <- mf.num[,c(zero.dx, TRUE)]==6
    is.miss <- rep(is.miss,rep(lp,length(is.miss)))
    right.marker[,  , c(zero.dx, TRUE)][is.miss] <- unc.pr
    if (method=="F2") {

        ## initialize first marker if necessary
        is.A.minus <- mf.num[,c(TRUE, zero.dx)]==4
        is.A.minus <- rep(is.A.minus,rep(3,length(is.A.minus)))
        is.a.minus <- mf.num[,c(TRUE, zero.dx)]==5
        is.a.minus <- rep(is.a.minus,rep(3,length(is.a.minus)))
        left.marker[,  , c(TRUE, zero.dx)][is.A.minus] <- dom.A.pr
        left.marker[,  , c(TRUE, zero.dx)][is.a.minus] <- dom.a.pr
        ## initialize last marker if necessary
        is.A.minus <-mf.num[,c(zero.dx, TRUE)]==4
        is.A.minus <- rep(is.A.minus,rep(3,length(is.A.minus)))
        is.a.minus <-mf.num[,c(zero.dx, TRUE)]==5
        is.a.minus <- rep(is.a.minus,rep(3,length(is.a.minus)))
        right.marker[,  , c(zero.dx, TRUE)][is.A.minus] <- dom.A.pr
        right.marker[,  , c(zero.dx, TRUE)][is.a.minus] <- dom.a.pr
    }
    for(i in 2:dmm[3]) {
        i.tx <- as.numeric(md.ind[i - 1])
        is.miss <- mf.num[,i]==6
        if(any(is.miss)) {
            should.be <- t(as.matrix(left.marker[, is.miss, i - 1])
                           ) %*% tx.mat[,  , i.tx]
            left.marker[, is.miss, i] <- t(should.be)
        }
        if (method == "F2") {
            is.A.minus <- mf.num[,i]==4
            is.a.minus <- mf.num[,i]==5
            if(any(is.A.minus)) {
                
                should.be <- t(t(as.matrix(left.marker[, is.A.minus, i - 1])
                                 ) %*% (tx.mat[,  , i.tx]*
                                        rep(dom.A.pr/unc.pr,rep(3,3))))

                left.marker[, is.A.minus, i ] <-
                    sweep(should.be, 2, apply(should.be, 2, sum.or.one), "/")
            }
            if (any(is.a.minus)) {
                should.be <- t(t(as.matrix(left.marker[, is.a.minus, i - 1])
                                 ) %*% (tx.mat[,  , i.tx]*
                                        rep(dom.a.pr/unc.pr,rep(3,3)))/unc.pr)

                left.marker[, is.a.minus, i ] <-
                    sweep(should.be, 2, apply(should.be, 2, sum.or.one), "/")
            }
        }
    }
    for(i in rev(1:(dmm[3] - 1))) {
        i.tx <- as.numeric(md.ind[i])
        is.miss <- mf.num[,i]==6
        if(any(is.miss)) {
            should.be <- t(as.matrix(right.marker[, is.miss, i + 1]
                                     )) %*% tx.mat[,  , i.tx]
            right.marker[, is.miss, i] <- t(should.be)
        }
        if (method == "F2") {
            is.A.minus <- mf.num[,i]==4
            is.a.minus <- mf.num[,i]==5
            if(any(is.A.minus)) {
                
                should.be <- t(t(as.matrix(right.marker[, is.A.minus, i + 1])
                                 ) %*% (tx.mat[,  , i.tx]*
                                        rep(dom.A.pr/unc.pr,rep(3,3))))

                right.marker[, is.A.minus, i ] <-
                    sweep(should.be, 2, apply(should.be, 2, sum.or.one), "/")
            }
            if (any(is.a.minus)) {
                should.be <- t(t(as.matrix(right.marker[, is.a.minus, i + 1])
                                 ) %*% (tx.mat[,  , i.tx]*
                                        rep(dom.a.pr/unc.pr,rep(3,3))))

                right.marker[, is.a.minus, i ] <-
                    sweep(should.be, 2, apply(should.be, 2, sum.or.one), "/")
            }
        }
    }
    for(i in 1:dmm[3]) {
        is.miss <- mf.num[,i]==6
        if(any(is.miss)) {
            should.be <- as.matrix((left.marker[, is.miss, i] * 
                                    right.marker[, is.miss, i])/unc.pr)
            
            center.marker[, is.miss, i] <-
                sweep(should.be, 2, 
                      apply(should.be, 2, sum.or.one), "/")
        }
        if (method == "F2") {
            is.A.minus <- mf.num[,i]==4
            is.a.minus <- mf.num[,i]==5
            if(any(is.A.minus)) {
                mpy <- ifelse(dom.A.pr==0.0,0.0,(1/dom.A.pr/unc.pr^2))
                should.be <- as.matrix((left.marker[, is.A.minus, i] * 
                                        right.marker[, is.A.minus, i])*mpy)
                

                center.marker[, is.A.minus, i ] <-
                    sweep(should.be, 2, apply(should.be, 2, sum.or.one), "/")
            }
            if (any(is.a.minus)) {
                mpy <- ifelse(dom.a.pr==0.0,0.0,(1/dom.a.pr/unc.pr^2))
                should.be <- as.matrix((left.marker[, is.a.minus, i] * 
                                        right.marker[, is.a.minus, i])*mpy)
                

                center.marker[, is.a.minus, i ] <-
                    sweep(should.be, 2, apply(should.be, 2, sum.or.one), "/")
            }
        }
    }

    aperm(center.marker, c(2, 3, 1))
}
