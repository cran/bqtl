"residuals.linear.bayes" <-
    function(object,...)
    predict.linear.bayes(object,...,return.resids=TRUE)

### since most of the code is the same in residuals and predict, why
### duplicate it?

"fitted.linear.bayes" <-
    function(...)
    predict.linear.bayes(...)

"predict.linear.bayes"<-
    function(object,newdata=lb.call$ana.obj,return.resids=FALSE,...)
{
    coefs <- coef(object)
    lb.call <- match.call(linear.bayes,object$call)
    x <- eval(lb.call$x)
    partial <- eval(lb.call$partial)
    if (x[[1]]!="~"){          #not a formula:try to go on...
        if (return.resids)
            stop("cannot create residuals without formula in linear.bayes()")
        newdata <- eval(substitute(newdata))
        preds <- as.matrix(newdata)%*%coefs
        preds <- preds - mean(preds)    #center the results
        return(preds)
    }
    ana.obj <- eval(newdata)
    scope <- ana.obj$reg.names
    
### take this from varcov:
    bqtl.specials <- c("configs","locus") # dont use: add or dom so
                                        # swap, hk wont jam up    
    reg.terms <- terms(x,specials=bqtl.specials)
    if (any(attr(reg.terms,"order")>1))
        stop("cannot use interactions in this formula")
    covs <-
        rhs.bqtl(reg.terms,ana.obj,bqtl.specials,NULL,scope,
                 expand.specials=FALSE)
    response <-
        dimnames(attr(reg.terms,"factors"))[[1]][attr(reg.terms,"response")]
    if (!is.null(partial)) {
        p.terms <- terms(partial,specials=bqtl.specials)
        if (any(attr(p.terms,"order")>1))
            stop("cannot use interactions in this formula")
        cntl <-
            rhs.bqtl(p.terms,ana.obj,bqtl.specials,NULL,scope,
                     expand.specials=FALSE)
        dummy.formula <- eval(parse(text=paste(response,"~",covs,"+",cntl)))
    }
    else {
        cntl <- NULL
        dummy.formula <- eval(parse(text=paste(response,"~",covs)))
    }
    dummy.terms <- terms(dummy.formula)
    mf <- model.frame(dummy.terms,ana.obj$data, na.action=na.omit)
    if (is.null(partial)) {
        x <- model.matrix(dummy.terms,mf)[,-1]
        y <- model.extract(mf,"response")
        preds  <- x%*%matrix(coefs,nc=1)
        ## can't always get an intercept from linear.bayes, hence:
        preds<-preds + mean(y) - mean(preds)
    }
    else
    {
        
        cntl.out <- eval(parse(text=paste("~ . -", cntl)))
        x.formula <- update(dummy.formula,
                            eval(parse(text=paste("~ . -", cntl)))  )
        x <- model.matrix(terms(x.formula),mf)[,-1]
        y <- model.extract(mf,"response")
        z <- model.matrix(terms(eval(parse(text=paste("~", cntl)))),
                          mf)[,-1,drop=FALSE]
        pred.z <- y - lsfit(z,y)$resid
        preds <- x%*%matrix(coefs,nc=1) + pred.z
        preds <- preds + mean(y) - mean(preds)
    }
    if (return.resids){
        y - preds
    }
    else
        preds
}

