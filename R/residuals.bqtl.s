"residuals.bqtl"<-
    function(object,...)
{
    coefs <- coef(object)
    coef.names <- names(coefs)[-1]        # drop intercept
    bqtl.call <- match.call(bqtl,object$call)  
    bqtl.terms <- terms(eval(bqtl.call$reg.formula))
    data.obj <- bqtl.call[["ana.obj"]]
    newdata <- substitute(data.obj$data)
    model.fr <-
        model.frame(bqtl.terms,data=eval(newdata),
                    na.action=na.omit)
    model.mat <- model.matrix(delete.response(bqtl.terms),model.fr)
    response <- model.extract(model.fr,"response")
    
    coef.order <-  1 +                   #assume intercept is first
        match(dimnames(model.mat)[[2]][-1],coef.names)
    
    resids <-   response -model.mat %*% as.matrix( coefs[c(1,coef.order)] )
    resids
    
}
