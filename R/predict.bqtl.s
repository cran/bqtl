"predict.bqtl"<-
    function(object,newdata)
{
    coefs <- coef(object)
    coef.names <- names(coefs)[-1]        # drop intercept
    bqtl.call <- match.call(bqtl,object$call)
    reg.form <- eval(bqtl.call$reg.formula)
###    bqtl.terms <- terms(bqtl.call$reg.formula)
    if (missing(newdata)) {
        data.obj <- bqtl.call[["ana.obj"]]
        newdata <- substitute(data.obj$data)
        bqtl.terms <- terms(reg.form)
        model.fr <-
            model.frame(bqtl.terms,data=eval(newdata),
                        na.action=na.omit)
    }
    else {
        bqtl.terms <- delete.response(terms(reg.form))
        if (class(newdata) == "analysis.object")
            model.fr <- model.frame(bqtl.terms,data=eval(newdata$data),
                                    na.action=na.omit)
        else
            model.fr <- model.frame(bqtl.terms,data=eval(newdata),
                                    na.action=na.omit)
    }
    
    model.mat <- model.matrix(delete.response(bqtl.terms),model.fr)
    
    coef.order <-  1 +                   #assume intercept is first
        match(dimnames(model.mat)[[2]][-1],coef.names)

    fit <-  model.mat %*% as.matrix( coefs[c(1,coef.order)] )
    dimnames(fit)[[2]] <- "fitted.value"
    fit
}
"fitted.bqtl"<-
    predict.bqtl
