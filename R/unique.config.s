"uniq.config"<-
    function(swap.obj)
{
    configs <- apply(swap.obj$config, 2:3, function(x)
                     paste(sort(x), collapse = "."))
    unique.con <- unique(configs)
    conf.match <- match(configs, unique.con)
    uniq.mat <-
        matrix(swap.obj$config, nr =
               dim(swap.obj$config)[1])[, ! duplicated(conf.match), drop = FALSE]
    list(uniq = uniq.mat, match = conf.match)
}

