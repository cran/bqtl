"make.marker.numeric" <- function(marker.frame,
                                  level.names=NULL )
{
### convert factor or character frames into numeric
### for use in make.state.matrix
  if (is.null(level.names)) {
    level.names <- c("AA","Aa","aa","A-","a-","--")}
  else
    if ( length(level.names) != 6 )
      stop("length(level.names) must = 6 if none NULL value is used")

  char.frame <- unlist(lapply(marker.frame,as.character))
  char.frame[is.na(char.frame)] <- level.names[6]
  char.frame[char.frame=="NA"] <- level.names[6]
  level.numeric <- seq(length(level.names))
  names(level.numeric) <- level.names
  dm <- dim(marker.frame)

  res <- matrix(level.numeric[char.frame],ncol=dm[2],
                dimnames=list(NULL,names(marker.frame)))
  if (any(is.na(res))) {
    stop(paste("unrecognized genotype(s):",unique(res[is.na(res)])))
  }
  res
  
}
