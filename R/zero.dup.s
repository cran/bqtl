"zero.dup"<-
function(x, dig = 6)
{
	ifelse(duplicated(signif(x, dig)), 0, x)
}
