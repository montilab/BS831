#####
#This function converts a color string from colors() to hexidecimal string
#####
#' @export
to.hex<-function(x){
	cols<-col2rgb(x)
	red<-cols[1]
	green<-cols[2]
	blue<-cols[3]
	return(rgb(red, green, blue, maxColorValue = 255))
}
