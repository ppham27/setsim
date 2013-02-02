bisect <-
function(fn,a,b,tol=0.00001,error=0.0001,decimal=5,input,VZ,c_value){
	fn1 <- function(x){fn(x,...)}
	Ans <- 0
	diff <- 10
	old_p <- (a+b)/2
	if (fn(a,input,VZ,c_value)==0){Ans <- a}
	else if (fn(b,input,VZ,c_value)==0){Ans <- b}
	else{
    while (diff >= tol){
      if ((b-a)/2 < error){Ans <- old_p}
      else if (fn(a,input,VZ,c_value)*fn(old_p,input,VZ,c_value) > 0){a <- old_p}
      else{b <- old_p}
      new_p <- (a + b)/2
      diff <- abs(fn(new_p, input, VZ, c_value) - fn(old_p,input, VZ, c_value))
      old_p <- new_p
      Ans <- old_p}
    round(Ans, decimal)}
}