

z.score.single.cells <- function ( ma ) {
	i = 0
	ret <- t(
			apply(ma,1, function (x) {
						i = i+1
						n <- which(x==0)
						if ( length(x) - length(n) > 1 ){
							x[-n] <- scale(as.vector(t(x[-n])))
						}
						else {
							x[] = -20
						}
						x[n] <- -20
						x}
			)
	)
	ret
}

