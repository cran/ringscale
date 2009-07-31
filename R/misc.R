gauss.pdf <- function(mu = 0, sigma = 1) {
	# returns the probability density function for
	# the normal distribution with parameters 'mu'
	# and 'sigma'
	
	function(x) {
		1/(sigma*sqrt(2*pi)) * exp(-0.5*((x-mu)/(sigma))^2)
	}
}

mask.circle <- function(width = 11, height = width, radius.outer = 10, radius.inner = 0, x0 = floor((width+1)/2), y0 = floor((height+1)/2)) {
	# creates a circular Boolean mask embedded
	# in a rectangular mask
	
	x <- (1:width) - x0
	y <- (1:height) - y0
	
	mask <- sqrt(t(matrix(x, width, height))^2 + (matrix(y, height, width))^2)
	mask <- (mask < radius.outer) & (mask >= radius.inner)
}

tcltk.path.get <- function() {
	# tcl/tk choose path dialogue
	
	if (interactive()){
		require(tcltk)
		tclvalue(tkchooseDirectory())
	}
}

tcltk.file.get <- function() {
	# tcl/tk open file dialogue
	
	if (interactive()){
		require(tcltk)
		tclvalue(tkgetOpenFile(filetypes = "{{Ringscale Sessions} {.RData}} {{All files} *}"))
	}
}

graphical <- function() {
	# true, if X11 and tclktk are available
	all(interactive(), capabilities("X11"), capabilities("tcltk"))
}

progress.Bar <- function(percent, bar.name = "", bar.char = "-", bar.length = 32, bar.left = "|", bar.right = bar.left, bar.end = "") {
	# draws a progress bar for the console
	
	if (interactive()) {
		cat(bar.name, bar.left, rep(bar.char, round(bar.length*percent)), rep(" ", bar.length - round(bar.length*percent)), bar.right, " ", fill(floor(100*percent), to = 3, with = " "), "%", bar.end, "\r", sep = "")
		flush.console()
	}
}

progress.Meter <- function(count, symbols = c("-", "\\", "|", "/", "-", "\\", "|", "/"), length = 1, left = " ", right = " ") {
	# creates an "animated" progress meter
	
	paste(left, symbols[((count %% length(symbols)) + 1):(min(length(symbols),(count %% length(symbols)) + length))], right, sep = "")
}

strlen <- function(x) {
	# return the length of a string
	
	length(unlist(strsplit(as.character(x), split = "")))
}

fill <- function(x, to = 4, with = "0") {
	# fills the string with 'with' from the
	# left to a given length 'to'
	# example: fill("17", to = 4, with = "0") is "0017"
	
	x.out <- NULL
	length(x.out) <- length(x)
	for (i in 1:length(x)) {
		remaining <- to - strlen(x[i])
		filled <- as.character(x[i])
		fills <- paste(rep(with, max(0,remaining)), collapse = "")
		x.out[i] <-paste(fills, filled, sep = "")
	}
	x.out
}

minsec <- function(t) {
	# converts seconds into a string of minutes and seconds
	# example: minsec(1704) is "28m24s"
	
	t.min <- floor(t/60)
	t.sec <- floor(t - 60*t.min)
	return(paste(t.min, "m", fill(t.sec, 2), "s", sep = ""))
}

instrument.get <- function(header, key, name, symbol, unit) {
	# needed to read in instrument specific values, such as
	# diameter, observation wavelength, ...
	
	# try to find it in the FITS header first
	value <- as.numeric(header.find(header, key))
	
	# if we didn't find it in the header, ask the user
	if ((length(value) == 0 || is.na(value)) && interactive()) {
		repeat {
			input <- as.numeric(readline(paste("Could not determine ", name, ", please enter. ", symbol, " [", unit, "] = ", sep = "")))
			if (!is.na(input)) {
				value <- input
				break
			}
		}
	} else
		cat(symbol, " = ", value, " ", unit, "\n", sep = "")
	
	value
}

header.find <- function(header, key) {
	# find the value for the key 'key' in the
	# specified FITS header
	
	header[2*which(header[((1:length(header)) %% 2) == 1] == key)]
}

speckle.Diameter <- function(D = 8.2, lambda = 2.2, pixscale = 0.013270) {
	# returns the speckle diameter [px]
	# D: telescope diameter [m]
	# lambda: wavelength [micro m]
	# pixscale: pixel scale [arcsec/px]
	
	# first zero of J_1 (Bessel function) divided by pi
	const.Bessel <- 1.219670
	
	# pixel scale [rad]
	pixscale.rad <- pixscale * 4.84813681109535e-6
	
	# D.speckle = 2 * theta_1 = 2 * 1.22 * lambda/D
	2 * const.Bessel * lambda * 1e-6 / D / (pixscale.rad)
}