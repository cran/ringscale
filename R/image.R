projections.argmax <- function(img) {
	# give the argmax of the projection of an
	# image onto the coordinate axes
	
	stopifnot(is.matrix(img))
	c(which.max(apply(img, 1, sum, na.rm = TRUE)), which.max(apply(img, 2, sum, na.rm = TRUE)))
}

img.centroid <- function(img) {
	# calculates the centroid of an image
	
	img.dim <- dim(img)
	X <- img.dim[1]
	Y <- img.dim[2]
	
	weighted.sum <- c(0,0)
	
	for (x in 1:X)
		for (y in 1:Y)
			weighted.sum <- weighted.sum + c(x,y) * img[x,y]
	
	weighted.sum/sum(img, na.rm = TRUE)
}

star.locate <- function(img, radius = 5, method.approx = projections.argmax, method.exact = img.centroid) {
	# locate the star in an image, using 'method.approx'
	# for approximate location and 'method.exact' for a
	# more exact estimation of a small region around the
	# approximate location
	
	pos.approx <- method.approx(img)
	
	# cut out region around the approximate location...
	region <- img[(pos.approx[1]-radius):(pos.approx[1]+radius), (pos.approx[2]-radius):(pos.approx[2]+radius)]
	# ... and pass it to the more exact, but local method
	pos.approx - (radius + 1) + method.exact(region)
}

circle.add <- function(pos.rel, radius.rel = 0.005, line.width = 1, color = "green") {
	# add a circle to the current image;
	# position and radius are relative to
	# the image size (i.e. in [0,1])
	
	op <- par(lwd = line.width)
	symbols(x = pos.rel[1], y = pos.rel[2], circles = radius.rel, fg = color, inches = FALSE, add=TRUE)
	par(op)
}

img.astro <- function(img, pixscale = NA, delta.arcsec = 0.1, colormap = "heat", distance.pc = NA, name = NULL) {
	# for the visualization of the results; has some useful features for
	# astronomical images (distances in arcsec and AU, for example)
	
	img.dim <- dim(img)
	X <- img.dim[1]
	Y <- img.dim[2]
	
	# color maps
	if 		(colormap == "bw") { colors <- c("black", "white") }
	else if	(colormap == "wb") { colors <- c("white", "black") }
	else if	(colormap == "heat") { colors <- c("black", "orangered", "orange", "yellow", "white") }
	else if (colormap == "bbryw") { colors <- c("black", "blue", "red", "yellow", "white") }
	
	image(img, xlab = "x [px]", ylab = "y [px]", axes = FALSE, col = colorRampPalette(colors)(100))
	if (!is.null(name))
		title(name)
	
	# label axes with pixel coordinates
	k <- ceiling(X/10)
	axis(1, rev(seq(from = X - 1, to = 0, by = -k))/(X - 1), rev(seq(from = X, to = 1, by = -k)))
	axis(2, rev(seq(from = Y - 1, to = 0, by = -k))/(Y - 1), rev(seq(from = Y, to = 1, by = -k)))
	
	# draw the scale [arcsec] in the bottom left corner
	if (!is.na(pixscale) && delta.arcsec > 0) {
		n <- delta.arcsec / pixscale
		k <- ceiling(X/6/n)
		lines(c(0.05, 0.05 + k*n/X), c(0.07, 0.07), col = "green", lwd = 2)
		text(0.13, 0.10, paste(k*n*pixscale, " arcsec", sep = ""), col = "green")
		# also draw the scale [AU], of distance of object is available
		if (!is.null(distance.pc) && !is.na(distance.pc))
			text(0.13, 0.04, paste(round(k*n*pixscale*distance.pc,1), " AU", sep = ""), col = "green")
	}
}

filter.FFT <- function(x,fir,a=1,b=0) {
	# filter via Fast Fourier Transform
	
	g <- matrix (0, nrow=nrow(x), ncol=ncol(x))
	N <- nrow(x)
	M <- ncol(x)
	n <- nrow(fir)
	m <- ncol(fir)
	n2 <- n%/%2
	m2 <- m%/%2
	g [1:(n-n2), 1:(m-m2)] <- fir [(n2+1):n, (m2+1):m]
	g [(N-n2+1):N, 1:(m-m2)] <- fir [1:n2, (m2+1):m]
	g [1:(n-n2), (M-m2+1):M] <- fir [(n2+1):n, 1:m2]
	g [(N-n2+1):N, (M-m2+1):M] <- fir [1:n2, 1:m2]
	z <- fft (fft(x) * fft(g), inverse=TRUE)
	y <- Re(z) / length(x)
	return (a*y+b)
}

filter.gauss.1D <- function(radius = 5, sigma = radius) {
	# creates a one.dimensional Gaussian filter mask
	# of size 2*radius + 1
	
	(gauss.pdf(mu = 0, sigma = sigma))((-radius:radius)/sigma)
}

img.lowpass <- function(img, sigma = 1, radius = 3) {
	# lowpass filter an image with a two-dimensional
	# Gaussian filter mask
	
	# create a 2D Gaussian mask out of an 1D mask
	mask.gauss <- filter.gauss.1D(radius, sigma) %*% t(as.matrix(filter.gauss.1D(radius, sigma)))
	
	# filter via Fast Fourier Transform
	filter.FFT(img, mask.gauss)
}

clip <- function(img, min = 0, max = 100) {
	# clips an image according to the given
	# min and max values
	
	img[img < min] <- min
	img[img > max] <- max
	img
}

img.scale <- function(img, percent = 90) {
	# automatically clips the image, with min and max
	# values based on the 'percent' parameter
	
	q <- quantile(x = img, probs = c(100 - percent, percent + 100)/2/100, na.rm = TRUE, type = 3)
	gray.lower <- q[1]
	gray.upper <- q[2]
	
	list(image = clip(img, gray.lower, gray.upper), lower = gray.lower, upper = gray.upper)
}

roi <- function(img, radius) {
	# given a single image 'img', this function creates
	# a region of interest (ROI) image, wherein the star
	# is centered
	
	stopifnot(is.matrix(img))

	# locate the star
	pos <- star.locate(img)
	
	# we take integer values here, because for subpixel
	# accuracy we would have to resample the image
	pos.round <- round(pos)
	
	img.dim <- dim(img)
	X <- img.dim[1]
	Y <- img.dim[2]
	x <- pos.round[1]
	y <- pos.round[2]
	
	# this is where the final ROI will go
	roi <- array(data = NA, dim = c(2*radius+1,2*radius+1))
	
	# this is what we cut out from the image (be aware,
	# that this may be smaller than the final ROI, because
	# of the image borders)
	roi.cut <- img[max(1, x-radius):min(X, x+radius), max(1, y-radius):min(Y, y+radius)]
	
	roi.dim <- dim(roi.cut)
	U <- roi.dim[1]
	V <- roi.dim[2]
	u <- x - max(1, x-radius)
	v <- y - max(1, y-radius)
	
	# insert the cut-out region into the ROI-skeleton
	roi[(radius + 1 - u):(radius + U - u), (radius + 1 - v):(radius + V - v)] <- roi.cut
	
	roi
}

rois.fromFITS <- function(filename, radius) {
	# given a FITS filename, this function reads out the cube
	# and returns the region of interest (ROI) image(s)
	
	stopifnot(file.exists(filename))
	
	# one FITS cube may contain more than one image
	count.image <- 0
	
	# read out FITS image data
	cube <- readFITS(filename)$imDat
	
	# for FITS cubes with 2 dimensions (i.e. one image per cube)
	if (is.na(dim(cube)[3]))
		dim(cube)[3] <- 1
	
	Z <- dim(cube)[3]
	
	# this is the skeleton for the final ROI images
	rois <- array(dim = c(2*radius + 1, 2*radius + 1, Z))
	
	# cut out the region of interest from the original image
	for (img in 1:Z)
		rois[,,(count.image <- count.image + 1)] <- roi(cube[,,img], radius = radius)
	
	rois
}