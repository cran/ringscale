\name{ringscale}
\alias{ringscale}
\alias{analyze}
\alias{evaluate}

\alias{circle.add}
\alias{clip}
\alias{fill}
\alias{filter.FFT}
\alias{filter.gauss.1D}
\alias{gauss.pdf}
\alias{graphical}
\alias{header.find}
\alias{img.astro}
\alias{img.centroid}
\alias{img.lowpass}
\alias{img.scale}
\alias{instrument.get}
\alias{mask.circle}
\alias{minsec}
\alias{progress.Bar}
\alias{progress.Meter}
\alias{projections.argmax}
\alias{ringscale.batch}
\alias{roi}
\alias{roi.candidates}
\alias{rois.fromFITS}
\alias{speckle.Diameter}
\alias{star.locate}
\alias{strlen}
\alias{tcltk.file.get}
\alias{tcltk.path.get}

\title{Ringscale Method For The Detection Of Faint Companions}
\description{Implementation of the \dQuote{Ringscale} method as proposed
	in the student research project \dQuote{DETECTION OF FAINT COMPANIONS AROUND
	YOUNG STARS IN SPECKLE PATTERNS OF VLT/NACO CUBE MODE IMAGES BY MEANS OF
	POST-PROCESSING} at the Friedrich-Schiller-University of Jena.}
\usage{
analyze(filenames = character(0),
	result.filename = paste(format(Sys.time(),
	"\%Y-\%m-\%d-\%H\%M\%S"), ".RData", sep = ""),
	evaluate = TRUE,
	header.skip = FALSE)
evaluate(filename = character(0), support.min = 0.5)
}
\arguments{
	\item{filenames}{Filenames of the FITS files to be analyzed.}
	\item{result.filename}{Filename of the result of \code{analyze()}.}
	\item{evaluate}{If \code{TRUE}, \code{evaluate()} is called for the appropriate file just after the
	analysis finished.}
	\item{header.skip}{If \code{TRUE}, the FITS-header will not be read out at the beginning of the analysis. You may try this if you encounter performance problems or if you know that your header contains no information which is needed by Ringscale (such as telescope diameter, wavelength et cetera). You will have to enter these by hand then.}
	\item{filename}{Filename of the result of \code{analyze()} to be evaluated.}
	\item{support.min}{Minimum support value for which candidates are listed and marked in \code{evaluate()}.}
}
\details{The \pkg{ringscale} package basically provides two \dQuote{public}
	functions: \code{analyze()} and \code{evaluate()}.
	
	\code{analyze()} applies the \dQuote{Ringscale} method as proposed in the
	paper of the student research project on FITS files. Note, that these FITS
	files do not necessarily have to be pre-processed, so raw VLT/NACO images
	are applicable as input. If \code{analyze()} is called without parameters,
	a \pkg{tcltk}-dialog will open and let the user select the input files.
	Otherwise, the filenames of the processed FITS files have to be given as
	character vector, including the entire path (see examples). The results of
	the analysis are saved in the file \code{result.filename}. This allows to
	have a look at the results at a later time again very quickly, as the
	analysis itself may take some time. If \code{evaluate = TRUE},
	\code{evaluate()} is called for the appropriate file just after the
	analysis finished.
	
	The \code{evaluate()} function evaluates the results provided by
	\code{analyze()}. Therefore, the input consists of only one file,
	the result of \code{analyze()}. The filename can be given as
	parameter. If omitted, a \pkg{tcltk}-dialog will let the user select it
	just like for \code{analyze()}. \code{evaluate()} shows the found
	companions in a hypothesis-map, where the possible companions are marked
	with circles according to their rank and support value (the lower the rank,
	the bigger the radius and the bigger the support, the bigger the line width
	of the circles). Additionally, a list of all found candidates is returned.
	\code{evaluate()} only records hypotheses which at least have a support
	value of \code{support.min} to prevent false positive detections.
	}
\references{
Haase, Daniel (2009): \emph{Detection of Faint Companions around Young Stars in Speckle Patterns of VLT/NACO Cube Mode Images by Means of Post-Processing}. Student research project, Friedrich-Schiller-University of Jena.

\url{http://ringscale.haase-zm.de}}
\author{Daniel Haase}

\examples{
# if analyze() is called without parameters,
# a dialog will let you select the FITS files;
# the result is saved in the working directory
# with a filename composed of the current time,
# as for instance "2009-07-30-134037.RData":
\dontrun{analyze()}

# the filename of the result can be changed by:
\dontrun{analyze(result.filename = "/home/daniel/results/myResult.dat")}

# the filenames of the FITS files may also be
# given as character vector; in this case, no
# dialog will appear:
\dontrun{fn <- list.files("/var/fits", pattern = "\\.[fF][iI][tT][sS]$")}
\dontrun{analyze(filenames = fn)}

# this also applies for the evaluation; if called without
# parameters, a dialog will appear:
\dontrun{evaluate()}

# otherwise, use
\dontrun{evaluate(filename = "/home/daniel/results/myResult.dat")}

# the candidate list returned by evaluate() can be viewed by
\dontrun{candidates <- evaluate()}
\dontrun{print(candidates)}
}

\keyword{methods}