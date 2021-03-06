% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Denoise.R
\name{denoiseImage}
\alias{denoiseImage}
\title{Use a given dictionary D to denoise image.}
\usage{
denoiseImage(I_noise, D, sigma, stepsize = 1)
}
\arguments{
\item{I_noise}{The image to be denoised. In form of matrix.}

\item{D}{D is the dictionary used in Y=DA to denoise.}

\item{sigma}{Noise level.}

\item{stepsize}{(optional) The stepsize when splicting the image. Default is 1}
}
\value{
The denoised image in for of a matrix.
}
\description{
The noisy image is split into sqrt(m) by sqrt(m)
patches. Each patch is vectorized into a column of matrix Y.
Using the given D, the sparse coding A_hat in Y=DA is obtained.
Then Y_denoise=DA_hat. The final
denoised image is reconstruncted on the denoised patches.
}
\details{
See \url{https://arxiv.org/abs/1605.07870}
}
\examples{
I = lena_crop #use a smaller image as an example
## add noise
sigma=20;
I_noise=AddNoise(I,sigma)
## use ODCT dictionary
D0=ODCT(64,100)
## denoise
\dontrun{
I_clean=denoiseImage(I_noise,D0,sigma)
}
}
