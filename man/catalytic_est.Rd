% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/single_kp.R
\name{catalytic_est}
\alias{catalytic_est}
\title{Estimation of single catalytic constant}
\usage{
catalytic_est(method = T, time, species, enz, subs, MM, catal,
  nrepeat = 10000, jump = 1, burning = 0, catal_m = 1,
  catal_v = 1e+06)
}
\arguments{
\item{method}{method selection: T=TQ model, F=SQ model(default = T)}

\item{time}{observed time interval}

\item{species}{observed trajectory of product}

\item{enz}{enzyme concentration}

\item{subs}{substrate concentration}

\item{MM}{true value of MM constant}

\item{catal}{initial value of catalytic constant}

\item{nrepeat}{total number of iteration (default=10000)}

\item{jump}{length of distance (default =1)}

\item{burning}{lenth of burning period (default =0)}

\item{catal_m}{prior mean of gamma prior (default =1)}

\item{catal_v}{prior variance of gamma prior (default =1e+6)}
}
\value{
A vector of posterior samples of catalytic constant
}
\description{
The function estimates catalytic constant using input data set
and enzyme, substrate concentration, and MM constant.
}
\details{
The function catalytic_est generates a set of Monte Carlo
simulation samples from posterior distribution of catalytic
constant of enzyme kinetics model. Because the function considers
catalytic constant as a parameter to be estimated, the user should
input three constants of enzyme concentration, substrate concentration
and Michaelis-Menten constant.
}
\examples{
data("Chymo_low")
time1=Chymo_low[,1]
species1=Chymo_low[,2]
Chymotrypsin.low<-catalytic_est(method=TRUE,time=time1,species=species1,enz=4.4e+7,subs=4.4e+7
                 ,MM=4.4e+8, catal=0.1,nrepeat = 10000)
}
