% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{simCYD14_twostage}
\alias{simCYD14_twostage}
\title{Simulated CYD14 Two-Phase Data Set}
\format{A data frame simulated to mimic a phase 3 dengue vaccine efficacy trial with 10250 rows and the following 40 variables:
\describe{
  \item{X}{Row number}
  \item{STUDYID}{Study name}
  \item{SUBJID}{Subject ID}
  \item{VACC}{Vaccination status: vaccine=1; placebo=0}
  \item{AGEYRS}{Subject age at enrollment}
  \item{SEX}{M=male; F=female}
  \item{Male}{Binary indicator: 1=male, 0=female}
  \item{COUNTRY}{Country code}
  \item{IDN-MEX}{Country code indicator variables}
  \item{VCD}{Any VCD (virally confirmed dengue)=1, no VCD=0}
  \item{VCD_Sero1-VCD_Sero4}{Indicator of VCD specific serotype 1 - serotype 4}
  \item{baseline_PRNT_Sero1-baseline_PRNT_Sero4}{Baseline serotype PRNT titers for serotypes 1-4}
  \item{baseline_PRNT_SeroAverage}{Geometric mean of baseline PRNT titers}
  \item{M13_PRNT_Sero1-M13_PRNT_Sero4}{M13 serotype PRNT titers for serotypes 1-4}
  \item{M13_PRNT_SeroAverage}{Geometric mean of M13 PRNT titers}
  \item{M13_MNv2_Sero1-M13_MNv2_Sero4}{M13 serotype MNv2 titers for serotypes 1-4}
  \item{M13_MNv2_SeroAverage}{Geometric mean of M13 MNv2 titers}
  \item{weight}{Observation level sampling weight for sampling to phase-two}
  \item{delta}{Indicator of sampling to phase-two sample}
}}
\usage{
data(simCYD14_twostage)
}
\description{
Data simulated to mimic the CYD14 Tetravalent Dengue Vaccine Trial conducted by Sanofi-Pasteur
}
\references{
Capeding et al. (2014) The Lancet 384:1358--1365
}
\keyword{datasets}
