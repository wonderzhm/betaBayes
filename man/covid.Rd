\name{covid}
\alias{covid}

\title{COVID-19 County Level Data}
\description{
A county level COVID-19 dataset in US. It is of interest to examine the association between several county-level characteristics and the cumulative numbers of confirmed cases and deaths. County-level characteristics are based on the 2018 ACS 5-year estimates. 
}
   
\usage{data(covid)}

\format{
  \tabular{ll}{
    FIPS:\tab FIPS county code\cr
    PopE:\tab total population \cr
    MaleP: \tab percentage of people who are male \cr
    WhiteP: \tab percentage of people who are white \cr
    BlackP:\tab percentage of people who are black or African American\cr
    Age65plusP:\tab  percentage of people who are 65 years and over\cr
    PovertyP:\tab percentage of people whose income in the past 12 months is below poverty \cr
    RUCC_2013:\tab 2013 Rural Urban Continuum Code, with a higher value indicating a more rural county\cr
    State:\tab two-letter state abbreviation code \cr
    deaths: \tab cumulative number of deaths as of October 13, 2020 \cr
    cases: \tab cumulative number of confirmed cases as of October 13, 2020 \cr
  }
}

\examples{
data(covid)
head(covid)
}

