% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tte_icd.R
\name{create_tte_icd}
\alias{create_tte_icd}
\title{Convert temporal ICD-10 codes into usable information for Cox regression}
\usage{
create_tte_icd(
  dataframe_ID_BL,
  icd_table_long,
  tbl_mortality_contact,
  icd_codes = c("I21"),
  variable_name = NULL,
  censor_time = 5
)
}
\arguments{
\item{dataframe_ID_BL}{A dataframe with a participant ID column. Dataframe gets returned with ICD-10 information is added to it}

\item{icd_table_long}{A long table with ICD-10 information, including the ID, ICD-10 code and date on which it was diagnosed}

\item{tbl_mortality_contact}{A dataframe containing participant ID, baseline visit date and last contact/death date columns.}

\item{icd_codes}{One or multiple ICD-10 codes of interest}

\item{variable_name}{Desired column names. Defaults to concatenated icd_codes}

\item{censor_time}{Time at which events are censored}
}
\value{
Returns dataframe_ID_BL, with the outcome (event yes/no + time to event/censoring) and prevalence at baseline added.
}
\description{
Convert temporal ICD-10 codes into usable information for Cox regression
}
