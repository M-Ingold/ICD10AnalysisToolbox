
#' ICD-10-CM to Clinical Classifications Software Refined (CCSR) mapping
#'
#' @description
#' Mapping of ICD-10-CM codes to Clinical Classifications Software Refined
#' (CCSR) codes, as well as descriptions and categories,
#' adapted from CCSR v2025.1
#'
#' @format ## `DXCCSR`
#' A data frame with 90,768 rows and 6 columns:
#' \describe{
#'   \item{X.ICD.10.CM.CODE.}{ICD-10 code}
#'   \item{ICD.10.CM.CODE.DESCRIPTION}{ICD-10 description}
#'   \item{X.Default.CCSR.CATEGORY.IP.}{Default CCSR ID}
#'   \item{Default.CCSR.CATEGORY.DESCRIPTION.IP}{Default CCSR description}
#'   \item{X.CCSR.CATEGORY.1-5}{Assigned category or categories, descending in order of relevance}
#'   \item{Rationale.for.Default.Assignment}{Rationale for default CCSR code assignment}

#' }
#'
#' @source US Agency for Healthcare Research and Quality <https://hcup-us.ahrq.gov/toolssoftware/ccsr/dxccsr.jsp>
#' (Fiscal Year 2025 version)
#'
#' @references HCUP Clinical Classifications Software Refined (CCSR) for ICD-10-CM diagnoses, v2025.1. Healthcare Cost and Utilization Project (HCUP). Agency for Healthcare Research and Quality, Rockville, MD. hcup-us.ahrq.gov/toolssoftware/ccsr/dxccsr.jsp. Accessed January 2025.
"DXCCSR"
