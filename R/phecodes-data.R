#' ICD-10-CM to PheCode mapping
#'
#' @description
#' Mapping of ICD-10-CM codes to Phecodes, as well as descriptions and categories,
#' adapted from Phecode Map 1.2
#'
#' @format ## `phecodes`
#' A data frame with 90,768 rows and 6 columns:
#' \describe{
#'   \item{ICD}{ICD-10 code}
#'   \item{Flag}{ICD version (always 10)}
#'   \item{ICDString}{ICD-10 description}
#'   \item{Phecode}{Phecode as number}
#'   \item{PhecodeString}{Phecode description}
#'   \item{PhecodeCategory}{Disease category}
#' }
#'
#' @source PheWAS Resources <https://phewascatalog.org/phewas/_w_7b1256e92bb84559bf1c493fd51cf540/data/Phecode_map_v1_2_icd9_icd10cm.csv.zip>
#' (downloaded June 2026)
#'
#' @references Denny, J., Bastarache, L., Ritchie, M. et al. Systematic comparison of phenome-wide association study of electronic medical record data and genome-wide association study data. Nat Biotechnol 31, 1102–1111 (2013).
#' (\doi{10.1038/nbt.2749})
"phecodes"



