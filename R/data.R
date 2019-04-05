#' lag_dat
#' 
#' Simulated Sequential Data
#'
#' Data were simulated based on code at:
#' https://stephens999.github.io/fiveMinuteStats/simulating_discrete_chains_1.html
#'
#' @docType data
#'
#' @usage data(lag_dat)
#'
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'   \item{seqvar}{A vector containing sequential data--3 values: A, B, and C}
#'   \item{group}{group membership for the seqvar column--groups are 1 and 2}
#'   \item{linenum}{line numbers for the seqvar column}
#'}
#' @keywords datasets
#'
#' @references https://stephens999.github.io/fiveMinuteStats/simulating_discrete_chains_1.html
#'
#' @examples
#' data(lag_dat)
#'
#' LagTables::trprobs(lag_dat, lagvar = "seqvar")
#' LagTables::trprobs(lag_dat, lagvar = "seqvar", laggroup="group")
#' LagTables::lagmodels(lag_dat, lagcol="seqvar", lagnum=2)
"lag_dat"
