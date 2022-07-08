#' @import data.table
#' @importFrom magrittr %>% set_rownames extract
#' @importFrom foreach %do% foreach %dopar%
#' @importFrom doRNG %dorng%
#' @importFrom fpCompare %==% %!=%
# #' @importFrom utils combn
#' @importFrom stats cor sd
NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
