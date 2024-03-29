#' Daily changes in the S&P 500 index
#' 
#' @description 
#' This dataset contains the quantised daily changes \ifelse{html}{\out{x<sub>i</sub>}}{x_i} in the Standard & Poor’s index price, 
#' from January 2, 1928 until October 7, 2016. The price changes are quantised to 7 values as follows: 
#' If the change between two successive trading days (i \ifelse{html}{\out{-}}{-} 1) and i is smaller than 
#' \ifelse{html}{\out{-3&percnt;}}{-3\%}, \ifelse{html}{\out{x<sub>i</sub>}}{\eqn{x_i}} is set equal to 0; if the change is between \ifelse{html}{\out{-3&percnt;}}{-3\%} and \ifelse{html}{\out{-2&percnt;}}{-2\%}, 
#' \ifelse{html}{\out{x<sub>i</sub>&equals;1}}{\eqn{x_i = 1}}; for 
#' changes in the intervals (\ifelse{html}{\out{-2&percnt;}}{-2\%}, \ifelse{html}{\out{-1&percnt;}}{-1\%}], (\ifelse{html}{\out{-1&percnt;}}{-1\%}, 1\ifelse{html}{\out{&percnt;}}{\%}], 
#' (1\ifelse{html}{\out{&percnt;}}{\%}, 2\ifelse{html}{\out{&percnt;}}{\%}], and (2\ifelse{html}{\out{&percnt;}}{\%}, 3\ifelse{html}{\out{&percnt;}}{\%}]
#' \ifelse{html}{\out{x<sub>i</sub>}}{\eqn{x_i}} is set 
#' equal to 2,3,4 and 5, respectively; and for changes greater than 3\ifelse{html}{\out{&percnt;}}{\%}, \ifelse{html}{\out{x<sub>i</sub>}}{\eqn{x_i}} = 6.  
#'  
#' @docType data
#' @format An object of class \code{"character"}. 
#'
#' @keywords datasets
#'
#' @references {Yahoo! finance.}
#' (\href{https://finance.yahoo.com}{yahoo_finance})
#' @examples 
#' BCT(SP500, 10)
#' 
"SP500"