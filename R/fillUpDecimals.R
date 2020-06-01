#' @name fillUpDecimals
#' @aliases fillUpDecimals
#' @title Formats decimal part of numbers.
#' @description The function adds zeros after the decimal point of each element of a numeric vector so that all elements have the same number of digits after the decimal point.
#' 
#' @usage fillUpDecimals(numVec, howManyToFill, fill)
#' 
#' @param numVec A vector of numbers in a numeric or character format.
#' @param howManyToFill The target number of digits after the decimal point. If \code{howManyToFill} is \code{NULL} (the default), then the function uses the maximum number of digits after the decimal point among the element of \code{numVec}.
#' @param fill What character to insert as a filler to reach the target number of digits (characters) after the decimal point. The default is \code{'0'}.
#' 
#' @return A character string with equal number of decimal places after the decimal point.
#' @details The function adds zeros (or other characters) after the decimal point of each element of \code{numVec} so that all elements have the same number of decimal places after the decimal point. Lower and upper case letters are treated as digits (see the examples). \code{NA} values are returned unchanged. It is not recommended to include elements with characters other than digits and letters into \code{numVec} because the function will not work as described for these elements.
#' 
#' @examples
#' fillUpDecimals(c("2", "3.4", "A", NA))
#' fillUpDecimals(c("2", "3.4", "A.5", NA), howManyToFill = 3)
#' fillUpDecimals(c("2", "3.4", "A", NA), fill = "X")
#' fillUpDecimals(c("2", "3.4", "A", NA), howManyToFill = -3)
#' fillUpDecimals(c("2", "3.4", "A.zx", NA), fill = "?")
#' ### It is not recommended to include elements
#' ### with characters other than digits and letters
#' ### because the function will not work as described
#' ### for these elements
#' fillUpDecimals(c("2", "3.4", "A.zx", NA, "#", "#a"), fill = "?")
#' 
#' @importFrom grDevices dev.size gray
#' @importFrom graphics mtext par plot polygon text
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#' @keywords decimal digits
#' @export
###################################################################
######### Adds zeros in decimal places if needed
###################################################################
fillUpDecimals <- function(numVec, howManyToFill=NULL, fill = "0"){
  ### howManyToFill - if NULL then it fill up to the number of maximum decimal points
  ###   for example, if numStr = c(0.33, 34, 4.1) then it fills up to two decimal points
  ###   and the result will be a character string wiht c("0.33", "34.00", "4.10")
  # numStr = as.character(as.numeric(numVec))
  numStr = as.character(numVec)
  noWholeNum = gsub("^[-+ ]*[0-9a-zA-Z]*$", "", numStr, perl=TRUE)
  decimChar = gsub("^[-+ ]*[0-9a-zA-Z]*\\.", "", noWholeNum)
  if(is.null(howManyToFill)){
    howManyToFill = max(nchar(decimChar[!is.na(decimChar)]))
  }
  if(howManyToFill != 0){
    numToFill = pmax(0, howManyToFill - nchar(decimChar))
    fills = sapply(numToFill, function(x) {if(is.na(x)) "" else paste(rep(fill, x), collapse="")} )
    dPoints = sapply(nchar(noWholeNum), function(x){if(!is.na(x) & x==0) "." else ""})
    res = paste(numStr, dPoints, fills, sep="")
    res[res == "NA"] = NA
  }else{
    res = numStr
  }
  res
}
