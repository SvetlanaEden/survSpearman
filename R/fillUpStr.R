#' @name fillUpStr
#' @aliases fillUpStr
#' @title Adds characters to the head or tail of each element of the character vector.
#' @description Adds characters to the head or tail of each element of the character vector to make all elements the same length.
#' 
#' @usage fillUpStr(strVec, howManyToFill, where, fill)
#' 
#' @param strVec A vector of character strings.
#' @param howManyToFill The target number of characters in each element of \code{strVec}. If \code{howManyToFill} is \code{NULL} or less than the number of characters in the longest element of \code{strVec}, then the number of characters in the longest element is used.
#' @param where Where to place the characters: at the beginning of each element (\code{'head'}, the default) or at the end of each element (\code{'tail'})
#' @param fill What character to add. The default is a blank space.
#' 
#' @return A character string with equal number of characters in each element.
#' @details The function adds characters to the head or tail of each element of the character vector, \code{strVec}, so that all elements have the same number of characters. \code{NA} values are returned unchanged. If \code{howManyToFill} is \code{NULL} or less than the number of characters in the longest element of \code{strVec}, then the number of characters in the longest element is used.
#' 
#' @examples
#' fillUpStr(c("A", "aaa", NA, "bb", "4.5"), where="tail", howManyToFill = 4, fill = "_")
#' fillUpStr(c("A", "aaa", NA, "bb", "4.5"), where="tail", howManyToFill = -3, fill = "_")
#' fillUpStr(c("A", "aaa", NA, "bb", "4.5"), where="head", howManyToFill = -3, fill = "*")
#' 
#' @author Svetlana K Eden, \email{svetlanaeden@gmail.com}
#' @keywords character vector
#' @export
###################################################################
######### Adds characters to the head or tail of each element of the character vector
###################################################################
fillUpStr <- function(strVec, howManyToFill=NULL, where="head", fill = " "){
  ### "strVec" to fill
  ### howManyToFill - if NULL then it fill up to the number of maximum characters
  ### "where" can be either "head" or "tail"
  if(!any(where %in% c("head", "tail"))) stop("Argument 'where' can be either 'head' or 'tail'\n")
  if(is.null(howManyToFill)){
    howManyToFill = diff(range(nchar(strVec), na.rm = TRUE))
  }
  # numToFill = max(cnchar(strVec), na.rm = TRUE) - nchar(strVec)
  numToFill = max(max(nchar(strVec), na.rm = TRUE), howManyToFill) - nchar(strVec)
  numToFill[is.na(numToFill)] = 0
  fills = sapply(numToFill, function(x)paste(rep(fill, x), collapse=""))
  if (where == "head"){
    res = paste(fills, strVec, sep="")
  }else{
    res = paste(strVec, fills, sep="")
  }
  res[is.na(strVec)] = NA
  res
}
