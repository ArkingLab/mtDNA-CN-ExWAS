# scientific_10x ----------------------------------------
#' Format numbers to scientific notation (10^#)
#' 
#' Formats numbers to scientific notation using the form #.# x 10^# for display
#' on axes. Leading zeros and '+' signs are removed from the output and spaces
#' are inserted to match all exponents to uniform length for uniform display.
#' @param values A numeric vector
#' @param digits A single integer value specifying the number of digits to
#'   display after the decimal, trailing zeroes will be preserved
#' @return An expression
#' @export
#'
#' @examples
#' x <- 1:3
#' y <- c(0.1, 100, 1000)
#' 
#' # Base Plotting
#' originalMargins <- par()$mar
#' plot(y ~ x, 
#'      axes = F, 
#'      par(mar = c(5, 5, 4, 2) + 0.1), ylab = NA)
#' axis(1)
#' axis(2, 
#'      at = y, 
#'      labels = scientific_10x(y), 
#'      las = 2)
#' par(mar = originalMargins)
#' 
#' # ggplot2
#' ggplot2::qplot(x, y) + 
#'   ggplot2::scale_y_continuous(labels = scientific_10x)
#' ggplot2::qplot(x, y) + 
#'   ggplot2::scale_y_continuous(breaks = y, 
#'                               labels = scientific_10x(y, digits = 2))
#' 

scientific_10x <- function(values, digits = 1) {
  if(!is.numeric(values)){
    stop("values must be numbers")
  }
  if(grepl("^\\d{2}$", digits)){
    stop("digits must a one or two digit whole number")
  }
  
  x <- sprintf(paste0("%.", digits, "e"), values)
  
  x <- gsub("^(.*)e", "'\\1'e", x)
  
  longestExponent <- max(sapply(gregexpr("\\d{1,}$", x), attr, 'match.length'))
  zeroTrimmed <- ifelse(longestExponent > 2,
                        paste0("\\1", paste(rep("~", times = longestExponent-1), collapse = "")),
                        "\\1")
  x <- gsub("(e[+|-])[0]", zeroTrimmed, x)
  
  x <- gsub("e", "~x~10^", x)
  
  if(any(grepl("\\^\\-", x))){
    x <- gsub("\\^\\+", "\\^~~", x)
  } else {
    x <- gsub("\\^\\+", "\\^", x)
  }
  # return this as an expression
  # parse(text=x)
  return(x)
} 
