# produces a sequence from the minimum to the maximum of x and same length as x

nseq <- function(x, len = length(x)) {seq(min(x, na.rm = TRUE), 
                                          max(x, na.rm = TRUE), length = len)}

