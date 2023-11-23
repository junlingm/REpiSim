Prior = R6::R6Class(
  "Prior",
  
  private = list(
    dist = NULL,
    .call = NULL,
    args = NULL,
    .parms = NULL
  ),
  
  public = list(
    initialize = function(dist) {
      if (!is.call(dist))
        stop("dist must be a function call")
      l = as.list(dist)
      l = c(l[[1]], NA, l[-1], log=TRUE)
      private$.call = as.call(l)
      if (!is.call(dist)) 
        stop("dist must be a function call")
      e = Expression$new(dist)
      private$.parms = e$parms
    },
    
    log.density = function(x, par) {
      call = private$.call
      call[2] = x
      eval(call, envir = if (length(private$.parms)==0) NULL else as.list(par))
    }
  ),
  
  active = list(
    parms = function() {
      private$.parms
    }
  )
)