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
      dist["log"]=TRUE
      private$.call = dist
      if (!is.call(dist)) 
        stop("dist must be a function call")
      e = Expression$new(dist)
      private$.parms = e$parms
    },
    
    log.density = function(x, par) {
      eval(private$.call, envir = if (length(private$.parms)==0) NULL else as.list(par))
    }
  ),
  
  active = list(
    parms = function() {
      private$.parms
    }
  )
)