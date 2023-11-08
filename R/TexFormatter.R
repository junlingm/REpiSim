#' A S3 method to merge two lists
#' @param a a list
#' @param b a list
#' @return a merged list 
#' @export
merge.list = function(a, b) {
  if (length(b) == 0) return(a)
  if (length(a) == 0) return(b)
  nb = names(b)
  for (i in 1:length(b)) {
    n = nb[[i]]
    if (!is.null(n) && n != "") {
      a[n] = b[n] 
    } else if (!(b[[i]] %in% a)) {
      a = c(a, b[i])
    }
  }
  a
}

#' R6 class to typeset R expressions as latex equations
#' 
#' @docType class
#' @export
TexFormatter = R6Class(
  "TexFormatter",
  private = list(
    # the brackets to be used to enclose a matrix
    matrix.bracket = c("[", "]"),
    # the default symbols
    default = list(
      sin = "\\sin",
      cos = "\\cos",
      tan = "\\tan",
      cot = "\\cot", 
      sec = "\\sec", 
      csc = "\\csc",
      alpha = "\\alpha",
      beta = "\\beta",
      gamma = "\\gamma",
      delta = "\\delta",
      epsilon = "\\varepsilon",
      zeta = "\\zeta",
      eta = "\\eta",
      theta = "\\theta",
      iota = "\\iota",
      kappa = "\\kappa",
      lambda = "\\lambda",
      mu = "\\mu",
      nu = "\\nu",
      ksi = "\\ksi",
      pi = "\\pi",
      sigma = "\\sigma",
      tau = "\\tau",
      phi = "\\phi",
      chi = "\\chi",
      psi = "\\psi",
      omega = "\\omega",
      Gamma = "\\Gamma",
      Delta = "\\Delta",
      Theta = "\\Theta",
      Lambda = "\\Lambda",
      Pi = "\\Pi",
      Sigma = "\\Sigma",
      Phi = "\\Phi",
      Psi = "\\Psi",
      Omega = "\\Omega",
      infty = "\\infty"
    ),
    
    # returns the tex code for x 
    substitute = function(x, keep = FALSE) {
      r = self$symbols[[x]]
      if (is.null(r)) r = private$default[[x]]
      if (!is.null(r)) r else if (keep) x else NULL
    },
    
    # if x is a bracketed expression, remove the brackets
    remove.bracket = function(x) {
      if (is.call(x) && x[[1]] == "(")
        private$remove.bracket(x[[2]]) else x
    }
  ),

  public = list(
    #' @field symbols a named list that matches R names (names in the list)
    #' to latex symbols (values in the list). Trig functions, greek letters
    #' and infinity are altomatically defined. They do not need to be defined 
    #' in this list.
    symbols = list(),

    #' @description constructor
    #' @param symbols a named list that matches R names (names in the list)
    #' to latex symbols (values in the list). Trig functions, greek letters
    #' and infinity are automatically defined. They do not need to be defined 
    #' in this list. But they may be redefined here
    #' @examples 
    #' tex = TexFormatter$new(II="[I\gets I]")
    #' tex$typeset(quote(sin(alpha*I)/II))
    initialize = function(symbols=list()) {
      self$symbols = symbols
    },
    
    #' @description set/change the tex symbols of R names
    #' 
    #' @param symbols a named list that matches R names (names in the list)
    #' to latex symbols (values in the list). Trig functions, greek letters
    #' and infinity are automatically defined. They do not need to be defined 
    #' in this list. But they may be redefined here
    #' @return an invisible self to chain methods
    #' @examples
    #' tex = TexFormatter$new(II="[I\gets I]")
    #' tex$set.symbols(list(alpha="a")$typeset(quote(sin(alpha*I)/II))
    set.symbols = function(symbols) {
      self$symbols = merge(self$symbols, symbols)
      invisible(self)
    },

    #' description typeset an R expression
    #' @param R the R expression to be typeset.
    #' @return a character holding the latex commands
    #' @examples
    #' tex = TexFormatter$new(II="[I\gets I]")
    #' tex$typeset(quote(sin(alpha*I)/II))
    typeset = function(R) {
      if (is.call(R)) {
        l = as.list(R)
        switch(
          as.character(l[[1]]),
          `'` = paste0(self$typeset(l[[2]]), "'"),
          `+` = if (length(l) == 2) self$typeset(l[[2]]) else
            paste(self$typeset(l[[2]]), self$typeset(l[[3]]), sep = " + "),
          `-` = if (length(l) == 2) {
            paste("-", self$typeset(l[[2]]))
          } else paste(self$typeset(l[[2]]), "-", self$typeset(l[[3]])),
          `*` = paste(self$typeset(l[[2]]), self$typeset(l[[3]])),
          `/` = paste0("\\frac{", self$typeset(l[[2]]), "}{", 
                       self$typeset(private$remove.bracket(l[[3]])), "}"),
          `^` = paste0(self$typeset(l[[2]]), "^{", 
                       self$typeset(private$remove.bracket(l[[3]])), "}"),
          `(` = paste0("(", self$typeset(l[[2]]), ")"),
          paste0(self$typeset(l[[1]]), "(", 
                do.call(paste, c(lapply(l[-1], self$typeset), list(sep=", "))),
                ")")
          )
      } else if (is.name(R)) {
        s = private$substitute(R)
        if (!is.null(s)) s else {
          tex = ""
          m = regexec("([a-z]+)(?:([_^])?(.+))*", R, ignore.case=TRUE, perl=TRUE)[[1]]
          if (m[[1]] == -1 || (m[[3]] == 0 && m[[4]] == 0)) as.character(R) else {
            len = attr(m, "match.length")
            mid = if (m[[3]] == 0) "_" else
              substring(R, m[[3]], m[[3]] + len[[3]] - 1)
            paste0(
              tex, 
              private$substitute(substring(R, m[[2]], m[[2]] + len[[2]] - 1), TRUE),
              mid,
              "{", private$substitute(substring(R, m[[4]], m[[4]] + len[[4]] - 1), TRUE), "}"
            )
           }
        }
      } else if (is.matrix(R)) {
        paste(
          paste0("\\left", private$matrix.bracket[[1]]),
          paste0("\\begin{array}{", paste(rep("c", ncol(R)), collapse=""), "}"),
          paste(
            apply(R, 1, function(row) {
              paste(sapply(row, self$typeset), collapse=" & ")
            }),
            collapse = "\\\\\n"
          ),
          "\\end{array}",
          paste0("\\right", private$matrix.bracket[[2]]),
          sep = "\n"
        )
      } else as.character(R)
    },
    
    #' @description enclose the tex commands in brackets
    #' @param ... the tex commands to enclose
    #' @param closure a character vector of the brackets to be used. It defaults to `c("{","}")`.
    #' If the vector is length 1, then it is used as both the opening and closing bracket.
    #' @return the enclosed tex command
    #' @examples 
    #' tex = TexFormatter$new()
    #' tex$enclose("\\bf This is a test")
    #' tex$enclose(tex$typeset.equation(quote(alpha==beta)), 
    #'   closure=c("\\begin{align}", "\\end{align}"))
    enclose = function(..., closure=c("{", "}")) {
      if (length(closure) == 0) {
        closure = c("", "")
      } else if (length(closure) == 1)
        closure = c(closure, closure)
      paste0(closure[[1]], ..., closure[[2]])
    },
    
    #' @description generate a latex command
    #' @param name the name of the latex command.
    #' @param ... the arguments of the latex command.
    #' @param argument a character vector as an alternative way to provide arguments for the latex command
    #' @param option the optional argument to latex command (the one that is enclosed in square brackets.)
    #' @return a character holding the latex command
    #' @examples 
    #' tex=TexFormatter$new()
    #' tex$command("textcolor", "blue", "this text is blue")
    #' tex$command("documentclass", argument="article", option="12pt")
    command = function(name, ..., argument=NULL, option = NULL) {
      tex = paste0(
        "\\", name,
        if (!is.null(option)) self$enclose(option, closure=c("[", "]")) else "", 
        if (!is.null(argument)) self$enclose(argument) else ""
      )
      l = c(...)
      arguments = if (length(l) == 0) " " else
        paste(paste0("{", c(...), "}"), collapse="")
      paste0(tex, arguments)
    },
    
    #' @description generate a latex environment
    #' @param name the name of the latex environment
    #' @param ... the contents of the latex environment
    #' @param argument a named list to provide arguments for the latex environment
    #' @param option the optional argument to latex environment (the one that is enclosed in square brackets.)
    #' @return a character holding the latex environment
    #' @examples
    #' tex = TexFormatter$new()
    #' tex$environment("frame", argument="frame title", "line 1", "line 2")
    environment = function(name, ..., argument=NULL, option = NULL) {
      paste(
        paste0(
          self$command("begin", name),
          if (!is.null(option)) self$enclose(option, closure=c("[", "]")) else "", 
          if (!is.null(argument)) self$enclose(argument) else ""
        ), 
        ..., 
        self$command("end", name),
        sep="\n"
      )
    },

    #' @description typeset a system of equations 
    #' 
    #' This method typesets a system of equations that can be used as the
    #' content of a latex align or eqnarray environments. 
    #' @param equations a list of equations or inequalities. The left and right sides 
    #' of the equations may be connected using "<-" or "==".
    #' @param align can be one of "left", "right", "both", which specifies if the
    #' ampersand (&) is placed on the left, right or both sides of the equal (or other 
    #' comparison) sign. The default value is left. 
    #' @return a character holding the latex commands. 
    #' @details Note that the latex environment itself (e.g., \\begin{align}\\end{align}) is not returned.
    #' @examples 
    #' SIR = Model$new(
    #'   title = "SIR",
    #'   S ~ -beta*S*I/N, # the dS/dt equation
    #'   I ~ beta*S*I/N - gamma*I, # the dI/dt equation
    #'   R ~ gamma*I, # the dR/dt equation
    #'   N = S + I + R # the total population N
    #' )
    #' tex = TexFormatter$new()
    #' tex$typeset.equation(SIR$equations$equations)
    #' tex$typeset.equation(SIR$equations$where)
    typeset.equation = function(equations, align=c("left", "right", "both")) {
      equation <- function(eq, align) {
        if (is.call(eq)) {
          eq = as.list(eq)
        }
        if (!is.list(eq))
          stop("eq should be either a lsit of a call")
        op = as.character(eq[[1]])
        if (! op %in% c("<-", "->", "<<-", "->>", "=", "==", ">", "<", ">=", "<="))
          stop("invalid equation")
        sep = switch(
          op,
          `=` = "=",
          `==` = "=",
          `<-` = "=",
          `->` = "=",
          `->>` = "=",
          `<<-` = "=",
          `>` = ">",
          `<` = "<",
          `>=` = "\\leq ",
          `<=` = "\\geq "
        )
        sep = switch(
          align,
          none = sep,
          left = paste("&", sep),
          right = paste(sep, "&"),
          both = paste("&", sep, "&")
        )
        paste0(self$typeset(eq[[2]]), sep, self$typeset(eq[[3]]))
      }
      align = match.arg(align[[1]], c("none", "left", "right", "both"))
      if (length(equations) == 0) {
        NULL
      } else if (length(equations) == 1) {
        if (is.list(equations)) equations = equations[[1]]
        equation(equations, align = "none")
      } else {
        paste(
          sapply(equations, equation, align=align),
          collapse="\\\\\n"
        )
      }
    },

    #' @description a text representation of the class
    format = function() {
      paste(c(
        "a tex typesetter for R expressions",
        lapply(names(self$symbols), function(s) {
          paste0("  ", s, ": ", self$symbols[[s]])
        })
      ), collapse="\n")
    }
  )
)
