#' declare a regular polygon for the shape of a compartment in a flow char.
#' @param n the number of sides for the polygon.
#' @return a character holding the tikz option for a node
#' @export
Regular.polygon <- function(n) {
  paste0("regular polygon,regular polygon sides=", n)
}

#' declare a regular star for the shape of a compartment in a flow char.
#' @param n the number of sides for the star
#' @return a character holding the tikz option for a node
#' @export
Star <- function(n, ratio=NULL) {
  s = paste0("star,star points=", n)
  if (is.null(ratio)) s else paste0(s, ",star point ratio=", ratio)
}

Node.shapes = c("rectangle", "circle", "ellipse", "split circle",
                "diamond", "crossout", "strikeout", "forbidden sign")

#' declare a node for a flowchart
#'
#' A node in a flowchart represents a compartment.
#' 
#' @param name the name of the compartment
#' @param label the latex label for the compartment
#' @param at a numeric vector giving the coordinate of the node
#' @param shape a character giving the shape of the compartment
#' @param rounded.corners if the shape should use rounded corners
#' @param fill an R color for filling the shape
#' @param ... other properties for the node
#' @return a character holding the tikz command for a node
#' @export
Node <- function(name, label, at,
                 shape=Node.shapes, 
                 rounded.corners = TRUE,
                 fill=NULL, ...) {
  if (!is.numeric(at) && length(at) != 2)
    stop("invalid coordinates (", at, ")")
  s = if (is.null(shape)) c() else {
    m = pmatch(shape[[1]], Node.shapes)
    if (is.na(m)) shape else paste0("draw,", Node.shapes[[m]])
  }
  if (rounded.corners) s = c(s, "rounded corners")
  if (!is.null(fill)) s = c(s, paste0("fill=", fill))
  paste0("\\node[", paste(c(s, ...), collapse=","), "] at (",
         at[1], ",", at[2], ") (", name, ") {$", label, "$};")
}

Label.directions = c("above", "below", "left", "right")

#' specify the properties of a label in a flowchart
#' 
#' A label is the rate of a transition.
#' 
#' @param pos a value in the interval [0,1] specifying the position
#' of the label on the arc, 0 is at the starting point, 1 is at the
#' end point. Default to 0.5.
#' @param direction the relative position of the label to the arc
#' @return a character holding the tikz command for the property of the label
#' @export
Label <- function(pos=NULL, direction=Label.directions) {
  args = c()
  if (!is.null(pos)) {
    if (!is.numeric(pos) || pos < 0 || pos > 1)
      stop("pos must be a number between 0 and 1")
    args = c(args, paste0("pos=", pos))
  }
  if (!is.null(direction))
    args = c(args, match.arg(direction, Label.directions))
  paste(args, collapse=",")
}

#' specify the properties of an arc in a flowchart
#' 
#' An arc represents a transition.
#' 
#' @param from the name of the compartment that this flow originates.
#' @param to the name of the compartment that this flow goes to.
#' @param rate the rate of the flow, i.e., the label. An R expression
#' @param label the property of the label, speficied by the Label function.
#' @param edge ghe property for the arc, such as "bend left" (or "(") or "bend right" (or ")").
#' @return a character holding the tikz command for the property of the label
#' @export
Arc <- function(from, to, rate, label=NULL, edge=NULL) {
  if (!is.null(edge)) {
    if (edge == "(") edge = "bend left" else 
      if (edge == ")") edge = "bend right"
  }
  if (is.numeric(from) && length(from) >= 2) 
    from = paste(from[1:2], collapse=",")
  if (is.numeric(to) && length(to) >= 2)
    to = paste(to[1:2], collapse=",")
  paste0("\\path[->] (", from, ") edge[", edge, 
         "] node[", label, "] {\\small$", rate, "$} (", to, ");")
}

#' generate the flow chart for a compartmental model
#' @param model an object of Compartmental class
#' @param tex an object of TexFormatter class
#' @param ... the tikz commands. Compartment descriptions and 
#' transition descriptions must be given by named arguments, which names
#' must be the name for the compartment or transition, and the value is 
#' a description given by Node and Arc functions.
#' @param preambles a character holding a sequence of latex commands as latex preambles.
#' @return a character holding a tikzpicture latex environment.
#' @details When calling the Node function, the name and label of the compartment
#' does not need to be passed. The name it is given as the name of the argument.
#' The label is generated using tex. Similarly, when calling Arc, 
#' the from, to and rate need not be specified.
#' @examples
#' # an SIR model
#' SIR = Compartmental$new(S, I, R, title="SIR")
#' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
#' SIR$transition(I->R ~ gamma*I, name="recovery")
#' flowchart(SIR, TexFormatter$new(),
#'   S=Node(at=c(0,0)),
#'   I=Node(at=c(2,0)),
#'   R=Node(at=c(4,0))
#' )
#' @export
flowchart = function(model, tex, ..., preambles=NULL) {
  args = as.list(substitute(list(...)))[-1]
  compartments = model$compartments
  transitions = model$transitions
  a = NULL
  if (length(args) > 0) for (i in 1:length(args)) {
    n = names(args[i])
    v = args[[i]]
    if (n %in% compartments) {
      l = if (!is.call(v) || v[[1]] != "Node") list(at=c(0,0)) else 
        as.list(v)[-1]
      l[["name"]] = n
      l[["label"]] = tex$typeset(as.name(n))
      args[n] = do.call("Node", l)
    } else if (!is.null(transitions[[n]])) {
      l = if (!is.call(v) || v[[1]] != "Arc") list() else as.list(v)[-1]
      e = transitions[[n]]
      if (is.null(l[["from"]])) l[["from"]] = e$from
      if (is.null(l[["to"]])) l[["to"]] = e$to
      if (is.null(l[["rate"]])) l[["rate"]] = tex$typeset(e$rate)
      args[n] = do.call("Arc", l)
    } else if (n == "") {
      a = c(a, v)
    } else stop("Undefined compartment or transition ", n)
  }
  
  missing = Filter(
    function(x) !is.null(x),
    lapply(
      compartments,
      function(C) if (is.null(args[[C]])) C else NULL
    )
  )
  if (length(missing) > 0)
      stop(
        "The position of the compartment",
        if (length(missing) > 1) " " else "s ",
        paste(missing, collapse=", "),
        if (length(missing) > 1) " is " else " are ",
        "not specified"
      )

  chart = c(
    list("tikzpicture"), 
    args[compartments],
    lapply(
      transitions,
      function(e) {
        if (!is.null(args[[e$name]])) args[[e$name]] else 
          Arc(from=e$from, to=e$to, rate=tex$typeset(e$rate))
      }
    )
  )
  paste(c(preambles, do.call(tex$environment, chart)), collapse ="\n")
}

#' change the extension of a file
#' @param file a character giving the file path
#' @param ext.from the extension to change from
#' @param ext.to the extension to change to
#' @details if the file has the extension given by ext.from, then 
#' the extension is changed. Otherwise, the extension given by ext.to
#' is appended to file.
#' @examples
#' # both of the following gives "test.pdf"
#' change.extension("test.tex", ".tex", ".pdf")
#' change.extension("test", ".tex", ".pdf")
#' @export
change.extension <- function(file, ext.from, ext.to) {
  if (is.null(file)) return(NULL)
  f = tolower(file)
  if (substring(f, nchar(f)-nchar(ext.from) + 1) == tolower(ext.from))
  paste0(substring(file, 1, nchar(f)-nchar(ext.from)), ext.to)
}

#' generate a PDF file from the given latex commands
#' 
#' @param file the path to the to-be-generated pdf file.
#' @param ... the latex commands.
#' @param fontsize the fontsize in latex points.
#' @param preambles the latex preambles
#' @details the document class is "standalone", the tikz package is
#' automatically loaded.
#' @examples
#' # an SIR model
#' SIR = Compartmental$new(S, I, R, title="SIR")
#' SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
#' SIR$transition(I->R ~ gamma*I, name="recovery")
#' tikz.pdf("SIR.pdf", flowchart(SIR, TexFormatter$new(),
#'   S=Node(at=c(0,0)),
#'   I=Node(at=c(2,0)),
#'   R=Node(at=c(4,0))
#' ))
#' @export
tikz.pdf <- function(file, ..., fontsize=12, preamble=NULL) {
  if (!require(tools, quietly = TRUE)) {
    stop("The tools package must be installed before using tikz.pdf.")
  }
  tex = TexFormatter$new()
  file = change.extension(file, ".pdf", ".tex")
  
  doc = c(
    tex$command("documentclass", "standalone", option="tikz"),
    preamble,
    tex$environment(
      "document",
      tex$command(
        "fontsize", 
        paste0(fontsize, "pt"), 
        round(fontsize*1.2,1)
      ),
      tex$command("selectfont"),
      ...
    )
  )
  doc = paste(doc, collapse="\n")
  if (!is.null(file)) sink(file)
  cat(doc, "\n")
  if (!is.null(file)) sink()
  if (!is.null(file)) {
    dir = getwd()
    setwd(dirname(tools::file_path_as_absolute(file)))
    tryCatch(
      tools::texi2pdf(file, clean=TRUE),
      error = function(e) NULL)
    setwd(dir)
    file.remove(file)
  }
}

if (exists("TEST") && is.logical(TEST) && TEST) {
  SIR = Compartmental$new(S, I, R, title="SIR")
  SIR$transition(S->I ~ beta*S*I/N, N=S+I+R, name="infection")
  SIR$transition(I->R ~ gamma*I, name="recovery")
  tikz.pdf("test.pdf", flowchart(SIR, TexFormatter$new(),
   S=Node(at=c(0,0)),
   I=Node(at=c(2,0)),
   R=Node(at=c(4,0))
  ))
}
