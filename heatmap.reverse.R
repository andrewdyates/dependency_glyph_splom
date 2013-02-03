## http://hosho.ees.hokudai.ac.jp/~kubo/Rdoc/library/gplots/html/heatmap.2.html
x

# dendrogram control
Rowv = TRUE
Colv = if (symm) "Rowv" else TRUE
distfun = dist
hclustfun = hclust
dendrogram = c("both", "row", "column", "none")
symm = FALSE

# data scaling
scale = c("none", "row", "column"),
na.rm = TRUE,

# image plot
revC = identical(Colv, "Rowv"),
add.expr,

# mapping data to colors
breaks,
symbreaks = min(x < 0, na.rm = TRUE) || scale != "none",

# colors
col = "heat.colors",

# block separation
colsep,
rowsep, 
sepcolor = "white",
sepwidth = c(0.05, 0.05),

# cell labeling
cellnote,
notecex = 1, 
notecol = "cyan",
na.color = par("bg"),

# level trace
trace = c("column", "row", "both", "none"),
tracecol = "cyan",
hline = median(breaks),
vline = median(breaks),
linecol = tracecol,

# Row / Column labeling
margins = c(5, 5),
ColSideColors,
RowSideColors,
cexRow = 0.2 + 1/log10(nr), 
cexCol = 0.2 + 1/log10(nc),
labRow = NULL,
labCol = NULL,

# color key + density info
key = TRUE,
keysize = 1.5,
density.info = c("histogram", "density", "none"),
denscol = tracecol,
symkey = min(x < 0, na.rm = TRUE) || symbreaks,
densadj = 0.25,

# plot labels
main = NULL, 
xlab = NULL,
ylab = NULL,

# plot layout
lmat = NULL,
lhei = NULL,
lwid = NULL,
...

# ------------------------------
    # linear scaling function
    scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
    }
    # return value
    retval <- list()
    # get parameter default values
    scale <- if (symm && missing(scale)) 
        "none"
    else match.arg(scale)
    dendrogram <- match.arg(dendrogram)
    trace <- match.arg(trace)
    density.info <- match.arg(density.info)
    if (length(col) == 1 && is.character(col)) 
        col <- get(col, mode = "function")
    if (!missing(breaks) && (scale != "none")) 
        warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
    if (is.null(Rowv) || is.na(Rowv)) 
        Rowv <- FALSE
    if (is.null(Colv) || is.na(Colv)) 
        Colv <- FALSE
    else if (Colv == "Rowv" && !isTRUE(Rowv)) 
        Colv <- FALSE
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
        stop("`x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
        stop("`x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
        stop("`margins' must be a numeric vector of length 2")
    ### If no cellnotes, make a matrix of empty strings
    if (missing(cellnote)) 
        cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
    if (!inherits(Rowv, "dendrogram")) {
        if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
            c("both", "row"))) {
            if (is.logical(Colv) && (Colv)) 
                dendrogram <- "column"
            else dedrogram <- "none"
            warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting row dendogram.")
        }
    }
    if (!inherits(Colv, "dendrogram")) {
        if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
            c("both", "column"))) {
            if (is.logical(Rowv) && (Rowv)) 
                dendrogram <- "row"
            else dendrogram <- "none"
            warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
                dendrogram, "'. Omitting column dendogram.")
        }
    }
    if (inherits(Rowv, "dendrogram")) {
        ddr <- Rowv
        rowInd <- order.dendrogram(ddr)
    }
    # end
