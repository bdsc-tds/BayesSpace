#' @export
#' @importFrom magrittr `%>%`
#' @importFrom dplyr arrange distinct
pairwiseComp <- function(data, func.comp, func.data = NULL, to.mat = TRUE) {
  if (!is.null(func.data)) {
    data <- sapply(
      data,
      function(x) func.data(x),
      simplify = FALSE
    )
  }

  ret <- do.call(
    rbind,
    lapply(
      names(data),
      function(x) {
        do.call(
          rbind,
          lapply(
            names(data),
            function(y) {
              data.frame(
                var1 = x,
                var2 = y,
                val = func.comp(data[[x]], data[[y]])
              )
            }
          )
        )
      }
    )
  ) %>%
    arrange(
      var1, var2
    )

  if (to.mat) {
    ret <- matrix(
      ret$val,
      byrow = TRUE,
      nrow = sqrt(dim(ret)[1]),
      dimnames = list(
        row = distinct(ret, var1)$var1,
        col = distinct(ret, var2)$var2
      )
    )
  }

  ret
}

#' @export
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @importFrom grid grid.rect gpar grid.text
plotHeatmap <- function(data, name,
                        col.func = colorRamp2(c(0, 1), c("white", "red")),
                        plot.triangle = TRUE, plot.diag = FALSE,
                        anno.cells = TRUE, cluster_rows = FALSE,
                        cluster_columns = FALSE, ...) {
  cell.func <- NULL
  if (any(c(plot.triangle, anno.cells))) {
    cell.func <- function(j, i, x, y, w, h, fill) {
      if (i >= j || (!plot.triangle && i < j)) {
        grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))

        if (anno.cells) {
          grid.text(sprintf("%.2f", data[i, j]), x, y, gp = gpar(fontsize = 10))
        }
      }
    }
  }

  if (plot.triangle && !plot.diag) {
    data <- data[-1, -dim(data)[2]]
  }

  Heatmap(
    data,
    name = name,
    rect_gp = gpar(type = "none"),
    col = col.func,
    cell_fun = cell.func,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    row_names_side = "left",
    column_names_rot = 45,
    ...
  )
}
