## This file created by Hassan Pazira

#' @export

b.diag <- function(..., fill = 0) {
  len <- ...length()
  if ((len) == 0L)
    stop("Enter individual matrices or a list of matrices.")
  else if (len > 1L)
    matrices <- list(...)
  else if (!is.list(firstb <- ..1))
    matrices <- list(firstb)
  else if (length(firstb) == 1L)
    matrices <- firstb
  else matrices <- firstb
  if (is.list(matrices) &
      !all(unlist(lapply(matrices, is.matrix))))
    matrices <- lapply(matrices, as.matrix)
  if (length(sapply(matrices, function (x) is.matrix(x)) +
              sapply(matrices, function (x) is.list(x))) > 1 &
      (any((sapply(matrices, function (x) is.matrix(x)) +
            sapply(matrices, function (x) is.list(x))) == 1 )) &
      (!all((sapply(matrices, function (x) is.matrix(x)) +
            sapply(matrices, function (x) is.list(x))) == 1 )))
      stop("Enter only individual matrices or one list of matrices1.")
  if (length(matrices) == 1 & is.null(dim(matrices[[1]])))
    matrices[[1]] <- as.matrix(matrices[[1]])
  if (noname <- !any(unlist(lapply(lapply(matrices, colnames), is.null))))
    col_names <- unlist(lapply(matrices, colnames))
  dimensions <- sapply(matrices, dim)
  total_rows <- sum(dimensions[1, ])
  total_cols <- sum(dimensions[2, ])
  result <- matrix(fill, nrow = total_rows, ncol = total_cols)
  current_row <- current_col <- 1
  for (i in seq_along(matrices)) {
    rows <- dimensions[1, i]
    cols <- dimensions[2, i]
    result[current_row:(current_row + rows - 1),
           current_col:(current_col + cols - 1)] <- as.matrix(matrices[[i]])
    current_row <- current_row + rows
    current_col <- current_col + cols
  }
  if (noname) colnames(result) <- rownames(result) <- col_names
  if (is.character(result))
    result <- matrix(as.numeric(result), nrow = nrow(result), ncol = ncol(result))
  return(result)
}
