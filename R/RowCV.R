RowCV <- function(x) {
  sqrt((rowSums(((x - rowMeans(x))^2)/(dim(x)[2] - 1)))/(rowMeans(x)+1e-10) * 100)
}