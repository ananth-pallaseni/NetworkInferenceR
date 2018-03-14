
infer_network <- function(filepath, method="PIDC") {
  data <- t(read.table(filepath))
  if (method == "PIDC") {
    network <- infer_pidc_network(data)
  }
  else if (method == "PUC") {
    network <- infer_puc_network(data)
  }
  else if (method == "MI") {
    network <- infer_mi_network(data)
  }
  else {
    errstr <- paste("Unrecognised method argument:", method, "Should be one of: PIDC, PUC or MI")
    print(errstr)
  }

  return(network[order(-network[,3]),])
}
