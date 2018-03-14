
infer_network <- function(filepath, method="PID") {
  data <- t(read.table(filepath))
  if (method == "PID") {
    network <- infer_pid_network(data)
  }
  else if (method == "PUC") {
    network <- infer_puc_network(data)
  }
  else if (method == "MI") {
    network <- infer_mi_network(data)
  }
  else {
    errstr <- paste("Unrecognised method argument:", method, "Should be one of: PID, PUC or MI")
    print(errstr)
  }

  return(network[order(-network[,3]),])
}
