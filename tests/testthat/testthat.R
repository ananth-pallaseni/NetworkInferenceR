library(testthat)
library(NetworkInferenceR)

context("NetworkInferenceR")


datadir <- "../data/"
datapath <- paste(datadir, "yeast1_10_data.txt", sep = "")

mi_net <- infer_network(datapath, "MI")
mi_truth <- read.table(paste(datadir, "mi.txt", sep = ""))
test_that("MI is correct", { for (i in c(1, 5, 10, 20, 40)) {
 expect_equal(mi_net[i, 3], mi_truth[i * 2, 3])
}
})

puc_net <- infer_network(datapath, "PUC")
puc_truth <- read.table(paste(datadir, "puc.txt", sep = ""))
test_that("PUC is correct", { for (i in c(1, 5, 10, 20, 40)) {
  expect_equal(puc_net[i, 3], puc_truth[i * 2, 3])
}
})

pid_net <- infer_network(datapath, "PID")
pid_truth <- read.table(paste(datadir, "pidc.txt", sep = ""))
test_that("PID is correct", { for (i in c(1, 5, 10, 20, 40)) {
  expect_equal(pid_net[i, 3], pid_truth[i * 2, 3])
}
})
