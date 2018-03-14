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

pidc_net <- infer_network(datapath, "PIDC")
pidc_truth <- read.table(paste(datadir, "pidc.txt", sep = ""))
test_that("PIDC is correct", { for (i in c(1, 5, 10, 20, 40)) {
  expect_equal(pidc_net[i, 3], pidc_truth[i * 2, 3])
}
})
