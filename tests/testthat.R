# This file is part of the standard setup for testthat.
# It is recommended that you do not modify it.
#
# Where should you do additional test configuration?
# Learn more about the roles of various files in:
# * https://r-pkgs.org/tests.html
# * https://testthat.r-lib.org/reference/test_package.html#special-files

library(testthat)
library(ecurve)

# Setting the seed to a fixed value
# to debug potential issues with some tests.
# The `nlm()` optimization occasionally fails
# to find a good enough MLE, which fails the test.
set.seed(1234)

test_check("ecurve")
