# Tests for PathwaySpace-class constructor
test_that("Check PathwaySpace-class constructor", {
  data("gtoy1", package = "RGraphSpace")
  pts <- buildPathwaySpace(gtoy1, nrc = 30)
  expect_true(is(pts, "PathwaySpace"))
})
