# Tests for PathwaySpace-class methods
test_that("Check PathwaySpace-class methods", {
  data("gtoy1", package = "RGraphSpace")
  pts <- buildPathwaySpace(gtoy1, nrc = 30)
  expect_true(is(pts, "PathwaySpace"))
})
