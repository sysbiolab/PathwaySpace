# Tests for PathwaySpace-class methods
test_that("Check PathwaySpace-class methods", {
  data("gtoy1", package = "PathwaySpace")
  pts <- buildPathwaySpace(gtoy1, nrc = 30)
  expect_true(is(pts, "PathwaySpace"))
})
