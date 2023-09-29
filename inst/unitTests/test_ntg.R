# Unit tests fot PathwaySpace-class methods
test_pts <- function(){
  data("gtoy1", package = "PathwaySpace")
  pts <- buildPathwaySpace(gtoy1, nrc = 30)
  checkTrue(is(pts, "PathwaySpace"))
}
