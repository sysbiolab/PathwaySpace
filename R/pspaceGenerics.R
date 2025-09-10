
setGeneric("circularProjection", 
  function(ps, pdist = 0.15, ...) {
    standardGeneric("circularProjection")
  }, package = "PathwaySpace"
)

setGeneric("polarProjection", 
  function(ps, pdist = 0.5, ...) {
    standardGeneric("polarProjection")
  }, package = "PathwaySpace"
)

setGeneric("silhouetteMapping", 
  function(ps, pdist = 0.05, ...) {
    standardGeneric("silhouetteMapping")
  }, package = "PathwaySpace"
)

setGeneric("summitMapping", 
  function(ps, maxset = 30, minsize = 30, 
    threshold = 0.5, ...) {
    standardGeneric("summitMapping")
  }, package = "PathwaySpace"
)

setGeneric("plotPathwaySpace", 
  function(ps, ...) {
    standardGeneric("plotPathwaySpace")
  }, package = "PathwaySpace"
)

setGeneric("getPathwaySpace", 
  function(ps, what = "status") {
    standardGeneric("getPathwaySpace")
  }, package = "PathwaySpace"
)

setGeneric("vertexSignal", function(x, value) standardGeneric("vertexSignal"))
setGeneric("vertexSignal<-", function(x, value) standardGeneric("vertexSignal<-"))

setGeneric("vertexDecay", function(x, value) standardGeneric("vertexDecay"))
setGeneric("vertexDecay<-", function(x, value) standardGeneric("vertexDecay<-"))


