
setGeneric("circularProjection", 
    function(pts, kns = 8, pdist = 0.15, rescale = TRUE,
        verbose = TRUE, decay_fun = weibullDecay, 
        knn = deprecated(), ...) {
        standardGeneric("circularProjection")
    }, package = "PathwaySpace"
)

setGeneric("polarProjection", 
    function(pts, kns = 8, pdist = 0.5, rescale = TRUE, 
        theta = 180, directional = FALSE, verbose = TRUE, 
        decay_fun = weibullDecay, 
        knn = deprecated(), ...) {
        standardGeneric("polarProjection")
    }, package = "PathwaySpace"
)

setGeneric("silhouetteMapping", 
    function(pts, baseline = 0.01, pdist = 0.05, verbose = TRUE) {
        standardGeneric("silhouetteMapping")
    }, package = "PathwaySpace"
)

setGeneric("summitMapping", 
    function(pts, maxset = 30, minsize = 30, threshold = 0.5, 
        verbose = TRUE, segm_fun = summitWatershed, ...) {
        standardGeneric("summitMapping")
    }, package = "PathwaySpace"
)

setGeneric("plotPathwaySpace", 
    function(pts, colors = pspace.cols(), trim.colors = c(3, 2, 1, 2, 3), 
        bg.color = "grey85", theme.name = c("th0", "th1", "th2", "th3"),
        title = "PathwaySpace", font.size = 1, font.color = "white",
        xlab = "Pathway coordinates 1", ylab = "Pathway coordinates 2", 
        zlab = "Density", zlim = NULL, slices = 25, add.grid = TRUE, 
        grid.color = "white", add.contour = TRUE, contour.color = "white", 
        label.summits = TRUE, marks = FALSE, mark.size = 3, 
        mark.color = "white", mark.padding = 0.5, 
        mark.line.width = 0.5, use.dotmark = FALSE) {
        standardGeneric("plotPathwaySpace")
    }, package = "PathwaySpace"
)

setGeneric("getPathwaySpace", 
    function(pts, what = "status") {
        standardGeneric("getPathwaySpace")
    }, package = "PathwaySpace"
)

setGeneric("vertexSignal", function(pts, value) standardGeneric("vertexSignal"))
setGeneric("vertexSignal<-", function(pts, value) standardGeneric("vertexSignal<-"))

setGeneric("vertexWeight", function(pts, value) standardGeneric("vertexWeight"))
setGeneric("vertexWeight<-", function(pts, value) standardGeneric("vertexWeight<-"))


