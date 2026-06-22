
# Tests for circularProjection and polarProjection methods
#
# These codify behavior that was previously verified only by ad-hoc, manual
# checks -- in particular, a beta=0 topology-leak regression was once
# silently reintroduced during a refactor and only caught by manual review.
# These tests exist so that kind of regression is caught automatically.

test_that("circularProjection produces a valid projection", {
  data("gtoy1", package = "RGraphSpace")
  ps <- buildPathwaySpace(gtoy1, nrc = 100)
  vertexSignal(ps) <- 1
  ps <- circularProjection(ps, verbose = FALSE)
  expect_true(is(ps, "PathwaySpace"))
  expect_true(sum(ps@projection@signal) > 0)
})

test_that("polarProjection produces a valid projection", {
  data("gtoy2", package = "RGraphSpace")
  ps <- buildPathwaySpace(gtoy2, nrc = 200)
  vertexSignal(ps) <- 1
  ps <- polarProjection(ps, beta = 5, verbose = FALSE)
  expect_true(is(ps, "PathwaySpace"))
  expect_true(sum(ps@projection@signal) > 0)
})

test_that("polarProjection at beta=0 is independent of node topology", {
  # At beta=0, every node should project an identical, signal-only circle,
  # regardless of its degree or its incident edge lengths. This was the
  # original topology-leak bug: edge length scaling was not neutralized at
  # beta=0, so nodes with different degrees produced different-sized
  # circles even with identical signal.
  data("gtoy2", package = "RGraphSpace")
  ps <- buildPathwaySpace(gtoy2, nrc = 200)
  vertexSignal(ps) <- 1
  ps <- polarProjection(ps, beta = 0, verbose = FALSE)
  expect_equal(sum(ps@projection@signal), 378.405, tolerance = 1e-3)
})

test_that("isolated nodes project identically in directional and non-directional modes", {
  # Isolated nodes (no edges) are deliberately given the same area-matched
  # treatment regardless of 'directional', so this sum should be identical
  # to the non-directional beta=0 case above.
  data("gtoy2", package = "RGraphSpace")
  ps <- buildPathwaySpace(gtoy2, nrc = 200)
  vertexSignal(ps) <- 1
  ps_dir <- polarProjection(ps, beta = 0, directional = TRUE, verbose = FALSE)
  expect_equal(sum(ps_dir@projection@signal), 378.405, tolerance = 1e-3)
})

test_that("polarProjection aggregate.fun methods give expected results", {
  data("gtoy2", package = "RGraphSpace")
  ps <- buildPathwaySpace(gtoy2, nrc = 200)
  vertexSignal(ps) <- c(1, 2, 0.5, 3, 1.5, 0.8)
  
  expected <- c(
    mean      = 8.706840,
    wmean     = 17.413679,
    log.wmean = 17.413679,
    exp.wmean = 17.413679
  )
  for (m in names(expected)) {
    z <- polarProjection(ps, beta = 5, aggregate.fun = signalAggregation(m),
      verbose = FALSE)
    expect_equal(sum(z@projection@signal), expected[[m]], tolerance = 1e-4,
      label = paste("polarProjection with aggregate.fun =", m))
  }
})

test_that("circularProjection aggregate.fun methods give expected results", {
  data("gtoy2", package = "RGraphSpace")
  ps <- buildPathwaySpace(gtoy2, nrc = 200)
  vertexSignal(ps) <- c(1, 2, 0.5, 3, 1.5, 0.8)
  
  expected <- c(
    mean      = 178.607140,
    wmean     = 357.046364,
    log.wmean = 357.046355,
    exp.wmean = 357.048541
  )
  for (m in names(expected)) {
    z <- circularProjection(ps, aggregate.fun = signalAggregation(m),
      verbose = FALSE)
    expect_equal(sum(z@projection@signal), expected[[m]], tolerance = 1e-4,
      label = paste("circularProjection with aggregate.fun =", m))
  }
})

test_that("polarProjection falls back correctly for k-truncation and custom aggregators", {
  # k < vcount and non-default aggregate.fun both bypass the fast
  # aggregation path; this confirms that fallback path still runs and
  # produces a sane, non-degenerate result rather than erroring out.
  data("gtoy2", package = "RGraphSpace")
  ps <- buildPathwaySpace(gtoy2, nrc = 200)
  vertexSignal(ps) <- c(1, 2, 0.5, 3, 1.5, 0.8)
  
  z_k <- polarProjection(ps, beta = 5, k = 2, verbose = FALSE)
  expect_true(sum(z_k@projection@signal) > 0)
  
  z_custom <- polarProjection(ps, beta = 5,
    aggregate.fun = function(x) median(x), verbose = FALSE)
  expect_true(sum(z_custom@projection@signal) > 0)
})

test_that("signalAggregation() defaults to the 'mean' method", {
  f <- signalAggregation()
  expect_equal(attr(f, "name"), "meanSignal")
})
