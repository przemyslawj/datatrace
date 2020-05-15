library(testthat)

context('Tests for calcPlaceField')

test_that('test place field correct', {
  x =     c(1, 2, 3, 4, 5, 5, 4, 3)
  y =     c(1, 1, 1, 1, 1, 1, 1, 1)
  trace = c(0, 0, 0, 1, 1, 0, 1, 0)
  binnedTrace = trace + 1
  
  pf = calcPlaceField(bin_x=x, bin_y=y, nbins_x=6, nbins_y=2, 
                      trace=trace, binnedTrace=binnedTrace, minOccupancy=1, 
                      kernelSize=0, gaussianVar=0)
  expect_equal(pf$field[1:6, 1], c(0, 0, 0, 1, 0.5, NA))
  expect_gt(pf$spatial.information, 0.8)
  expect_equal(c(1, 1, 2, 2, 2, 0), pf$occupancy[1:6, 1])
})
  
test_that('test place field correct and smoothed', {
  x =     c(1, 2, 3, 4, 5, 5, 4, 3)
  y =     c(1, 1, 1, 1, 1, 1, 1, 1)
  trace = c(0, 0, 0, 1, 0, 0, 1, 0)
  binnedTrace = trace + 1
  
  pf = calcPlaceField(bin_x=x, bin_y=y, nbins_x=6, nbins_y=2, 
                      trace=trace, binnedTrace=binnedTrace, minOccupancy=1, 
                      kernelSize=3, gaussianVar=1)
  expect_equal(pf$field[1:2, 1], c(0, 0))
  expect_true(all(is.na(pf$field[, 2])))
  expect_gt(pf$field[3, 1], 0.31)
  expect_gt(pf$field[4, 1], 0.45)
  expect_gt(pf$field[5, 1], 0.37)
  expect_gt(pf$spatial.information, 0.2)
})

context('Tests calculated place field metrics ')
test_that('test cell.spatial.info succeeds', {
  df = data.frame(
    cell_id = 0,
    bin.y = 1,
    bin.x = c(1, 2, 3, 4, 5, 5, 4, 3),
    trace = c(0, 0, 0, 1, 1, 0, 1, 0)
  )
  df$response_bin = df$trace + 1
  df$timestamp = seq_along(df$trace)
  df$nevents = df$trace
  result = cell.spatial.info(df, 6, 2, min.occupancy.sec=1, bin.hz=1, kernel.size=0, gaussian.var=0)
  
  expect_true(all(result$field[1:5,1] >= 0))
  expect_true(is.na(result$field[6,1]))
  expect_gt(result$cell_info$mutual.info, 0.6)
  expect_gt(result$cell_info$spatial.information, 0.5)
})

test_that('test cell.spatial.info discards unoccupied', {
  df = data.frame(
    cell_id = 0,
    bin.y = 1,
    bin.x = c(1, 2, 3, 4, 5, 5, 4, 3),
    trace = c(0, 0, 0, 1, 1, 0, 1, 0)
  )
  df$response_bin = df$trace + 1
  df$timestamp = seq_along(df$trace)
  df$nevents = df$trace
  result = cell.spatial.info(df, 6, 2, min.occupancy.sec = 1, bin.hz = 2, generate.plots = TRUE,
                             kernel.size=0, gaussian.var=0)
  
  expect_true(all(result$field[3:5,1] >= 0))
  expect_true(all(is.na(result$field[1:2,1])))
  expect_gt(result$cell_info$mutual.info, 0.5)
  expect_gt(result$cell_info$spatial.information, 1.0)
})
