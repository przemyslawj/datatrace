library(testthat)

context('Tests for calcPlaceField')

test_that('test place field correct', {
  xy =    c(1, 2, 3, 4, 5, 5, 4, 3)
  trace = c(0, 0, 0, 1, 1, 0, 1, 0)
  binnedTrace = trace + 1
  
  pf = calcPlaceField(xy, 6, trace, binnedTrace, c(5, length(trace)), 0, 2, 1)
  expect_equal(pf$field[1:6], c(0, 0, 0, 1, 0.5, NaN))
  expect_gt(pf$spatial.information, 0.5)
  expect_equal(c(1, 1, 2, 2, 2, 0), pf$occupancy[1:6])
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
  result = cell.spatial.info(df, 6, 2, min.occupancy.sec = 1, bin.hz = 1)
  
  field = result$field %>% .to.matrix(6, 2)
  expect_true(all(field[2:6,2] >= 0))
  expect_true(is.na(field[1,2]))
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
  result = cell.spatial.info(df, 6, 2, min.occupancy.sec = 1, bin.hz = 2, generate.plots = TRUE)
  
  field = result$field %>% .to.matrix(6, 2)
  expect_true(all(field[4:6,2] >= 0))
  expect_true(all(is.na(field[1:3,2])))
  expect_gt(result$cell_info$mutual.info, 0.6)
  expect_gt(result$cell_info$spatial.information, 0.5)
})
