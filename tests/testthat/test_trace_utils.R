context('Test for chunk shuffle')
test_that('chunkShuffle succeeds', {
  trace = 1:10
  trialEnds = c(6, 12)
  res = chunkShuffle(trace, trialEnds, 2)
  expect_true(all(res <= 10))
  expect_true(all(res >= 1))
  expect_true(all(res[1:6] <= 6))
  expect_true(all(res[7:10] >= 7))
  expect_true(all(diff(res)[c(1,3,5,7,9)] == 1))
})

context('Test for randomShift')
test_that('randomShift succeeds', {
  trace = 1:10
  trialEnds = c(6, 10)
  res = randomShift(trace, trialEnds, 2)
  expect_true(all(res[1:6] <= 6))
  expect_true(all(res[7:10] >= 7))
  expect_true(all(abs(res[1:6] - trace[1:6]) >= 2))
  expect_true(all(abs(res[7:10] - trace[7:10]) >= 1))
})

context('Tests for the filter.running')

test_that('filter.running succeeds for one trial', {
  running.vel.thr = 2
  epoch.dur.ms = 100
  fr = 20
  running.pos = seq(0, by=running.vel.thr*8 / fr, length.out=4)
  df = data.frame(
    id=1:20,
    trial_id = '1',
    cell_id = 1,
    exp_title = 'trial',
    timestamp = seq(0, by=1000/fr, length.out=20),
    x = c(running.pos, rep(4.0, 6), running.pos + 2.0, rep(10.0, 6)),
    y = 1,
    velocity = c(rep(running.vel.thr * 4, 5), rep(0, 5), rep(running.vel.thr * 4, 5), rep(0, 5))
  )
  running.index = isRunning(df, 2, running.vel.thr, epoch.dur.ms)
  #df.filtered = filter.running(df, min.run.velocity=2, mean.run.velocity=running.vel.thr, window.dur.ms=epoch.dur.ms)
  df.filtered = df[which(running.index),]
  expect_equal(nrow(df.filtered), 10)
  expect_equal(df.filtered$id, c(1:5, 11:15))
}) 

test_that('test detect.events succeeds', {
  df = data.frame(deconv_trace=c(1:10, c(0, 12:20)), 
                  animal='a', 
                  cell_id=rep(1:2,each=10), 
                  date='2019-01-01')
  df = data.table(df)
  detect.events(df, deconv.threshold=0.11)
  expect_equal(df$is.event[1], FALSE)
  expect_equal(df$is.event[11], FALSE)
  expect_equal(sum(df$is.event), 18)
})

test_that('bin.time.space succeeds', {
  df = data.frame(
    trace=1:10,
    timestamp=seq(0, 900, 100),
    x=seq(0, 90, 10),
    y=seq(90, 0, -10),
    trial=1,
    trial_id='test_1',
    animal='A',
    cell_id=1,
    exp_title='trial',
    is.event=0,
    dist=0,
    date='2019-01-01'
  )
  res = bin.time.space(df, 5, 5,
                       bin.quantile.fractions = c(0.5, 1.0), 
                       binned.var='trace', 
                       timebin.dur.msec = 200)
  expect_equal(nrow(res), 5)
  expect_equal(res$time_bin, 0:4)
  expect_equal(res$timestamp, seq(50, 850, 200))
  expect_equal(res$bin.x, 0:4)
  expect_equal(res$bin.y, 4:0)
  expect_equal(res$bin.xy, c(5, 9, 13, 17, 21))
  expect_equal(res$trace, seq(1.5, 9.5, 2))
})