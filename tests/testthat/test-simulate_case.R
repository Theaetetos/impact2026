test_that("case lands on the deck's scale", {
  s <- simulate_case(case_params(), seed = 1)
  eng <- s$members$engaged_month
  expect_equal(sum(eng <= 2, na.rm = TRUE), 25, tolerance = 0.2)   # slide 24
  expect_equal(sum(!is.na(eng)), 102, tolerance = 0.1)             # ~0.6% of 17k
  expect_true(abs(sum(s$members$buying) - 17000) < 2000)
})

test_that("truth and matching agree on average", {
  skip_on_cran()
  est <- do.call(rbind, lapply(1:20, function(sd)
    estimate_case(simulate_case(case_params(n_post = 36L), seed = sd), 35L)))
  truth <- est$estimate[est$method == "Truth"]
  pg <- est$estimate[est$method == "Matching + G-computation"]
  expect_lt(abs(mean(pg - truth)), 3 * sd(pg - truth) / sqrt(20))
  pp <- est$estimate[est$method == "Pre/post"]
  expect_gt(mean(pp - truth), 3)            # regression to the mean
  did <- est$estimate[est$method == "Difference-in-differences"]
  expect_lt(mean(did - truth), -1)          # prioritisation bias
})
