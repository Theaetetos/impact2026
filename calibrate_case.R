# Calibrate the simulated self-funded case and produce the Part 3
# deck numbers.
#
# Run from the project root after devtools::load_all(). Writes
# results to data-raw/output/. Takes a few minutes on one core; set
# `reps` lower to iterate.

devtools::load_all()

reps    <- 1000
params  <- case_params(n_post = 36L)   # 3 years post-launch
cutoffs <- c(2L, 11L, 35L)       # analysis at months 3, 12, 36

run_one <- function(seed) {
  print(seed)
  print(Sys.time())
  s <- simulate_case(params, seed = seed)
  out <- do.call(
    rbind,
    lapply(cutoffs, function(k) estimate_case(s, cutoff = k))
  )
  # estimate_case() returns NULL below min_treated (see
  # methodological-decisions notes); when every cutoff for this seed
  # is skipped, out is NULL too -- `out$seed <- seed` on a bare NULL
  # would silently create a one-element list instead of a 0-row
  # frame, which then breaks the do.call(rbind, ...) below on
  # whichever seed hits this first. Drop the replication instead.
  if (is.null(out)) return(NULL)
  out$seed <- seed
  out
}

mc <- do.call(rbind, lapply(seq_len(reps), run_one))
dir.create("data-raw/output", showWarnings = FALSE)
saveRDS(mc, "data-raw/output/monte_carlo.rds")

# ----- Operating characteristics by method and analysis month -----
truth <- mc[mc$method == "Truth", c("seed", "cutoff", "estimate")]
names(truth)[3] <- "truth"
mc2 <- merge(mc[mc$method != "Truth", ], truth)
oc <- do.call(
  rbind,
  lapply(
    split(mc2, list(mc2$method, mc2$cutoff)),
    function(d){
      # Only est_matched() ever sets skip_reason (see its roxygen
      # note); for every other method this is a no-op filter. bias/
      # sd/coverage are computed on successful rows only -- n plus
      # the two skip counts is how many replications were attempted
      # in total (reps, less any all-methods NULL at n_treated == 0,
      # which never reach mc2 at all), so a reader can see exactly
      # how much of the nominal replication count each number rests
      # on rather than have skips silently thin the denominator.
      ok <- is.na(d$skip_reason)
      # A handful of "successful" (skip_reason NA) rows for Pre/post
      # and Difference-in-differences still come back with a
      # non-finite CI: at n_treated == 1 there are zero residual
      # degrees of freedom for the cluster-robust/HC1 SE, so
      # lower/upper are NaN even though est_prepost()/est_did() never
      # flagged anything (they have no skip_reason mechanism -- that
      # belongs to est_matched() alone; see its roxygen note). Left
      # unguarded, a single NaN poisons mean() for the whole
      # method x cutoff cell, silently NAing out coverage/ci_width
      # built from 300+ otherwise-good rows. This is a genuine
      # analysis-scope fix (don't let one bad row erase everything),
      # not a new skip-tracking system for these two methods --
      # deliberately not building a skip_reason equivalent for
      # Pre/post or DID here. undefined_ci counts how many of the
      # "successful" rows this happened to, so it's visible rather
      # than silently dropped.
      ci_ok <- ok & is.finite(d$lower) & is.finite(d$upper)
      data.frame(
        method = d$method[1],
        month = d$cutoff[1] + 1,
        n = sum(ok),
        n_skipped_too_few_engaged =
          sum(d$skip_reason == "too_few_engaged", na.rm = TRUE),
        n_skipped_implausible_fit =
          sum(d$skip_reason == "implausible_fit", na.rm = TRUE),
        n_ci_undefined = sum(ok) - sum(ci_ok),
        truth = mean(d$truth[ok], na.rm = TRUE),
        mean_estimate = mean(d$estimate[ok], na.rm = TRUE),
        bias = mean(d$estimate[ok] - d$truth[ok], na.rm = TRUE),
        sd = sd(d$estimate[ok], na.rm = TRUE),
        ci_width = mean(
          d$upper[ci_ok] - d$lower[ci_ok],
          na.rm = TRUE
        ),
        coverage = mean(
          d$lower[ci_ok] <= d$truth[ci_ok] &
            d$truth[ci_ok] <= d$upper[ci_ok],
          na.rm = TRUE
        )
      )
    }
  )
)
print(oc, digits = 3)
write.csv(
  oc,
  "data-raw/output/operating_characteristics.csv",
  row.names = FALSE
)

# - Forest plot: one typical realisation vs. the Monte Carlo spread -
# Pick the "case" seed as the realisation whose estimates are jointly
# closest to the Monte Carlo means (a typical draw, not the most
# flattering one), and disclose the seed and the selection rule on
# the slide.
m12 <- mc2[mc2$cutoff == 11L, ]   # typical at month 12
# A skipped matching row (rare by month 36, but recorded rather than
# dropped -- see est_matched()) has an NA estimate; leaving it in
# would NA out that method's mean/sd for every seed via tapply(), not
# just the one seed that hit the skip. Excluding the seed entirely
# (rather than just its NA row) avoids a subtler bias: a seed missing
# one method's term would otherwise be summing fewer squared terms in
# `dist` below and could win `which.min()` for that reason alone.
incomplete <- unique(m12$seed[!is.na(m12$skip_reason)])
m12 <- m12[!(m12$seed %in% incomplete), ]
means <- tapply(m12$estimate, m12$method, mean)
sds   <- tapply(m12$estimate, m12$method, sd)
dist <- tapply(
  seq_len(nrow(m12)),
  m12$seed,
  function(i){
    sum(
      ((m12$estimate[i] - means[m12$method[i]])
       / sds[m12$method[i]])^2
    )
  }
)
case_seed <- as.integer(names(which.min(dist)))
message("Typical case seed: ", case_seed)

forest <- function(k, file, title) {
  case <- mc[mc$seed == case_seed & mc$cutoff == k, ]
  ord <- rev(
    c("Pre/post", "ITS", "CITS", "DiD", "P+G")
  )
  case <- case[match(ord, case$method), ]
  tru <- mc$estimate[
    mc$seed == case_seed & mc$cutoff == k & mc$method == "Truth"
  ]

  # case_seed is only guaranteed skip-free at cutoff 35 (that's what
  # it was chosen on); at month 3 the same seed's matching row may
  # well be a recorded skip (NA estimate/lower/upper -- likely, since
  # month 3 is exactly where n_treated < 20 is common). range() and
  # plot()/segments() below all need na.rm/NA-tolerant handling so
  # that shows up as a missing point in this panel rather than a
  # crash; note in words if it happened so nobody's stuck reading the
  # skip out of raw pixel positions.
  skipped <- case$method[!is.na(case$skip_reason)]
  if (length(skipped)) {
    message(
      "forest(",
      k,
      "): case_seed ",
      case_seed,
      " has a recorded skip for: ",
      paste(skipped, collapse = ", "),
      " (", case$skip_reason[!is.na(case$skip_reason)], ")"
    )
  }

  grDevices::png(file, 1600, 900, res = 200)
  op <- graphics::par(mar = c(4, 13, 2, 1))
  xl <- range(c(case$lower, case$upper, 0, 1.25), na.rm = TRUE)
  graphics::plot(
    case$estimate,
    seq_along(ord),
    xlim = xl,
    yaxt = "n",
    pch = 19,
    xlab = "Savings, PMPM", ylab = "", main = title
  )
  graphics::segments(
    case$lower,
    seq_along(ord),
    case$upper,
    seq_along(ord)
  )
  graphics::axis(2, seq_along(ord), ord, las = 1)
  graphics::abline(v = tru, lty = 2, col = "firebrick")
  graphics::abline(v = 1.25, lty = 3, col = "grey40")
  graphics::abline(v = 0, col = "grey80")
  graphics::par(op)
  grDevices::dev.off()
}
forest(
  2L,
  "data-raw/output/forest_month03.png",
  "Month 3, one simulated case"
)
forest(
  11L,
  "data-raw/output/forest_month12.png",
  "Month 12, one simulated case"
)
forest(
  35L,
  "data-raw/output/forest_month36.png",
  "Month 36, one simulated case"
)
