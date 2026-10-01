# All estimators report savings in the deck's PMPM convention (slide
# 15): savings per engaged member-year, divided by block
# member-months per engaged member at the design engagement rate
# (about 2,000). Positive = savings.
#
# Estimators only see what an analyst would see: they never read
# columns of `members` that start with "." or the `true_savings`
# matrix.

pmpm_divisor <- function(sim) {
  p <- sim$params
  sum(sim$members$buying) * 12 / (p$slots_per_month * 12)
}

# Mean monthly cost for each member over a window of columns
window_mean <- function(cost, rows, from, to) {
  vapply(seq_along(rows), function(k) {
    if (to[k] < from[k]) return(NA_real_)
    mean(cost[rows[k], from[k]:to[k]])
  }, numeric(1))
}

result_row <- function(method, cutoff, est_month, se_month,
                       n_treated, div, crit = 1.96) {
  data.frame(
    method = method, cutoff = cutoff, n_treated = n_treated,
    estimate = est_month * 12 / div,
    lower = (est_month - crit * se_month) * 12 / div,
    upper = (est_month + crit * se_month) * 12 / div,
    # NA for every estimator except est_matched(), which can
    # populate it with why it was skipped instead of returning a row
    # -- kept here (rather than added only in est_matched()) so
    # every row, skipped or not, has identical columns and rbind()
    # never has to special-case one estimator's shape.
    skip_reason = NA_character_
  )
}

# A skipped Matching + G-computation attempt still returns a row --
# NA estimate/lower/upper, tagged with why -- rather than vanishing,
# so calibrate_case.R can report both kinds of skip explicitly instead
# of the row (and the information that it was ever attempted) just
# disappearing from monte_carlo.rds.
skip_row_matched <- function(cutoff, n_treated, reason) {
  data.frame(
    method = "P+G", cutoff = cutoff,
    n_treated = n_treated, estimate = NA_real_, lower = NA_real_,
    upper = NA_real_, skip_reason = reason
  )
}

#' True savings the estimators are aiming at expected savings per
#' engaged member-month over the observed post-engagement months
#' through `cutoff`, in PMPM. Uses the latent counterfactual rates.
#' @export
true_effect <- function(sim, cutoff = 11L) {
  m <- sim$members; mo <- sim$months
  e <- which(!is.na(m$engaged_month) & m$engaged_month <= cutoff)
  cols <- which(mo >= 0 & mo <= cutoff)
  ts <- sim$true_savings[e, cols, drop = FALSE]
  eng <- outer(m$engaged_month[e], mo[cols], `<=`)
  result_row(
    "Truth",
    cutoff,
    sum(ts[eng]) / sum(eng),
    0,
    length(e),
    pmpm_divisor(sim)
  )
}

# Engaged members and the analytic frame: index month, pre and post
# windows.
#
# treat is written as rep(1L, length(e)) rather than a bare 1L:
# data.frame() recycles a length-1 column against columns of length
# N > 0, but NOT against columns of length 0 -- "arguments imply
# differing number of rows: 0, 1" -- so with zero engaged members
# (e empty) this used to error before estimate_case() even got a
# chance to see n_treated and skip cleanly. Same fix in
# comparison_frame() below for the same reason (its "nonbuying" pool
# for est_did() can be empty even when engaged members exist).
treated_frame <- function(sim, cutoff) {
  m <- sim$members; mo <- sim$months
  e <- which(!is.na(m$engaged_month) & m$engaged_month <= cutoff)
  idx <- match(m$engaged_month[e], mo)
  data.frame(row = e, treat = rep(1L, length(e)), index_col = idx)
}

comparison_frame <- function(sim, cutoff, pool) {
  m <- sim$members
  mo <- sim$months
  keep <- !is.na(m$eligible_month) &
    m$eligible_month <= cutoff &
    (is.na(m$engaged_month) | m$engaged_month > cutoff)
  keep <- keep & switch(
    pool,
    nonbuying = !m$buying,
    all = TRUE
  )
  c_rows <- which(keep)
  # Pseudo-enrollment date: engaged members wait in the capacity
  # queue between trigger and engagement, and regress toward the mean
  # while they wait. Give each comparison member a delay drawn from
  # the engaged members' delays so both groups are indexed the same
  # distance from their trigger spike.
  e <- !is.na(m$engaged_month) & m$engaged_month <= cutoff
  delays <- m$engaged_month[e] - m$eligible_month[e]
  # if/else, not ifelse(): the test is a scalar, and ifelse() would
  # return only the first sampled delay, recycled to every
  # comparison member.
  d <- if (length(delays)) {
    delays[
      sample.int(length(delays), length(c_rows), replace = TRUE)
    ]
  } else {
    0L
  }
  idx_month <- m$eligible_month[c_rows] + d
  # as with engaged members, not yet enrolled -> out
  ok <- idx_month <= cutoff
  data.frame(
    row = c_rows[ok],
    treat = rep(0L, sum(ok)),
    index_col = match(idx_month[ok], mo)
  )
}

add_windows <- function(sim, fr, cutoff) {
  m <- sim$members
  mo <- sim$months
  cut_col <- match(cutoff, mo)
  fr$pre  <- window_mean(
    sim$cost,
    fr$row,
    fr$index_col - 12L,
    fr$index_col - 1L
  )
  fr$post <- window_mean(
    sim$cost,
    fr$row,
    fr$index_col,
    rep(cut_col, nrow(fr))
  )
  fr$n_post <- cut_col - fr$index_col + 1L
  fr$pre_ed <- vapply(
    seq_len(nrow(fr)),
    function(k)
    sum(
      sim$counts$ed[
        fr$row[k],
        (fr$index_col[k] - 12L):(fr$index_col[k] - 1L)
      ]
    ),
    numeric(1)
  )
  fr$group <- m$group[fr$row]
  fr$age <- m$age[fr$row]
  fr$female <- m$female[fr$row]
  fr$conditions <- m$conditions[fr$row]
  fr$risk_score <- m$risk_score[fr$row]
  fr$index_month <- mo[fr$index_col]
  # Each group's pre-launch PMPM: observable, and it carries the
  # group cost level that made some groups more likely to buy.
  g_pre <- tapply(
    rowMeans(sim$cost[, mo < 0, drop = FALSE]),
    m$group,
    mean
  )
  fr$group_pre_pmpm <- as.numeric(g_pre[as.character(fr$group)])
  fr
}

#' Pre/post on engaged members only
#' @export
est_prepost <- function(sim, cutoff = 11L) {
  fr <- add_windows(sim, treated_frame(sim, cutoff), cutoff)
  fit <- stats::lm(I(pre - post) ~ 1, data = fr, weights = n_post)
  result_row(
    "Pre/post",
    cutoff,
    stats::coef(fit)[[1]],
    sqrt(sandwich::vcovHC(fit, type = "HC1")[1, 1]),
    nrow(fr),
    pmpm_divisor(sim)
  )
}

#' Interrupted time series on the buying block's monthly PMPM
#'
#' Level shift at launch with a linear trend and seasonal terms,
#' Newey-West errors. The block-level shift is converted to the
#' per-engaged convention.
#' @export
est_its <- function(sim, cutoff = 11L) {
  mo <- sim$months
  cols <- which(mo <= cutoff)
  b <- sim$members$buying
  d <- data.frame(t = mo[cols], pmpm = colMeans(sim$cost[b, cols]))
  d$post <- as.integer(d$t >= 0)
  d$s1 <- sin(2 * pi * d$t / 12); d$c1 <- cos(2 * pi * d$t / 12)
  fit <- stats::lm(pmpm ~ t + post + s1 + c1, data = d)
  shift <- -stats::coef(fit)[["post"]]
  se <- sqrt(
    sandwich::NeweyWest(
      fit,
      lag = 3,
      prewhite = FALSE
    )["post", "post"]
  )
  # Block PMPM shift -> per engaged member-month: scale by block
  # members per engaged member-month observed through the cutoff
  n_b <- sum(b)
  eng <- sim$members$engaged_month
  eng_mm <- sum(
    pmax(0, cutoff - eng[b & !is.na(eng) & eng <= cutoff] + 1)
  )
  post_months <- sum(d$post)
  k <- n_b * post_months / max(eng_mm, 1)
  # Newey-West is anti-conservative with few post-launch points; use
  # a t critical value on the post-period degrees of freedom.
  result_row(
    "ITS",
    cutoff,
    shift * k,
    se * k,
    sum(!is.na(eng) & eng <= cutoff),
    pmpm_divisor(sim),
    crit = stats::qt(0.975, df = max(1, sum(d$post) - 1))
  )
}

#' Comparative interrupted time series: buying minus non-buying block
#' PMPM
#'
#' Same segmented model as [est_its()] on the monthly difference
#' between the buying and non-buying groups, which removes anything
#' that changed for the whole self-funded block at launch.
#' @export
est_cits <- function(sim, cutoff = 11L) {
  mo <- sim$months
  cols <- which(mo <= cutoff)
  b <- sim$members$buying
  d <- data.frame(
    t = mo[cols],
    gap = colMeans(sim$cost[b, cols]) - colMeans(sim$cost[!b, cols])
  )
  d$post <- as.integer(d$t >= 0)
  d$s1 <- sin(2 * pi * d$t / 12); d$c1 <- cos(2 * pi * d$t / 12)
  fit <- stats::lm(gap ~ t + post + s1 + c1, data = d)
  shift <- -stats::coef(fit)[["post"]]
  se <- sqrt(
    sandwich::NeweyWest(
      fit,
      lag = 3,
      prewhite = FALSE
    )["post", "post"]
  )
  n_b <- sum(b)
  eng <- sim$members$engaged_month
  eng_mm <- sum(
    pmax(0, cutoff - eng[b & !is.na(eng) & eng <= cutoff] + 1)
  )
  k <- n_b * sum(d$post) / max(eng_mm, 1)
  # Newey-West is anti-conservative with few post-launch points; use
  # a t critical value on the post-period degrees of freedom.
  result_row(
    "CITS",
    cutoff,
    shift * k,
    se * k,
    sum(!is.na(eng) & eng <= cutoff),
    pmpm_divisor(sim),
    crit = stats::qt(0.975, df = max(1, sum(d$post) - 1))
  )
}

#' Difference-in-differences: engaged vs eligible members in groups
#' that did not buy
#'
#' Returns `NULL` if the nonbuying comparison pool is empty --
#' `treat` would then have zero variation, its coefficient goes NA,
#' `vcovCL()` silently drops it, and `V["treat", "treat"]` errors on
#' the missing name. Same category as `estimate_case()`'s
#' zero-engaged guard: DID needs *some* comparison group to exist,
#' which is a hard requirement, not a chosen floor -- so this drops
#' the row via `NULL` (`estimate_case()`'s `rbind()` already
#' tolerates that) rather than inventing a skip-reason category for
#' it.
#' @export
est_did <- function(sim, cutoff = 11L) {
  fr <- rbind(
    treated_frame(sim, cutoff),
    comparison_frame(sim, cutoff, pool = "nonbuying")
  )
  if (sum(fr$treat == 0L) < 1L) return(NULL)
  fr <- add_windows(sim, fr, cutoff)
  fit <- stats::lm(
    I(pre - post) ~ treat,
    data = fr,
    weights = n_post
  )
  V <- sandwich::vcovCL(fit, cluster = ~group)
  result_row(
    "DiD",
    cutoff,
    stats::coef(fit)[["treat"]],
    sqrt(V["treat", "treat"]),
    sum(fr$treat),
    pmpm_divisor(sim)
  )
}

#' Propensity matching with G-computation
#'
#' Engaged members matched (2:1, nearest neighbour with replacement)
#' to eligible members who were not engaged by the cutoff, in any
#' group, on age, sex, conditions, risk score (the predicted cost
#' that set priority), prior cost, prior ED visits, index month, and
#' the group's pre-launch PMPM. G-computation on the matched sample
#' with a log-link (quasi-Poisson) treatment-by-covariate outcome
#' model; cluster-robust SE by the delta method.
#'
#' Below `min_treated` engaged members the outcome model's
#' covariates can't be identified (the common-slopes fallback alone
#' has 8 parameters), which surfaced as GLM non-convergence or a
#' non-conformable delta-method variance -- returns a
#' `skip_reason = "too_few_engaged"` row rather than attempting a
#' degenerate fit.
#'
#' A 200-seed sweep with the floor overridden down to 1 found that
#' below ~11 engaged, this isn't just a precision problem: 88% of
#' attempts errored outright at 1-5 engaged, and 1.9% of attempts at
#' 6-10 engaged returned a *finite, unflagged* estimate as large as
#' 3.7e45 -- MatchIt/glm() reported success (no warning, no NA
#' coefficient) on a fit that had no business succeeding. So a
#' second, independent check runs on every fit regardless of N: if
#' any fitted mu0/mu1 exceeds 10x the largest average monthly cost
#' actually observed in this replication's matched sample, or the
#' delta-method SE isn't finite, the row is tagged
#' `skip_reason = "implausible_fit"` instead of being trusted. 10x is
#' deliberately generous -- genuine blowups ran 1e6 to 1e56 past any
#' realistic ceiling -- but it means "no error was thrown" is never
#' treated as sufficient evidence of a trustworthy result on its own.
#' @export
est_matched <- function(sim, cutoff = 11L, min_treated = 20L) {
  tf <- treated_frame(sim, cutoff)
  n_treated <- nrow(tf)
  if (n_treated < min_treated) {
    return(skip_row_matched(cutoff, n_treated, "too_few_engaged"))
  }
  fit_result <- tryCatch({
    fr <- rbind(
      tf,
      comparison_frame(sim, cutoff, pool = "all")
    )
    fr <- add_windows(sim, fr, cutoff)
    fr <- fr[stats::complete.cases(fr[, c("pre", "post")]), ]
    mm <- MatchIt::matchit(
      treat ~ age + female + conditions + risk_score + pre + pre_ed
              + index_month + log(group_pre_pmpm),
      data = fr,
      method = "nearest",
      ratio = 2,
      replace = TRUE,
      estimand = "ATT"
    )
    md <- MatchIt::match.data(mm)
    md$w <- md$weights * md$n_post
    md$post <- pmax(md$post, 0)
    list(md = md)
  }, error = function(e) list(error = conditionMessage(e)))
  if (!is.null(fit_result$error)) {
    return(skip_row_matched(cutoff, n_treated, "implausible_fit"))
  }
  md <- fit_result$md
  # Log-link outcome model: costs are convex in risk, and the engaged
  # sit in the high-risk tail, so a linear model under-predicts their
  # untreated cost. Common slopes only (no treat interaction): the
  # DGP's dominant confounder (willingness) is unobserved, so a
  # treatment-by-covariate model buys protection against a
  # second-order risk (residual imbalance on *observed* covariates)
  # that matching has usually already handled, while paying for it
  # with a real one -- interacting treat with a covariate matching
  # already balanced to near-constant leaves it nearly collinear
  # with its own interaction term, and G-computation's treat=0
  # counterfactual then extrapolates through that near-singularity
  # once the interaction terms zero out and stop cancelling the
  # oversized main effects (log-link turns that linear-predictor
  # blowup into an astronomical fitted cost). Drop group_pre_pmpm
  # and index_month from this formula for the same reason -- both
  # are still in the matchit() distance model above, where they earn
  # their keep by achieving balance; once balance is achieved they
  # have ~no within-sample variance left to identify an outcome-model
  # coefficient safely. index_month specifically: once cutoff is
  # fixed it's an exact affine function of post-window length
  # (r = -1 with n_post, not an independent calendar-trend signal),
  # and its coefficient's sign flips in ~40% of replications with
  # mean |t| < 1.7 at every cutoff -- noise, not a load-bearing
  # correction, and comparison_frame() already targets balance on it
  # directly (each comparison member's pseudo-enrollment delay is
  # drawn from the engaged members' own delay distribution).
  rhs <- "age + female + conditions + risk_score + log1p(pre) + pre_ed"
  f_out <- stats::as.formula(paste("post ~ treat +", rhs))
  pred <- tryCatch({
    fit <- stats::glm(
      f_out,
      family = stats::quasipoisson(),
      data = md,
      weights = w
    )
    tr <- md[md$treat == 1, ]
    tt <- stats::delete.response(stats::terms(fit))
    X1 <- stats::model.matrix(tt, transform(tr, treat = 1L))
    X0 <- stats::model.matrix(tt, transform(tr, treat = 0L))
    b <- stats::coef(fit)
    mu1 <- exp(drop(X1 %*% b)); mu0 <- exp(drop(X0 %*% b))
    list(fit = fit, tr = tr, X1 = X1, X0 = X0, mu1 = mu1, mu0 = mu0)
  }, error = function(e) list(error = conditionMessage(e)))
  if (!is.null(pred$error)) {
    return(skip_row_matched(cutoff, n_treated, "implausible_fit"))
  }

  # Plausibility fence: see the roxygen note above for why "no error
  # was thrown" isn't sufficient on its own.
  ceiling <- 10 * max(md$post)
  implausible <- !is.finite(ceiling) || ceiling <= 0 ||
    any(!is.finite(pred$mu0)) || any(!is.finite(pred$mu1)) ||
    max(pred$mu0, pred$mu1) > ceiling
  if (implausible) {
    return(skip_row_matched(cutoff, n_treated, "implausible_fit"))
  }

  wt <- pred$tr$n_post / sum(pred$tr$n_post)
  # savings = untreated - treated
  est <- sum(wt * (pred$mu0 - pred$mu1))
  # delta-method gradient
  g <- colSums(wt * (pred$mu0 * pred$X0 - pred$mu1 * pred$X1))
  V <- sandwich::vcovCL(pred$fit, cluster = ~group)
  se <- sqrt(drop(t(g) %*% V %*% g))
  if (!is.finite(se)) {
    return(skip_row_matched(cutoff, n_treated, "implausible_fit"))
  }
  result_row(
    "P+G",
    cutoff,
    est,
    se,
    nrow(pred$tr),
    pmpm_divisor(sim)
  )
}

#' Run every estimator (and the truth) on one simulated case
#'
#' `true_effect()`, pre/post, ITS, CITS, and DID only break down at
#' literally zero engaged members (`true_effect()`'s target is
#' `sum(ts[eng]) / sum(eng)` with an empty `eng` -- 0/0, not a
#' precision problem to tune a floor against), so the only guard
#' here is that mathematical necessity, not a chosen threshold.
#' `est_matched()` has its own, much higher floor (`min_treated`,
#' default 20) plus a post-fit plausibility check, because it needs
#' far more data to identify its outcome model and -- unlike the
#' other five -- can return a finite, unflagged, wildly wrong number
#' well past that floor (see its roxygen note). Below zero engaged,
#' nothing downstream is even definable, so the whole cutoff returns
#' `NULL`; between 1 and matching's floor, this returns 5 real rows
#' plus one recorded matching skip rather than discarding all six.
#' @export
estimate_case <- function(sim, cutoff = 11L, min_treated_matched = 20L) {
  n_treated <- sum(
    !is.na(sim$members$engaged_month) &
      sim$members$engaged_month <= cutoff
  )
  if (n_treated < 1L) return(NULL)
  rbind(
    true_effect(sim, cutoff),
    est_prepost(sim, cutoff),
    est_its(sim, cutoff),
    est_cits(sim, cutoff),
    est_did(sim, cutoff),
    est_matched(sim, cutoff, min_treated = min_treated_matched)
  )
}
