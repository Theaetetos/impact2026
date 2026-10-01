#' Parameters for the simulated self-funded (ASO) case
#'
#' Defaults are calibrated to the Part 1 commercial crosswalk
#' (slides 14-15): relative effects of 15% fewer ED visits, 6.5%
#' fewer admissions, 20% fewer medication-related visits, and 1.0
#' added primary care visit a year, with the slide 15 unit costs. The
#' block is the 40 of 140 self-funded groups that bought the program
#' (about 17,000 of about 60,000 members).
#'
#' Rates are monthly. Month indices run from `-n_pre` to
#' `n_post - 1`, with the program launching at month 0.
#'
#' @param ... Named values that override the defaults.
#' @return A named list of parameters.
#' @export
case_params <- function(...) {
  p <- list(
    # Block structure
    n_groups        = 140L,
    n_buying        = 40L,
    mean_group_size = 430,
      # sd of log group cost level
    group_sd        = 0.15,
      # logit slope of buying on group cost level
    buy_selection   = 1.0,
    n_pre           = 24L,
    n_post          = 12L,

      # Event drivers: base monthly rate for an average member
    drivers   = c("ed", "adm", "med", "pcp"),
    base_rate = c(ed = 0.16, adm = 0.045, med = 0.03, pcp = 2.4) / 12,
    unit_cost = c(ed = 2800, adm = 30000, med = 3000, pcp = 150),

      # Persistent frailty (lognormal) and how strongly each driver
      # loads on it
    frailty_sd      = 1.05,
    frailty_loading = c(ed = 1.0, adm = 1.0, med = 1.0, pcp = 0.3),

      # Transient health shocks: AR(1) on the log rate, monthly.
      # Eligibility selects on these spikes -> regression to the mean.
    ar_rho        = 0.85,
    ar_sd         = 0.30,  # innovation sd
    shock_loading = c(ed = 1.0, adm = 0.9, med = 1.0, pcp = 0.2),
      # winter ED peak
    ed_season_amp = 0.08,

      # Something else that changed at launch, block-wide (e.g. a
      # virtual care benefit). Hits buyers and non-buyers alike;
      # biases ITS and pre/post.
    concurrent_ed_change = -0.04,

      # Eligibility trigger over trailing 12 months (commercial
      # version of Part 1)
    trigger_ed     = 3L,
      # 1 admission plus this many ED visits
    trigger_adm_ed = 2L,

      # Engagement: CHW capacity, prioritised by predicted cost (risk
      # score)
        # ~1.5 CHWs on 60k -> ~0.43 CHW on 17k, 20/month each
    slots_per_month  = 8.5,
    accept_intercept = 0.35,
        # logit slope on latent willingness
    accept_willing   = 1.0,
        # sd of noise between log frailty and risk score
    risk_score_noise = 0.5,

      # Latent willingness lowers post-trigger use for everyone who
      # has it, engaged or not. Unobserved, so matching cannot
      # remove it.
    willing_effect = 0.06,

      # Program effect (relative for ed/adm/med, additive
      # visits/month for pcp)
    effect      = c(ed = -0.15, adm = -0.065, med = -0.20),
    pcp_add     = 1.0 / 12,
    effect_ramp = c(1 / 3, 2 / 3, 1)
  )
  over <- list(...)
  bad <- setdiff(names(over), names(p))
  if(length(bad)){
    stop("Unknown parameter(s): ", paste(bad, collapse = ", "))
  }
  p[names(over)] <- over
  return(p)
}

#' Simulate one realisation of the self-funded case
#'
#' Simulates member-month event counts for every member of all 140
#' groups, eligibility on a trailing-12-month trigger,
#' capacity-limited engagement in the buying groups, and the program
#' effect. Because the counterfactual (untreated) rate is known for
#' every engaged member-month, the true savings are returned
#' alongside the data an analyst would see.
#'
#' @param params Output of [case_params()].
#' @param seed Optional RNG seed.
#' @return A list with `members` (one row per member; columns
#'   starting with `.` are latent and must not be used by
#'   estimators), `cost` (member x month matrix of allowed dollars on
#'   the four drivers), `counts` (list of member x month count
#'   matrices), `true_savings` (member x month matrix of expected
#'   dollars saved), `months`, and `params`.
#' @export
simulate_case <- function(params = case_params(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  p <- params
  months <- seq(-p$n_pre, p$n_post - 1L)
  n_t <- length(months)
  launch_col <- p$n_pre + 1L

  # --- Groups and members -------------------------------------------
  g_size <- pmax(50L, stats::rpois(p$n_groups, p$mean_group_size))
  g_level <- stats::rnorm(p$n_groups, 0, p$group_sd)
  # Higher-cost groups are more likely to buy: rank on a noisy logit
  # score
  buy_score <- (p$buy_selection * g_level / p$group_sd
                + stats::rlogis(p$n_groups))
  g_buy <- rank(-buy_score) <= p$n_buying

  group <- rep(seq_len(p$n_groups), g_size)
  n <- length(group)
  age <- sample(18:64, n, replace = TRUE)
  female <- stats::rbinom(n, 1, 0.5)
  z <- 0.5 * scale(age)[, 1] + sqrt(1 - 0.25) * stats::rnorm(n)
  log_frailty <- p$frailty_sd * z - p$frailty_sd^2 / 2
  conditions <- stats::rpois(n, exp(-0.6 + 0.7 * z))
  risk_score <- log_frailty + stats::rnorm(n, 0, p$risk_score_noise)
  willing <- stats::rnorm(n)

  # --- Storage ------------------------------------------------------
  counts <- lapply(
    stats::setNames(p$drivers, p$drivers),
    function(d) matrix(0L, n, n_t)
  )
  true_savings <- matrix(0, n, n_t)

  eligible_col <- rep(NA_integer_, n)
  engaged_col  <- rep(NA_integer_, n)
  offered_col  <- rep(NA_integer_, n)
  buying <- g_buy[group]

  shock <- stats::rnorm(n, 0, p$ar_sd / sqrt(1 - p$ar_rho^2))
  base_log <- outer(
    log_frailty,
    p$frailty_loading[p$drivers]
  ) + g_level[group]
  cumul_slot <- 0

  for (j in seq_len(n_t)) {
    m <- months[j]
    shock <- p$ar_rho * shock + stats::rnorm(n, 0, p$ar_sd)

    # Engagement offers at the start of each post month
    if (m >= 0) {
      slots <- (floor(cumul_slot + p$slots_per_month)
                - floor(cumul_slot))
      cumul_slot <- cumul_slot + p$slots_per_month
      pool <- which(
        buying &
          !is.na(eligible_col) &
          eligible_col <= j &
          is.na(offered_col)
      )
      pool <- pool[order(-risk_score[pool])]
      filled <- 0L
      for (i in pool) {
        if (filled >= slots) break
        offered_col[i] <- j
        tst <- stats::plogis(
          p$accept_intercept + p$accept_willing * willing[i]
        )
        if (stats::runif(1) < tst) {
          engaged_col[i] <- j
          filled <- filled + 1L
        }
      }
    }

    season <- 1 + p$ed_season_amp * cos(2 * pi * (m %% 12) / 12)
    post_trigger <- !is.na(eligible_col) & eligible_col <= j
    will_mult <- ifelse(
      post_trigger,
      exp(-p$willing_effect * willing),
      1
    )
    k <- ifelse(
      is.na(engaged_col) | engaged_col > j,
      0L,
      pmin(j - engaged_col + 1L, length(p$effect_ramp))
    )
    ramp <- c(0, p$effect_ramp)[k + 1L]

    for (d in p$drivers) {
      lam0 <- p$base_rate[[d]] * exp(
        base_log[, d] + p$shock_loading[[d]] * shock
      ) * will_mult
      if (d == "ed") {
        lam0 <- lam0 * season * ifelse(
          m >= 0,
          1 + p$concurrent_ed_change,
          1
        )
      }
      # if/else, not ifelse(): the test is a scalar, and ifelse()
      # returns a result the length of its test -- i.e. member 1's
      # rate recycled to all n members.
      lam1 <- if (d == "pcp") {
        lam0 + ramp * p$pcp_add
      } else {
        lam0 * (1 + ramp * p$effect[[d]])
      }
      counts[[d]][, j] <- stats::rpois(n, lam1)
      add <- (lam0 - lam1) * p$unit_cost[[d]]
      true_savings[, j] <- true_savings[, j] + add
    }

    # Eligibility, evaluated at month end from month -1 on (prevalent
    # pool at launch)
    if (m >= -1) {
      win <- max(1L, j - 11L):j
      ed12 <- rowSums(counts$ed[, win, drop = FALSE])
      adm12 <- rowSums(counts$adm[, win, drop = FALSE])
      hit <- (ed12 >= p$trigger_ed) |
        (adm12 >= 1 & ed12 >= p$trigger_adm_ed)
      new <- hit & is.na(eligible_col)
      eligible_col[new] <- max(j + 1L, launch_col)
    }
  }

  cost <- Reduce(
    `+`,
    Map(function(d) counts[[d]] * p$unit_cost[[d]], p$drivers)
  )

  col_to_month <- function(x) ifelse(is.na(x), NA_integer_, months[x])
  members <- data.frame(
    id = seq_len(n),
    group = group,
    buying = buying,
    age = age,
    female = female,
    conditions = conditions,
    risk_score = risk_score,
    eligible_month = col_to_month(eligible_col),
    offered_month  = col_to_month(offered_col),
    engaged_month  = col_to_month(engaged_col),
    .log_frailty = log_frailty,
    .willing = willing
  )
  list(members = members, cost = cost, counts = counts,
       true_savings = true_savings, months = months, params = p)
}
