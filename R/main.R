#' wired: Weighted Adaptive Prediction with Structured Dependence
#'
#' @param ts_set A matrix, or data frame of numeric time series.
#' @param future Integer scalar: forecast horizon used both for
#'   marginal models and for the dependence transform lag.
#' @param dates Vector of date values for the plot. Default: NULL.
#' @param mode Transformation to be applied to the time series: one of `"additive"`, `"multiplicative"`, `"log_multiplicative"`.
#' @param n_testing Integer; number of expanding-window evaluation points. Default: 30.
#' @param dep_metric Dependence estimator for the correlation prototype:
#'   `"kendall"`, `"spearman"` (rank-based; mapped to Gaussian/t correlation),
#'   or `"pearson"` (linear correlation).
#' @param corr_adapt Time-adaptation mode for correlation:
#'   - `"static"`: single correlation from all aligned history,
#'   - `"ewma"`: exponentially weighted correlation (fast-reacting),
#'   - `"rolling"`: correlation from the last `roll_window` rows,
#'   - `"regime"`: blend calm vs stress correlations using a stress score.
#' @param ewma_lambda Numeric in (0,1); higher values react faster in `"ewma"`. Effective memory is about 1/lambda.
#' @param roll_window Integer; rolling window size for `"rolling"` and as a fallback in `"regime"`. It is truncated to available rows if necessary.
#' @param shrink_alpha Numeric in (0,1); shrink correlation toward identity to stabilize inversion and PD repair.
#' @param copula Copula family: `"gaussian"` or `"t"`. The t-copula introduces symmetric tail dependence controlled by `t_df`.
#' @param t_df Degrees of freedom for the t-copula; must be > 2. Lower values increase tail dependence.
#' @param stress_fun Stress score used by `"regime"`:
#'   `"mean_abs"` = mean absolute transformed return per row;
#'   `"rms"` = root-mean-square per row.
#' @param calm_q,stress_q Numeric quantiles in (0,1) with `calm_q < stress_q`.
#'   Rows with stress lower than `calm_q` form the calm set; rows with stress greater than
#'   `stress_q` form the stress set. If either set is too small the method falls
#'   back to a rolling correlation.
#' @param stress_smooth Integer (greater than 1); length of a trailing moving average
#'   applied to the stress score to reduce noise.
#' @param stress_blend_k Positive scalar controlling logistic sharpness when
#'   blending calm/stress correlations at the latest stress value.
#'   Larger `k`, sharper switching.
#' @param seed Integer RNG seed used both for copula draws and mixture components.
#'   For strict reproducibility across runs/platforms, keep packages and R versions fixed.
#' @param u_eps Small positive number used to clip uniform copula draws away from
#'   0 and 1 to avoid quantile extremes or infinite transforms.
#' @param ... Additional arguments forwarded to internal functions.
#'
#'
#'
#' @return
#' A list with:
#' \describe{
#'   \item{res_by_h}{Named list \code{h1..hH} (one per horizon) of per-horizon fits and helpers.}
#'   \item{rfun_*}{Joint draw helpers: \code{rfun_trafo(n)} and \code{rfun_level(n)} return 3-D arrays
#'     \eqn{H \times n \times p} (transformed vs level scale), and \code{rfun_both(n)} returns
#'     \code{list(trafo=..., level=...)} with the same shapes.}
#'   \item{plot}{Recorded base R plot object.}
#'   \item{meta}{Wrapper-level settings and controls (e.g., \code{future}, \code{mode}, \code{n_testing},
#'     dependence/correlation and copula parameters, and regime-stress controls).}
#' }
#'
#'
#' @importFrom stats simulate cor filter rnorm dnorm qnorm pnorm rchisq pt quantile median sd approxfun approx ecdf density lm predict residuals mad runif
#' @importFrom utils tail head
#' @importFrom imputeTS na_kalman
#' @import graphics
#' @import grDevices
#' @importFrom MASS mvrnorm
#' @importFrom mc2d rpert dpert ppert qpert
#' @importFrom forecast auto.arima forecast
#' @importFrom quantreg rq
#'
#'
#'
#' @examples
#'\donttest{
#'  set.seed(1)
#'  n <- 300
#'  ts_set <- data.frame(
#'    A = 100 + cumsum(rnorm(n, 0, 1)),
#'    B =  80 + cumsum(rnorm(n, 0, 1))
#'  )
#'
#'  fitH <- wired(
#'    ts_set   = ts_set,
#'    future   = 2,
#'    mode     = "additive",
#'    n_testing = 2,
#'
#'    dep_metric = "spearman",
#'    corr_adapt = "rolling",
#'    roll_window = 40,
#'    copula = "gaussian",
#'    seed = 123,
#'
#'    n_crps_mc = 30,
#'    q_grid_size = 10
#'  )
#'
#'  draws_level <- fitH$rfun_level(5)
#'  print(dim(draws_level))
#'
#'  both <- fitH$rfun_both(5)
#'}
#'
#'
#' @export
wired <- function(
    ts_set,
    future,
    dates = NULL,
    mode = c("additive", "multiplicative", "log_multiplicative"),

    n_testing = 30,

    # pass-through dependence / copula controls
    dep_metric = c("kendall", "spearman", "pearson"),
    corr_adapt = c("static", "ewma", "rolling", "regime"),
    ewma_lambda = 0.15,
    roll_window = 60,
    shrink_alpha = 0.05,

    copula = c("gaussian", "t"),
    t_df = 7,

    stress_fun = c("mean_abs", "rms"),
    calm_q = 0.50,
    stress_q = 0.85,
    stress_smooth = 5,
    stress_blend_k = 8,

    seed = 123,
    u_eps = 1e-6,

    ...
) {

  ###check_feasibility(ts_set, future, n_testing)

  feas <- check_feasibility(ts_set, future, n_testing)

  if (!isTRUE(feas$feasible)) {
    stop(feas$message)
  }

  future <- feas$suggested_future
  n_testing <- feas$suggested_n_testing

  mode <- match.arg(mode, c("additive", "multiplicative", "log_multiplicative"))

  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  H <- as.integer(future)

  dep_metric <- match.arg(dep_metric)
  corr_adapt <- match.arg(corr_adapt)
  copula <- match.arg(copula)
  stress_fun <- match.arg(stress_fun)

  # Fit per-horizon
  res_by_h <- vector("list", H)
  for (h in seq_len(H)) {
    res_by_h[[h]] <- support_wired(
      ts_set = ts_set,
      future = h,
      mode = mode,

      n_testing = n_testing,

      dep_metric = dep_metric,
      corr_adapt = corr_adapt,
      ewma_lambda = ewma_lambda,
      roll_window = roll_window,
      shrink_alpha = shrink_alpha,

      copula = copula,
      t_df = t_df,

      stress_fun = stress_fun,
      calm_q = calm_q,
      stress_q = stress_q,
      stress_smooth = stress_smooth,
      stress_blend_k = stress_blend_k,

      seed = seed + 1000L * h,  # deterministic, but different per horizon
      u_eps = u_eps,

      ...
    )
  }
  names(res_by_h) <- paste0("h", seq_len(H))

  # infer p and names from h=1 (from a small draw to get colnames robustly)
  tmp <- res_by_h[[1]]$rfun_trafo(2)
  p <- ncol(tmp)
  nm <- colnames(tmp)
  if (is.null(nm)) nm <- paste0("S", seq_len(p))

  # helpers to allocate arrays
  alloc_arr <- function(n) {
    out <- array(NA_real_, dim = c(H, n, p))
    dimnames(out) <- list(paste0("h", seq_len(H)), NULL, nm)
    out
  }

  rfun_trafo <- function(n) {
    n <- as.integer(n)
    if (n <= 0) return(array(numeric(0), dim = c(0, 0, 0)))
    out <- alloc_arr(n)
    for (h in seq_len(H)) out[h, , ] <- res_by_h[[h]]$rfun_trafo(n)
    out
  }

  rfun_level <- function(n) {
    n <- as.integer(n)
    if (n <= 0) return(array(numeric(0), dim = c(0, 0, 0)))
    out <- alloc_arr(n)
    for (h in seq_len(H)) out[h, , ] <- res_by_h[[h]]$rfun_level(n)
    out
  }

  rfun_both <- function(n) {
    n <- as.integer(n)
    if (n <= 0) return(list(
      trafo = array(numeric(0), dim = c(0, 0, 0)),
      level = array(numeric(0), dim = c(0, 0, 0))
    ))

    out_t <- alloc_arr(n)
    out_l <- alloc_arr(n)

    for (h in seq_len(H)) {
      both_h <- res_by_h[[h]]$rfun_both(n)
      out_t[h, , ] <- both_h$trafo
      out_l[h, , ] <- both_h$level
    }

    list(trafo = out_t, level = out_l)
  }

  p_graph <- plot_wired(res_by_h, ts_set, probs = c(0.1, 0.25, 0.50, 0.75, 0.9), x = dates, history_n = nrow(ts_set))

  list(
    res_by_h = res_by_h,
    rfun_trafo = rfun_trafo,
    rfun_level = rfun_level,
    rfun_both = rfun_both,
    plot = p_graph,
    meta = list(
      future = H,
      mode = mode,
      n_testing = n_testing,
      dep_metric = dep_metric,
      corr_adapt = corr_adapt,
      ewma_lambda = ewma_lambda,
      roll_window = roll_window,
      shrink_alpha = shrink_alpha,
      copula = copula,
      t_df = t_df,
      stress_fun = stress_fun,
      calm_q = calm_q,
      stress_q = stress_q,
      stress_smooth = stress_smooth,
      stress_blend_k = stress_blend_k,
      u_eps = u_eps
    )
  )
}



#' @keywords internal
support_wired <- function(
    ts_set,
    future,
    mode = c("additive", "multiplicative", "log_multiplicative"),

    n_testing = 50,

    dep_metric = c("kendall", "spearman", "pearson"),
    corr_adapt = c("static", "ewma", "rolling", "regime"),
    ewma_lambda = 0.15,
    roll_window = 60,
    shrink_alpha = 0.05,

    copula = c("gaussian", "t"),
    t_df = 7,

    stress_fun = c("mean_abs", "rms"),
    calm_q = 0.50,
    stress_q = 0.85,
    stress_smooth = 5,
    stress_blend_k = 8,

    seed = 123,
    u_eps = 1e-6,

    ...
) {
  mode <- match.arg(mode, c("additive", "multiplicative", "log_multiplicative"))

  dep_metric <- match.arg(dep_metric)
  corr_adapt <- match.arg(corr_adapt)
  copula <- match.arg(copula)
  stress_fun <- match.arg(stress_fun)

  if (!requireNamespace("MASS", quietly = TRUE)) stop("Package 'MASS' is required (mvrnorm).")

  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)

  # ----- normalize input -----
  if (is.matrix(ts_set) || is.data.frame(ts_set)) {
    ts_set <- as.data.frame(ts_set)
    ts_set <- lapply(ts_set, as.numeric)
  } else if (is.list(ts_set)) {
    ts_set <- lapply(ts_set, as.numeric)
  } else stop("ts_set must be a list or a matrix/data.frame of series.")
  p <- length(ts_set)
  if (p < 2) stop("Provide at least 2 time series.")
  if (is.null(names(ts_set))) names(ts_set) <- paste0("S", seq_len(p))

  # ----- helpers -----
  clip_u <- function(u) pmin(pmax(u, u_eps), 1 - u_eps)

  make_corr_pd <- function(R) {
    R <- as.matrix(R)
    R <- (R + t(R)) / 2
    diag(R) <- 1
    eig <- eigen(R, symmetric = TRUE)
    vals <- eig$values
    vecs <- eig$vectors
    floor_val <- max(1e-10, 1e-8 * max(vals, na.rm = TRUE))
    vals2 <- pmax(vals, floor_val)
    R2 <- vecs %*% diag(vals2, nrow = length(vals2)) %*% t(vecs)
    R2 <- (R2 + t(R2)) / 2
    diag(R2) <- 1
    R2[R2 > 1] <- 1
    R2[R2 < -1] <- -1
    R2
  }

  shrink_corr <- function(R, alpha) {
    if (!is.finite(alpha) || alpha <= 0) return(R)
    alpha <- min(max(alpha, 0), 0.99)
    R2 <- (1 - alpha) * R + alpha * diag(1, nrow(R))
    diag(R2) <- 1
    R2
  }

  rank_to_rho <- function(M, metric) {
    if (metric == "pearson") return(M)
    if (metric == "kendall") return(sin(pi * M / 2))
    2 * sin(pi * M / 6) # spearman
  }

  ewma_corr <- function(Z, lambda) {
    lambda <- max(min(lambda, 1), 1e-6)
    S <- matrix(0, ncol(Z), ncol(Z))
    for (i in seq_len(nrow(Z))) {
      x <- as.numeric(Z[i, ])
      x[!is.finite(x)] <- 0
      S <- (1 - lambda) * S + lambda * (x %*% t(x))
    }
    d <- sqrt(pmax(diag(S), 1e-12))
    R <- S / (d %o% d)
    diag(R) <- 1
    R
  }

  rolling_corr <- function(Z, window, metric) {
    w <- min(as.integer(window), nrow(Z))
    Z2 <- tail(Z, w)
    M <- if (metric == "pearson") cor(Z2, use = "pairwise.complete.obs") else
      cor(Z2, method = metric, use = "pairwise.complete.obs")
    R <- rank_to_rho(M, metric)
    diag(R) <- 1
    R
  }

  stress_series <- function(Z, smooth_k = 5) {
    raw <- if (stress_fun == "mean_abs") rowMeans(abs(Z), na.rm = TRUE) else sqrt(rowMeans(Z^2, na.rm = TRUE))
    k <- max(1L, as.integer(smooth_k))
    if (k == 1L) return(raw)
    filt <- rep(1 / k, k)
    as.numeric(filter(raw, filt, sides = 1, method = "convolution"))
  }

  logistic_blend <- function(x, x0, k) 1 / (1 + exp(-k * (x - x0)))

  deprecated_draw_U <- function(n, R) {
    n <- as.integer(n)
    if (copula == "gaussian") {
      Z <- mvrnorm(n = n, mu = rep(0, p), Sigma = R)
      return(clip_u(pnorm(Z)))
    }
    df <- as.numeric(t_df)
    if (!is.finite(df) || df <= 2) stop("t_df must be > 2.")
    Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = R)
    g <- rchisq(n, df = df) / df
    Zt <- Y / sqrt(g)
    clip_u(pt(Zt, df = df))
  }

  draw_U <- function(n, R) {
    n <- as.integer(n)
    if (n <= 0) return(matrix(numeric(0), nrow = 0, ncol = p))

    if (copula == "gaussian") {
      Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = R)
      if (n == 1) Z <- matrix(Z, nrow = 1)  # <--- critical
      U <- stats::pnorm(Z)
      return(clip_u(U))
    }

    df <- as.numeric(t_df)
    if (!is.finite(df) || df <= 2) stop("t_df must be > 2.")
    Y <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = R)
    if (n == 1) Y <- matrix(Y, nrow = 1)    # <--- critical
    g <- stats::rchisq(n, df = df) / df
    Zt <- Y / sqrt(g)                       # vector recycling works row-wise when g length n
    U <- stats::pt(Zt, df = df)
    clip_u(U)
  }

  # ============================================================
  # 1) Fit TWO sets of marginals: transformed + (optional) levels
  # ============================================================
  marg_trafo <- vector("list", p)
  names(marg_trafo) <- names(ts_set)

  marg_level <- NULL
  if (!is.null(mode)) {
    marg_level <- vector("list", p)
    names(marg_level) <- names(ts_set)
  }

  for (j in seq_len(p)) {
    pf_t <- probabilistic_mixer(
      ts = ts_set[[j]],
      future = future,
      n_testing = n_testing,
      level_mode = NULL,
      ...
    )
    if (!is.list(pf_t) || is.null(pf_t$qfun) || is.null(pf_t$rfun) || is.null(pf_t$pfun)) {
      stop(sprintf("Invalid transformed pred_fun for '%s'.", names(ts_set)[j]))
    }
    marg_trafo[[j]] <- pf_t

    if (!is.null(mode)) {
      pf_l <- probabilistic_mixer(
        ts = ts_set[[j]],
        future = future,
        n_testing = n_testing,
        level_mode = mode,
        ...
      )
      if (!is.list(pf_l) || is.null(pf_l$qfun) || is.null(pf_l$rfun) || is.null(pf_l$pfun)) {
        stop(sprintf("Invalid level pred_fun for '%s'.", names(ts_set)[j]))
      }
      marg_level[[j]] <- pf_l
    }
  }

  # ============================================================
  # 2) Dependence estimation from transformed history
  # ============================================================

  Z_list <- vector("list", p)
  for (j in seq_len(p)) {
    z <- trafo(ts_set[[j]], lag = future, mode = mode)
    z <- as.numeric(z)
    z <- z[is.finite(z)]
    if (length(z) < 8) stop(sprintf("Transformed history too short for '%s'.", names(ts_set)[j]))
    Z_list[[j]] <- z
  }

  L <- min(vapply(Z_list, length, integer(1)))
  if (L < 8) stop("Not enough aligned transformed history to estimate dependence.")
  Z_mat <- do.call(cbind, lapply(Z_list, function(z) tail(z, L)))
  colnames(Z_mat) <- names(ts_set)

  R0 <- switch(
    corr_adapt,
    static = {
      M <- if (dep_metric == "pearson") cor(Z_mat, use = "pairwise.complete.obs") else
        cor(Z_mat, method = dep_metric, use = "pairwise.complete.obs")
      R <- rank_to_rho(M, dep_metric); diag(R) <- 1; R
    },
    ewma = {
      R <- ewma_corr(Z_mat, lambda = ewma_lambda); diag(R) <- 1; R
    },
    rolling = {
      R <- rolling_corr(Z_mat, window = roll_window, metric = dep_metric); diag(R) <- 1; R
    },
    regime = {
      s <- stress_series(Z_mat, smooth_k = stress_smooth)
      ok <- is.finite(s)
      Z_ok <- Z_mat[ok, , drop = FALSE]
      s_ok <- s[ok]
      if (nrow(Z_ok) < 10) stop("Not enough rows after stress smoothing for regime correlation.")

      q_calm <- quantile(s_ok, probs = calm_q, names = FALSE, type = 8)
      q_stress <- quantile(s_ok, probs = stress_q, names = FALSE, type = 8)

      calm_rows <- which(s_ok <= q_calm)
      stress_rows <- which(s_ok >= q_stress)

      if (length(calm_rows) < 6 || length(stress_rows) < 6) {
        R <- rolling_corr(Z_ok, window = min(roll_window, nrow(Z_ok)), metric = dep_metric)
      } else {
        Mc <- if (dep_metric == "pearson") cor(Z_ok[calm_rows, , drop = FALSE], use = "pairwise.complete.obs") else
          cor(Z_ok[calm_rows, , drop = FALSE], method = dep_metric, use = "pairwise.complete.obs")
        Ms <- if (dep_metric == "pearson") cor(Z_ok[stress_rows, , drop = FALSE], use = "pairwise.complete.obs") else
          cor(Z_ok[stress_rows, , drop = FALSE], method = dep_metric, use = "pairwise.complete.obs")

        Rc <- rank_to_rho(Mc, dep_metric); diag(Rc) <- 1
        Rs <- rank_to_rho(Ms, dep_metric); diag(Rs) <- 1

        s_last <- tail(s_ok, 1)
        s_mid <- (q_calm + q_stress) / 2
        w <- logistic_blend(s_last, x0 = s_mid, k = stress_blend_k)

        R <- (1 - w) * Rc + w * Rs
        diag(R) <- 1
      }
      R
    }
  )

  R <- shrink_corr(R0, shrink_alpha)
  R <- make_corr_pd(R)
  colnames(R) <- rownames(R) <- names(ts_set)

  # ============================================================
  # 3) Joint sampling in BOTH spaces (same U per scenario)
  # ============================================================
  deprecated_map_q <- function(U, marg_list) {
    n <- nrow(U)
    X <- matrix(NA_real_, nrow = n, ncol = p)
    colnames(X) <- names(ts_set)
    for (j in seq_len(p)) X[, j] <- as.numeric(marg_list[[j]]$qfun(U[, j]))
    X
  }

  map_q <- function(U, marg_list) {
    # normalize U to n x p matrix
    if (is.null(dim(U))) {
      # vector case, assume length p and n=1
      if (length(U) != p) stop("U vector length != p.")
      U <- matrix(U, nrow = 1)
    } else {
      # matrix case: ensure correct ncol
      if (ncol(U) != p) stop("U matrix must have p columns.")
    }

    n <- nrow(U)
    X <- matrix(NA_real_, nrow = n, ncol = p)
    colnames(X) <- names(ts_set)

    for (j in seq_len(p)) {
      X[, j] <- as.numeric(marg_list[[j]]$qfun(U[, j]))
    }
    X
  }

  rfun_trafo <- function(n) {
    set.seed(seed)
    U <- draw_U(n, R)
    map_q(U, marg_trafo)
  }

  rfun_level <- if (is.null(mode)) {
    rfun_trafo
  } else {
    function(n) {
      set.seed(seed)
      U <- draw_U(n, R)
      map_q(U, marg_level)
    }
  }

  rfun_both <- function(n) {
    set.seed(seed)
    U <- draw_U(n, R)
    out_trafo <- map_q(U, marg_trafo)
    out_level <- if (is.null(mode)) out_trafo else map_q(U, marg_level)
    list(trafo = out_trafo, level = out_level)
  }

  list(
    marginals_trafo = marg_trafo,
    marginals_level = marg_level,   # NULL if mode=NULL
    R = R,
    rfun_trafo = rfun_trafo,
    rfun_level = rfun_level,
    rfun_both = rfun_both,
    meta = list(
      future = future,
      mode = mode,
      dep_metric = dep_metric,
      corr_adapt = corr_adapt,
      ewma_lambda = ewma_lambda,
      roll_window = roll_window,
      shrink_alpha = shrink_alpha,
      copula = copula,
      t_df = t_df,
      stress_fun = stress_fun,
      calm_q = calm_q,
      stress_q = stress_q,
      stress_smooth = stress_smooth,
      stress_blend_k = stress_blend_k,
      aligned_history_length = L,
      u_eps = u_eps
    )
  )
}



#' @keywords internal

probabilistic_mixer <- function(
    ts,
    future,
    n_testing = 50,
    n_crps_mc = 2000,
    seed = 123,
    temperature = NULL,
    q_eps = 1e-4,
    q_grid_size = 4000,
    level_mode = NULL
) {
  stopifnot(is.numeric(ts), length(ts) >= 20)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)

  # ---------- helpers ----------
  crps_mc <- function(x, y) {
    x <- as.numeric(x)
    x <- x[is.finite(x)]
    m <- length(x)
    if (m < 5 || !is.finite(y)) return(NA_real_)
    x <- sort(x)
    term1 <- mean(abs(x - y))
    i <- seq_len(m)
    term2 <- sum((2 * i - m - 1) * x) / (m * m)
    term1 - term2
  }

  draw_pred <- function(pf, n, shift = 0L) {
    set.seed(seed + shift)
    pf$rfun(n)
  }

  softmax <- function(z) {
    z <- z - max(z)
    w <- exp(z)
    w / sum(w)
  }

  predict_next_crps <- function(crps, t_idx, K = 10) {
    stopifnot(length(crps) == length(t_idx))

    ok <- is.finite(crps) & is.finite(t_idx)
    x <- crps[ok]
    t <- t_idx[ok]

    n <- length(x)
    if (n == 0) return(Inf)
    if (n == 1) return(x[1])
    if (n == 2) return(x[2])  # last observed

    # use only the last K points
    k <- min(K, n)
    xk <- tail(x, k)
    tk <- tail(t, k)

    # ----- robust slope: Theil–Sen -----
    slopes <- numeric(0)
    for (i in 1:(k - 1)) {
      dt <- tk[(i + 1):k] - tk[i]
      valid <- dt != 0
      if (any(valid)) {
        slopes <- c(slopes, (xk[(i + 1):k] - xk[i])[valid] / dt[valid])
      }
    }

    beta <- if (length(slopes) > 0) median(slopes) else 0

    # ----- robust intercept via centering -----
    t0 <- median(tk)
    x0 <- median(xk)

    t_next <- new_idx
    pred <- x0 + beta * (t_next - t0)

    if (!is.finite(pred)) pred <- tail(xk, 1)
    as.numeric(pred)
  }

  # ---------- build testing windows ----------
  n <- length(ts)
  new_idx <- n

  max_train_end <- n - future
  test_idx <- round(seq.int(from = floor(n/3), to = max_train_end, length.out = n_testing))

  crps_nv <- crps_ar <-crps_eg <- crps_hb <- crps_dr <- crps_vs <- crps_rm <- crps_sq <- numeric(length(test_idx))

  # ---------- backtest ----------
  for (k in seq_along(test_idx)) {
    end <- test_idx[k]
    ts_tr <- ts[1:end]

    d_full <- trafo(ts[1:(end + future)], future, mode = level_mode)
    y_true <- tail(d_full, 1)

    pf_nv <- naive_predictor(ts_tr, future, level_mode = level_mode)
    pf_ar <- arima_predictor(ts_tr, future, level_mode = level_mode)
    pf_eg <- ewma_gaussian_predictor(ts_tr, future, level_mode = level_mode)
    pf_hb <- historical_bootstrap_predictor(ts_tr, future, level_mode = level_mode)
    pf_dr <- drift_resid_bootstrap_predictor(ts_tr, future, level_mode = level_mode)
    pf_vs <- vol_scaled_naive_predictor(ts_tr, future, level_mode = level_mode)
    pf_rm <- robust_median_mad_predictor(ts_tr, future, level_mode = level_mode)
    pf_sq <- shrunk_quantile_predictor(ts_tr, future, level_mode = level_mode)

    crps_nv[k] <- crps_mc(draw_pred(pf_nv, n_crps_mc, sample.int(1000 + k, 1)), y_true)
    crps_ar[k] <- crps_mc(draw_pred(pf_ar, n_crps_mc, sample.int(1000 + k, 1)), y_true)
    crps_eg[k] <- crps_mc(draw_pred(pf_eg, n_crps_mc, sample.int(1000 + k, 1)), y_true)
    crps_hb[k] <- crps_mc(draw_pred(pf_hb, n_crps_mc, sample.int(1000 + k, 1)), y_true)
    crps_dr[k] <- crps_mc(draw_pred(pf_dr, n_crps_mc, sample.int(1000 + k, 1)), y_true)
    crps_vs[k] <- crps_mc(draw_pred(pf_vs, n_crps_mc, sample.int(1000 + k, 1)), y_true)
    crps_rm[k] <- crps_mc(draw_pred(pf_rm, n_crps_mc, sample.int(1000 + k, 1)), y_true)
    crps_sq[k] <- crps_mc(draw_pred(pf_sq, n_crps_mc, sample.int(1000 + k, 1)), y_true)
  }

  # ---------- predictive CRPS → weights ----------

  pred_crps <- c(
    nv    = predict_next_crps(crps_nv, test_idx),
    ar    = predict_next_crps(crps_ar, test_idx),
    eg    = predict_next_crps(crps_eg, test_idx),
    hb    = predict_next_crps(crps_hb, test_idx),
    dr	  = predict_next_crps(crps_dr, test_idx),
    vs    = predict_next_crps(crps_vs, test_idx),
    rm	  = predict_next_crps(crps_rm, test_idx),
    sq    = predict_next_crps(crps_sq, test_idx)
  )

  if (is.null(temperature))
    temperature <- max(sd(pred_crps), 1e-8)

  weights <- softmax(-pred_crps / temperature)

  # ---------- fit final predictors ----------

  pf_NV <- naive_predictor(ts, future, level_mode = level_mode)
  pf_AR <- arima_predictor(ts, future, level_mode = level_mode)
  pf_EG <- ewma_gaussian_predictor(ts, future, level_mode = level_mode)
  pf_HB <- historical_bootstrap_predictor(ts, future, level_mode = level_mode)
  pf_DR <- drift_resid_bootstrap_predictor(ts, future, level_mode = level_mode)
  pf_VS <- vol_scaled_naive_predictor(ts, future, level_mode = level_mode)
  pf_RM <- robust_median_mad_predictor(ts, future, level_mode = level_mode)
  pf_SQ <- shrunk_quantile_predictor(ts, future, level_mode = level_mode)

  pfs <- list(nv = pf_NV, ar = pf_AR, eg = pf_EG, hb = pf_HB, dr = pf_DR, vs = pf_VS, rm = pf_RM, sq = pf_SQ)
  w <- weights

  # ---------- mixture functions ----------
  dfun <- function(x)
    w["nv"] * pf_NV$dfun(x) +
    w["ar"] * pf_AR$dfun(x) +
    w["eg"] * pf_EG$dfun(x) +
    w["hb"] * pf_HB$dfun(x) +
    w["dr"] * pf_DR$dfun(x) +
    w["vs"] * pf_VS$dfun(x) +
    w["rm"] * pf_RM$dfun(x) +
    w["sq"] * pf_SQ$dfun(x)

  pfun <- function(x)
    w["nv"] * pf_NV$pfun(x) +
    w["ar"] * pf_AR$pfun(x) +
    w["eg"] * pf_EG$pfun(x) +
    w["hb"] * pf_HB$pfun(x) +
    w["dr"] * pf_DR$pfun(x) +
    w["vs"] * pf_VS$pfun(x) +
    w["rm"] * pf_RM$pfun(x) +
    w["sq"] * pf_SQ$pfun(x)

  rfun <- function(n) {
    comps <- sample(names(w), n, replace = TRUE, prob = w)
    out <- numeric(n)
    for (nm in names(w)) {
      idx <- which(comps == nm)
      if (length(idx)) out[idx] <- pfs[[nm]]$rfun(length(idx))
    }
    out
  }

  # ---------- quantile via approxfun ----------
  lo <- min(
    pf_NV$qfun(q_eps), pf_AR$qfun(q_eps), pf_EG$qfun(q_eps), pf_HB$qfun(q_eps), pf_DR$qfun(q_eps), pf_VS$qfun(q_eps), pf_RM$qfun(q_eps), pf_SQ$qfun(q_eps),
    na.rm = TRUE
  )
  hi <- max(
    pf_NV$qfun(1 - q_eps), pf_AR$qfun(1 - q_eps), pf_EG$qfun(1 - q_eps), pf_HB$qfun(1 - q_eps), pf_DR$qfun(1 - q_eps), pf_VS$qfun(1 - q_eps), pf_RM$qfun(1 - q_eps), pf_SQ$qfun(1 - q_eps),
    na.rm = TRUE
  )

  grid_x <- seq(lo, hi, length.out = q_grid_size)
  grid_p <- pmin(pmax(pfun(grid_x), 0), 1)
  grid_p <- cummax(grid_p)

  keep <- c(TRUE, diff(grid_p) > 0)
  inv_cdf <- approxfun(grid_p[keep], grid_x[keep], rule = 2)

  qfun <- function(p) inv_cdf(pmin(pmax(p, 0), 1))

  # ---------- output ----------
  list(
    rfun = rfun,
    dfun = dfun,
    pfun = pfun,
    qfun = qfun,
    weights = weights,
    pred_crps = pred_crps
  )
}



#' @keywords internal

.emp_pred <- function(samples, density_n = 512, seed = NULL) {
  s <- as.numeric(samples)
  s <- s[is.finite(s)]
  if (length(s) < 20) stop("Too few finite samples to build empirical distribution.")
  s_sort <- sort(s)
  ec <- ecdf(s_sort)
  dens <- density(s_sort, n = density_n)

  list(
    rfun = function(n) {
      if (!is.null(seed)) set.seed(seed)
      sample(s_sort, size = as.integer(n), replace = TRUE)
    },
    dfun = function(x) approx(dens$x, dens$y, xout = as.numeric(x), rule = 2)$y,
    pfun = function(q) ec(as.numeric(q)),
    qfun = function(p) {
      p <- pmin(pmax(as.numeric(p), 0), 1)
      as.numeric(quantile(s_sort, probs = p, type = 8, names = FALSE))
    }
  )
}


#' @keywords internal
# Wrap a predictor on deltas/growth into a predictor on LEVELS using last(ts)
.apply_level_mode <- function(pred_fun, last_level, level_mode = NULL) {
  if (is.null(level_mode)) return(pred_fun)

  if (!is.character(level_mode) || length(level_mode) != 1)
    stop("level_mode must be one of: NULL, 'additive', 'multiplicative', 'log_multiplicative'.")

  level_mode <- match.arg(level_mode, choices = c("additive", "multiplicative", "log_multiplicative"))

  a <- as.numeric(last_level)
  if (!is.finite(a)) stop("last level is not finite; cannot apply level_mode.")

  if (level_mode == "additive") {
    # Y = a + X
    return(list(
      rfun = function(n) a + pred_fun$rfun(n),
      qfun = function(p) a + pred_fun$qfun(p),
      pfun = function(q) pred_fun$pfun(as.numeric(q) - a),
      dfun = function(y) pred_fun$dfun(as.numeric(y) - a) # Jacobian = 1
    ))
  }

  if (level_mode == "multiplicative") {
    # Y = a * (1 + X)  => X = Y/a - 1, dx/dy = 1/a
    if (a <= 0) stop("multiplicative level_mode requires last(ts) > 0.")
    return(list(
      rfun = function(n) a * (1 + pred_fun$rfun(n)),
      qfun = function(p) a * (1 + pred_fun$qfun(p)),
      pfun = function(q) pred_fun$pfun(as.numeric(q) / a - 1),
      dfun = function(y) pred_fun$dfun(as.numeric(y) / a - 1) * (1 / a)
    ))
  }

  # log_multiplicative:
  # X = log(Y / a)  (i.e., log-return), so Y = a * exp(X)
  # Inverse: X = log(Y/a), Jacobian: dx/dy = 1/y
  if (a <= 0) stop("log_multiplicative level_mode requires last(ts) > 0.")
  list(
    rfun = function(n) a * exp(pred_fun$rfun(n)),
    qfun = function(p) a * exp(pred_fun$qfun(p)),
    pfun = function(q) {
      q <- as.numeric(q)
      # For q<=0, probability is 0 when a>0 and exp(.)>0
      out <- numeric(length(q))
      ok <- q > 0
      out[!ok] <- 0
      out[ok] <- pred_fun$pfun(log(q[ok] / a))
      out
    },
    dfun = function(y) {
      y <- as.numeric(y)
      out <- numeric(length(y))
      ok <- y > 0
      out[!ok] <- 0
      # f_Y(y) = f_X(log(y/a)) * (1/y)
      out[ok] <- pred_fun$dfun(log(y[ok] / a)) * (1 / y[ok])
      out
    }
  )
}

#' @keywords internal
trafo <- function(ts, lag = 1, mode = c("multiplicative", "log_multiplicative", "additive"), stride = TRUE) {
  mode <- match.arg(mode)
  ts <- as.numeric(ts)

  if (!is.numeric(lag) || length(lag) != 1 || !is.finite(lag) || lag < 1) {
    stop("lag must be a positive integer.")
  }
  lag <- as.integer(lag)

  n <- length(ts)
  if (n <= lag) stop("ts must have length > lag.")

  a <- head(ts, -lag)
  b <- tail(ts, -lag)

  # core transforms
  trafo_ts <- switch(
    mode,
    additive = b - a,
    multiplicative = (b / a) - 1,
    log_multiplicative = log(b / a)
  )

  # clean non-finite results (e.g., division by 0, log of <=0)
  trafo_ts[!is.finite(trafo_ts)] <- NA_real_

  # impute only if there are missing values in the TRANSFORMED series
  # (note: original code checked anyNA(ts), which can miss NA created by transform)
  if (anyNA(trafo_ts)) {
    if (!requireNamespace("imputeTS", quietly = TRUE)) {
      stop("Package 'imputeTS' is required for NA imputation (na_kalman).")
    }
    trafo_ts <- na_kalman(trafo_ts)
  }

  # stride: select every 'lag' step anchored to the END (includes latest)
  if (isTRUE(stride)) {
    idx <- seq.int(from = length(trafo_ts), to = 1L, by = -lag)
    trafo_ts <- trafo_ts[rev(idx)]
  }

  trafo_ts
}



#' @keywords internal

naive_predictor <- function(ts, future, seed = 123, level_mode = NULL) {
  stopifnot(is.numeric(ts), length(ts) >= 3)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)

  if (!requireNamespace("mc2d", quietly = TRUE)) stop("Package 'mc2d' is required.")

  y <- trafo(ts, future, mode = level_mode, stride = TRUE)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  if (length(y) < 3) stop("Downsampled series too short.")

  a <- min(y)
  b <- max(y)
  m <- as.numeric(tail(y, 1)) # mode = most recent
  if (a > b) { tmp <- a; a <- b; b <- tmp }
  m <- max(min(m, b), a)

  if (abs(b - a) < .Machine$double.eps) {
    base <- list(
      rfun = function(n) rep(a, as.integer(n)),
      dfun = function(x) ifelse(abs(as.numeric(x) - a) < .Machine$double.eps, Inf, 0),
      pfun = function(q) as.numeric(as.numeric(q) >= a),
      qfun = function(p) rep(a, length(p))
    )
    return(.apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode))
  }

  base <- list(
    rfun = function(n) { set.seed(seed); rpert(as.integer(n), min = a, mode = m, max = b) },
    dfun = function(x) dpert(as.numeric(x), min = a, mode = m, max = b),
    pfun = function(q) ppert(as.numeric(q), min = a, mode = m, max = b),
    qfun = function(p) qpert(pmin(pmax(as.numeric(p), 0), 1), min = a, mode = m, max = b)
  )

  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}


#' @keywords internal

arima_predictor <- function(ts, h = 5, level = 0.95, n_sim = 5000, seed = 123, level_mode = NULL) {
  stopifnot(is.numeric(ts), length(ts) >= 10)
  stopifnot(length(h) == 1, is.finite(h), h >= 1)
  h <- as.integer(h)

  if (!requireNamespace("forecast", quietly = TRUE)) stop("Package 'forecast' is required.")

  y <- trafo(ts, h, mode = level_mode, stride = TRUE)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  if (length(y) < 10) stop("Downsampled series too short for ARIMA.")

  fit <- auto.arima(y)

  set.seed(seed)
  sim_draws <- tryCatch({
    replicate(as.integer(n_sim), simulate(fit, nsim = 1, future = TRUE)[1])
  }, error = function(e) NULL)

  if (!is.null(sim_draws)) {
    sim_draws <- as.numeric(sim_draws)
    sim_draws <- sim_draws[is.finite(sim_draws)]
  }

  if (!is.null(sim_draws) && length(sim_draws) >= 50) {
    base <- .emp_pred(sim_draws, seed = seed)
    return(.apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode))
  }

  fc <- forecast(fit, h = 1, level = level * 100)

  mu <- as.numeric(fc$mean[1])
  upper <- as.numeric(fc$upper[1, 1])
  lower <- as.numeric(fc$lower[1, 1])

  z <- qnorm(0.5 + level / 2)
  sigma <- (upper - lower) / (2 * z)
  if (!is.finite(sigma) || sigma <= 0) sigma <- sqrt(fit$sigma2)
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1e-8

  base <- list(
    rfun = function(n) rnorm(as.integer(n), mean = mu, sd = sigma),
    dfun = function(x) dnorm(as.numeric(x), mean = mu, sd = sigma),
    pfun = function(q) pnorm(as.numeric(q), mean = mu, sd = sigma),
    qfun = function(p) qnorm(pmin(pmax(as.numeric(p), 0), 1), mean = mu, sd = sigma)
  )

  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}


#' @keywords internal

quantile_predictor <- function(
    ts, future,
    quants = seq(0.02, 0.98, length.out = 49),
    n_draws = 1000,
    seed = 123,
    rq_method = "fn",
    jitter_ties = TRUE,
    level_mode = NULL
) {
  stopifnot(is.numeric(ts), length(ts) >= 10)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)

  if (!requireNamespace("quantreg", quietly = TRUE)) stop("Package 'quantreg' is required.")

  d_ts <- trafo(ts, future, mode = level_mode, stride = TRUE)
  d_ts <- d_ts[is.finite(d_ts)]
  if (length(d_ts) < 10) stop("Transformed series too short after trafo().")

  idx <- seq.int(from = length(d_ts), to = 1L, by = -future)
  d_ts_ds <- d_ts[rev(idx)]
  n <- length(d_ts_ds)
  if (n < 10) stop("Downsampled series too short for quantile regression.")

  quants <- sort(unique(as.numeric(quants)))
  quants <- quants[quants > 0 & quants < 1]
  if (length(quants) < 3) stop("Need at least 3 quantiles strictly inside (0,1).")

  y <- as.numeric(d_ts_ds)

  if (jitter_ties && any(duplicated(y))) {
    set.seed(seed)
    eps <- 1e-10 * sd(y)
    if (!is.finite(eps) || eps == 0) eps <- 1e-10
    y <- y + rnorm(n, 0, eps)
  }

  dat <- data.frame(y = y, t = seq_len(n))
  fit <- rq(y ~ t, tau = quants, data = dat, method = rq_method)

  q_pred <- as.numeric(predict(fit, newdata = data.frame(t = n + 1)))
  q_pred <- cummax(q_pred)

  set.seed(seed)
  u <- runif(as.integer(n_draws))
  samples <- approx(x = quants, y = q_pred, xout = u, rule = 2)$y
  samples <- as.numeric(samples)
  samples <- samples[is.finite(samples)]
  if (length(samples) < 20) stop("Too few finite samples produced; check quantile fit.")

  base <- .emp_pred(samples, seed = seed)
  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}



#' @keywords internal

ewma_gaussian_predictor <- function(ts, future, lambda = 0.1, min_n = 10, level_mode = NULL) {
  stopifnot(is.numeric(ts), length(ts) >= 10)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  stopifnot(is.finite(lambda), lambda > 0, lambda <= 1)
  future <- as.integer(future)

  y <- trafo(ts, future, mode = level_mode, stride = TRUE)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < min_n) stop("Downsampled series too short for EWMA.")

  mu <- y[1]
  v <- 0

  for (i in 2:n) {
    mu_prev <- mu
    e <- y[i] - mu_prev               # innovation (pre-update)
    v <- lambda * e^2 + (1 - lambda) * v
    mu <- mu_prev + lambda * e
  }

  sigma <- sqrt(max(v, .Machine$double.eps))
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1e-8

  base <- list(
    rfun = function(m) rnorm(as.integer(m), mean = mu, sd = sigma),
    dfun = function(x) dnorm(as.numeric(x), mean = mu, sd = sigma),
    pfun = function(q) pnorm(as.numeric(q), mean = mu, sd = sigma),
    qfun = function(p) qnorm(pmin(pmax(as.numeric(p), 0), 1), mean = mu, sd = sigma)
  )

  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}



#' @keywords internal

historical_bootstrap_predictor <- function(ts, future, n_draws = 5000, decay = NULL, seed = 123, level_mode = NULL) {
  stopifnot(is.numeric(ts), length(ts) >= 10)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)

  y <- trafo(ts, future, mode = level_mode, stride = TRUE)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  if (length(y) < 10) stop("Downsampled series too short for bootstrap.")

  set.seed(seed)
  if (is.null(decay)) {
    draws <- sample(y, size = as.integer(n_draws), replace = TRUE)
  } else {
    stopifnot(is.finite(decay), decay > 0)
    n <- length(y)
    age <- rev(seq_len(n) - 1L)  # 0 for last, larger in the past
    w <- exp(-age / decay)
    w <- w / sum(w)
    draws <- sample(y, size = as.integer(n_draws), replace = TRUE, prob = w)
  }

  base <- .emp_pred(draws, seed = seed)
  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}


#' @keywords internal

drift_resid_bootstrap_predictor <- function(ts, future, n_draws = 5000, seed = 123, min_n = 12, level_mode = NULL) {
  stopifnot(is.numeric(ts), length(ts) >= 10)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)

  y <- trafo(ts, future, mode = level_mode, stride = TRUE)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < min_n) stop("Downsampled series too short for drift+residual bootstrap.")

  t <- seq_len(n)
  fit <- lm(y ~ t)
  mu_next <- as.numeric(predict(fit, newdata = data.frame(t = n + 1)))
  resid <- as.numeric(residuals(fit))
  resid <- resid[is.finite(resid)]
  if (length(resid) < 5) resid <- y - mean(y)

  set.seed(seed)
  draws <- mu_next + sample(resid, size = as.integer(n_draws), replace = TRUE)

  base <- .emp_pred(draws, seed = seed)
  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}


#' @keywords internal

vol_scaled_naive_predictor <- function(ts, future, vol_window = 20, min_n = 10, level_mode = NULL) {
  stopifnot(is.numeric(ts), length(ts) >= 10)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)
  vol_window <- as.integer(vol_window)

  y <- trafo(ts, future, mode = level_mode, stride = TRUE)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < min_n) stop("Downsampled series too short for vol-scaled naive.")

  mu <- as.numeric(tail(y, 1))
  w <- min(vol_window, n)
  sigma <- sd(tail(y, w))
  if (!is.finite(sigma) || sigma <= 0) sigma <- sd(y)
  if (!is.finite(sigma) || sigma <= 0) sigma <- 1e-8

  base <- list(
    rfun = function(m) rnorm(as.integer(m), mean = mu, sd = sigma),
    dfun = function(x) dnorm(as.numeric(x), mean = mu, sd = sigma),
    pfun = function(q) pnorm(as.numeric(q), mean = mu, sd = sigma),
    qfun = function(p) qnorm(pmin(pmax(as.numeric(p), 0), 1), mean = mu, sd = sigma)
  )

  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}


#' @keywords internal

robust_median_mad_predictor <- function(ts, future, window = 50, dist = c("laplace", "normal"), min_n = 10, level_mode = NULL) {
  dist <- match.arg(dist)
  stopifnot(is.numeric(ts), length(ts) >= 10)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)
  window <- as.integer(window)

  y <- trafo(ts, future, mode = level_mode, stride = TRUE)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < min_n) stop("Downsampled series too short for robust predictor.")

  w <- min(window, n)
  yy <- tail(y, w)
  mu <- median(yy)
  mad <- mad(yy, constant = 1)
  if (!is.finite(mad) || mad <= 0) mad <- mad(y, constant = 1)
  if (!is.finite(mad) || mad <= 0) mad <- 1e-8

  if (dist == "normal") {
    sigma <- mad / 0.67448975
    sigma <- max(sigma, 1e-8)

    base <- list(
      rfun = function(m) rnorm(as.integer(m), mean = mu, sd = sigma),
      dfun = function(x) dnorm(as.numeric(x), mean = mu, sd = sigma),
      pfun = function(q) pnorm(as.numeric(q), mean = mu, sd = sigma),
      qfun = function(p) qnorm(pmin(pmax(as.numeric(p), 0), 1), mean = mu, sd = sigma)
    )

    return(.apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode))
  }

  # Laplace
  b <- mad / log(2)
  b <- max(b, 1e-8)

  rlap <- function(n) {
    u <- runif(as.integer(n)) - 0.5
    mu - b * sign(u) * log(1 - 2 * abs(u))
  }
  dlap <- function(x) (1 / (2 * b)) * exp(-abs(as.numeric(x) - mu) / b)
  plap <- function(q) {
    q <- as.numeric(q)
    ifelse(q < mu, 0.5 * exp((q - mu) / b), 1 - 0.5 * exp(-(q - mu) / b))
  }
  qlap <- function(p) {
    p <- pmin(pmax(as.numeric(p), 0), 1)
    ifelse(p < 0.5, mu + b * log(2 * p), mu - b * log(2 * (1 - p)))
  }

  base <- list(rfun = rlap, dfun = dlap, pfun = plap, qfun = qlap)
  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}


#' @keywords internal

shrunk_quantile_predictor <- function(
    ts, future,
    taus = c(0.1, 0.25, 0.5, 0.75, 0.9),
    n_draws = 2000,
    seed = 123,
    rq_method = "fn",
    jitter_ties = TRUE,
    min_n = 12,
    level_mode = NULL
) {
  stopifnot(is.numeric(ts), length(ts) >= 10)
  stopifnot(length(future) == 1, is.finite(future), future >= 1)
  future <- as.integer(future)

  if (!requireNamespace("quantreg", quietly = TRUE)) stop("Package 'quantreg' is required.")

  y <- trafo(ts, future, mode = level_mode, stride = TRUE)
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < min_n) stop("Downsampled series too short for shrunk quantile predictor.")

  taus <- sort(unique(as.numeric(taus)))
  taus <- taus[taus > 0 & taus < 1]
  if (length(taus) < 3) stop("Need at least 3 taus strictly inside (0,1).")

  if (jitter_ties && any(duplicated(y))) {
    set.seed(seed)
    eps <- 1e-10 * sd(y)
    if (!is.finite(eps) || eps == 0) eps <- 1e-10
    y <- y + rnorm(n, 0, eps)
  }

  dat <- data.frame(y = y, t = seq_len(n))
  fit <- rq(y ~ t, tau = taus, data = dat, method = rq_method)

  q_pred <- as.numeric(predict(fit, newdata = data.frame(t = n + 1)))
  q_pred <- cummax(q_pred)

  set.seed(seed)
  u <- runif(as.integer(n_draws))
  draws <- approx(x = taus, y = q_pred, xout = u, rule = 2)$y

  base <- .emp_pred(draws, seed = seed)
  .apply_level_mode(base, last_level = tail(ts, 1), level_mode = level_mode)
}
#'
#'
#' @keywords internal
#'
plot_wired <- function(
    res_by_h,
    ts_set,
    x = NULL,                       # optional date/time vector length(ts)
    x_step_days = 1,                # used only if x is Date/POSIXct and step can't be inferred
    probs = c(0.05, 0.25, 0.50, 0.75, 0.95),
    history_n = 200,
    main_prefix = "Forecast",
    xlab = NULL,
    ylab = "level",
    ylab_cex = 0.8,
    lwd_hist = 1.5,
    lwd_med = 2,
    col_hist = "black",
    col_med = "blue",
    band_col = "lightblue",
    band_alpha = 0.25,
    add_points = TRUE,
    pch_hist = 16,
    cex_hist = 0.4,
    # --- x coercion controls ---
    x_tz = "UTC",
    x_try_epoch = TRUE,
    epoch_threshold = 1e9           # ~2001-09-09 in seconds; above this might be epoch seconds
) {
  # ---- validation ----
  if (!is.list(res_by_h) || length(res_by_h) < 1) stop("res_by_h must be a non-empty list of horizon results.")
  H <- length(res_by_h)

  if (is.null(res_by_h[[1]]$marginals_level)) {
    stop("LEVEL bands requested, but res_by_h[[1]]$marginals_level is NULL. Run muvar_prob_copula with mode != NULL.")
  }

  # normalize ts_set + names
  if (is.matrix(ts_set) || is.data.frame(ts_set)) {
    ts_set <- as.data.frame(ts_set)
    ts_set <- lapply(ts_set, as.numeric)
  } else if (is.list(ts_set)) {
    ts_set <- lapply(ts_set, as.numeric)
  } else {
    stop("ts_set must be a list of numeric vectors, or a matrix/data.frame.")
  }
  p <- length(ts_set)
  if (p < 1) stop("No series in ts_set.")
  if (is.null(names(ts_set))) names(ts_set) <- paste0("S", seq_len(p))

  # ============================================================
  # Robust x coercion for correct date/time plotting
  # ============================================================
  coerce_x <- function(x) {
    if (is.null(x)) return(NULL)

    if (inherits(x, "POSIXlt")) x <- as.POSIXct(x, tz = x_tz)

    if (inherits(x, "Date") || inherits(x, "POSIXct")) return(x)

    if (is.character(x)) {
      x2 <- as.Date(x)
      if (!all(is.na(x2))) return(x2)
      x3 <- as.POSIXct(x, tz = x_tz)
      if (!all(is.na(x3))) return(x3)
      stop("x is character but could not be parsed as Date or POSIXct. Provide ISO dates or a proper Date/POSIXct vector.")
    }

    if (is.numeric(x)) {
      if (x_try_epoch && all(is.finite(x)) && median(x, na.rm = TRUE) > epoch_threshold) {
        return(as.POSIXct(x, origin = "1970-01-01", tz = x_tz))
      }
      return(as.numeric(x))
    }

    stop("x must be NULL, Date, POSIXct/POSIXlt, character dates, or numeric.")
  }

  x <- coerce_x(x)
  # IMPORTANT: do NOT as.vector(x) here, it strips Date/POSIXct class!

  # sanity: monotone increasing helps date plots
  if (!is.null(x) && length(x) >= 2) {
    dx <- diff(as.numeric(x))
    if (any(!is.finite(dx))) warning("x contains non-finite differences; plotting may be incorrect.")
    if (any(dx <= 0, na.rm = TRUE)) warning("x is not strictly increasing. Date/time axis may look wrong.")
  }

  # extend x into forecast horizon
  extend_x <- function(x_full, H) {
    x_last <- tail(x_full, 1)

    if (inherits(x_full, "Date")) {
      if (length(x_full) >= 2) {
        step <- as.integer(tail(x_full, 1) - x_full[length(x_full) - 1])
        if (!is.finite(step) || step == 0) step <- as.integer(x_step_days)
      } else step <- as.integer(x_step_days)
      return(x_last + seq_len(H) * step)
    }

    if (inherits(x_full, "POSIXct")) {
      if (length(x_full) >= 2) {
        step <- as.numeric(tail(x_full, 1) - x_full[length(x_full) - 1])
        if (!is.finite(step) || step == 0) step <- as.numeric(x_step_days) * 86400
      } else step <- as.numeric(x_step_days) * 86400
      return(x_last + seq_len(H) * step)
    }

    if (length(x_full) >= 2) {
      step <- tail(x_full, 1) - x_full[length(x_full) - 1]
      if (!is.finite(step) || step == 0) step <- 1
    } else step <- 1
    x_last + seq_len(H) * step
  }

  # probs sanity
  probs <- sort(unique(as.numeric(probs)))
  if (any(!is.finite(probs)) || any(probs <= 0) || any(probs >= 1)) stop("probs must be in (0,1).")
  if (!any(abs(probs - 0.5) < 1e-12)) stop("probs must include 0.50.")

  # alpha helper (base)
  alpha_col <- function(col, alpha) {
    rgbv <- col2rgb(col) / 255
    rgb(rgbv[1, ], rgbv[2, ], rgbv[3, ], alpha = alpha)
  }
  band_fill <- alpha_col(band_col, band_alpha)

  # pretty date axis helper
  draw_time_axis <- function(x_range) {
    if (inherits(x_range, "Date")) {
      at <- pretty(x_range)
      axis.Date(1, at = at, format = "%d-%b", las = 1)
      return(invisible())
    }
    if (inherits(x_range, "POSIXct")) {
      at <- pretty(x_range)
      axis.POSIXct(1, at = at, format = "%d-%b", tz = x_tz, las = 1)
      return(invisible())
    }
    # numeric: default axis
    axis(1)
    invisible()
  }

  # layout + more left margin
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  nrow <- ceiling(sqrt(p))
  ncol <- ceiling(p / nrow)
  par(mfrow = c(nrow, ncol), mar = c(3.5, 4.8, 2.5, 1.0))

  q_store <- vector("list", p)
  names(q_store) <- names(ts_set)

  for (j in seq_len(p)) {
    ts <- as.numeric(ts_set[[j]])
    Tn <- length(ts)
    if (Tn < 2) stop(sprintf("Series '%s' too short.", names(ts_set)[j]))

    # per-series x (supports differing lengths)
    if (is.null(x)) {
      x_full <- seq_len(Tn)
    } else {
      if (length(x) < Tn) stop("x length is shorter than at least one series length.")
      x_full <- x[seq_len(Tn)]  # keeps Date/POSIXct class
    }

    idx0 <- max(1L, Tn - as.integer(history_n) + 1L)
    ts_hist <- ts[idx0:Tn]
    x_hist  <- x_full[idx0:Tn]

    x_fc <- extend_x(x_full, H)

    # quantiles per horizon
    Q <- matrix(NA_real_, nrow = H, ncol = length(probs))
    for (h in seq_len(H)) {
      rh <- res_by_h[[h]]
      if (is.null(rh$marginals_level) || is.null(rh$marginals_level[[j]]) || is.null(rh$marginals_level[[j]]$qfun)) {
        stop(sprintf("Missing LEVEL qfun for series '%s' at horizon h=%d.", names(ts_set)[j], h))
      }
      Q[h, ] <- as.numeric(rh$marginals_level[[j]]$qfun(probs))
    }

    y_all <- c(ts_hist, as.numeric(Q))
    y_rng <- range(y_all[is.finite(y_all)], na.rm = TRUE)

    ttl <- paste0(main_prefix, " - ", names(ts_set)[j])
    ylab_i <- if (is.null(ylab)) names(ts_set)[j] else ylab
    xlab_i <- xlab
    if (is.null(xlab_i)) xlab_i <- if (inherits(x_full, "Date") || inherits(x_full, "POSIXct")) "date" else "t"

    # if time-based x, suppress default x-axis and draw a formatted one
    is_time_x <- inherits(x_full, "Date") || inherits(x_full, "POSIXct")

    plot(x_hist, ts_hist, type = "l",
         main = ttl, xlab = xlab_i, ylab = ylab_i,
         lwd = lwd_hist, col = col_hist,
         xlim = range(c(x_hist, x_fc)),
         ylim = y_rng,
         cex.lab = ylab_cex,
         xaxt = if (is_time_x) "n" else "s")

    if (is_time_x) draw_time_axis(range(c(x_hist, x_fc)))

    if (isTRUE(add_points)) points(x_hist, ts_hist, pch = pch_hist, cex = cex_hist, col = col_hist)

    # bands
    med_idx <- which.min(abs(probs - 0.5))
    med <- Q[, med_idx]

    lower_probs <- sort(probs[probs < 0.5])
    for (lp in lower_probs) {
      up <- 1 - lp
      li <- which.min(abs(probs - lp))
      ui <- which.min(abs(probs - up))
      lo <- Q[, li]
      hi <- Q[, ui]

      polygon(
        x = c(x_fc, rev(x_fc)),
        y = c(lo, rev(hi)),
        col = band_fill,
        border = NA
      )
    }

    lines(x_fc, med, col = col_med, lwd = lwd_med)
    segments(x_full[Tn], ts[Tn], x_fc[1], med[1], col = col_med, lwd = 1)

    q_store[[j]] <- list(probs = probs, Q = Q, x_fc = x_fc, x_hist = x_hist, ts_hist = ts_hist)
  }

  p <- recordPlot()
  return(p)
}

#'
#'
#' @keywords internal
#'
check_feasibility <- function(
    ts_set,
    future,
    n_testing,
    min_train_frac = 1/3,
    min_points_per_model = 20
) {
  # normalize input
  if (is.matrix(ts_set) || is.data.frame(ts_set)) {
    ts_set <- as.data.frame(ts_set)
    ts_set <- lapply(ts_set, as.numeric)
  } else if (is.list(ts_set)) {
    ts_set <- lapply(ts_set, as.numeric)
  } else {
    stop("ts_set must be a list, matrix, or data.frame of numeric series.")
  }

  lens <- vapply(ts_set, length, integer(1))
  n_min <- min(lens)

  if (n_min < 20) {
    return(list(
      feasible = FALSE,
      confirmed = FALSE,
      suggested_n_testing = NA_integer_,
      suggested_future = NA_integer_,
      message = sprintf(
        "Series are too short (minimum length = %d). Need at least 20 observations.",
        n_min
      )
    ))
  }

  future <- as.integer(future)
  n_testing <- as.integer(n_testing)

  # helper: feasibility check for a given (future, n_testing)
  is_feasible_combo <- function(fut, ntest, n) {
    if (!is.finite(fut) || fut < 1) return(FALSE)
    if (!is.finite(ntest) || ntest < 1) return(FALSE)

    max_train_end <- n - fut
    if (max_train_end < 1) return(FALSE)

    start_idx <- floor(n / 3)
    if (start_idx > max_train_end) return(FALSE)

    test_idx <- round(seq.int(
      from = start_idx,
      to = max_train_end,
      length.out = ntest
    ))

    # enough training observations in earliest split
    earliest_end <- min(test_idx)
    earliest_train_len <- earliest_end

    # enough transformed points after lag/stride for downstream predictors
    enough_history <- (n - fut) >= min_points_per_model

    earliest_train_len >= max(20, ceiling(n * min_train_frac)) && enough_history
  }

  # 1) original combination feasible -> confirm
  if (is_feasible_combo(future, n_testing, n_min)) {
    return(list(
      feasible = TRUE,
      confirmed = TRUE,
      suggested_n_testing = n_testing,
      suggested_future = future,
      message = sprintf(
        "Parameters are feasible as provided: future = %d, n_testing = %d.",
        future, n_testing
      )
    ))
  }

  # 2) first try reducing n_testing, keeping future fixed
  max_ntesting_try <- max(1L, n_testing)
  for (nt in seq.int(from = max_ntesting_try, to = 1L, by = -1L)) {
    if (is_feasible_combo(future, nt, n_min)) {
      return(list(
        feasible = TRUE,
        confirmed = FALSE,
        suggested_n_testing = nt,
        suggested_future = future,
        message = sprintf(
          paste(
            "Original parameters are not feasible.",
            "Suggested feasible combination: future = %d, n_testing = %d.",
            "Priority rule applied: reduced n_testing first."
          ),
          future, nt
        )
      ))
    }
  }

  # 3) last resort: reduce future, then n_testing
  for (fut in seq.int(from = future - 1L, to = 1L, by = -1L)) {
    for (nt in seq.int(from = n_testing, to = 1L, by = -1L)) {
      if (is_feasible_combo(fut, nt, n_min)) {
        return(list(
          feasible = TRUE,
          confirmed = FALSE,
          suggested_n_testing = nt,
          suggested_future = fut,
          message = sprintf(
            paste(
              "Original parameters are not feasible.",
              "Suggested feasible combination: future = %d, n_testing = %d.",
              "Priority rule applied: reduced n_testing first, then future as last resort."
            ),
            fut, nt
          )
        ))
      }
    }
  }

  # 4) nothing workable
  list(
    feasible = FALSE,
    confirmed = FALSE,
    suggested_n_testing = NA_integer_,
    suggested_future = NA_integer_,
    message = sprintf(
      paste(
        "No feasible combination found, even after reducing n_testing and future.",
        "Minimum series length is %d."
      ),
      n_min
    )
  )
}
#'
#'
#'
old_check_feasibility <- function(
    ts_set,
    future,
    n_testing,
    min_train = 20,
    context = "ts_set"
) {
  # ---- normalize ts_set to data.frame ----
  if (is.matrix(ts_set) | is.list(ts_set)) ts_set <- as.data.frame(ts_set)
  if (!is.data.frame(ts_set)) stop("ts_set must be a data.frame (or list or matrix coercible to data.frame).")
  if (ncol(ts_set) < 1) stop("ts_set has no columns.")
  if (nrow(ts_set) < 2) stop("ts_set is too short (nrow < 2).")

  future <- as.integer(future)
  n_testing <- as.integer(n_testing)
  min_train <- as.integer(min_train)

  if (!is.finite(future) || future < 1) stop("future must be >= 1.")
  if (!is.finite(n_testing) || n_testing < 0) stop("n_testing must be >= 0.")
  if (!is.finite(min_train) || min_train < 3) stop("min_train must be >= 3.")

  N <- nrow(ts_set)

  # ---- compute effective transformed length L ----
  L <- floor((nrow(ts_set) - future)/future)

  if (!is.finite(L) || L <= 0) {
    stop(sprintf(
      "%s is too short for future=%d (effective transformed length L=%s).",
      context, future, as.character(L)
    ))
  }

  # ---- feasibility check based on your rule ----
  if (n_testing == 0) {
    if (L < min_train) {
      max_future <- max(1L, N - min_train)
      msg <- c(
        sprintf("%s is too short for future=%d with n_testing=0.", context, future),
        sprintf("Need L >= min_train, but L=%d and min_train=%d.", L, min_train),
        sprintf("Suggested maximum future for min_train=%d is: %d", min_train, max_future)
      )
      stop(paste(msg, collapse = "\n"))
    }
    return(invisible(list(
      ok = TRUE,
      N = N, L = L,
      min_win = L,
      max_future = max(1L, N - min_train),
      max_n_testing = NA_integer_
    )))
  }

  min_win <- floor(L / n_testing)

  if (min_win >= min_train) {
    # also return useful maxima
    max_future <- max(1L, N - min_train * n_testing)
    max_n_testing <- max(1L, floor(L / min_train))
    return(invisible(list(
      ok = TRUE,
      N = N, L = L,
      min_win = min_win,
      max_future = max_future,
      max_n_testing = max_n_testing
    )))
  }

  # ---- suggestions (solve the inequality) ----
  # Need: floor((N - future)/future/n_testing) >= min_train  ~ N / (1 + min_train*n_testing) >= future
  max_future <- N / (1 + min_train * n_testing)
  max_future <- max(1L, as.integer(max_future))

  max_n_testing <- floor(L / min_train)
  max_n_testing <- max(0L, as.integer(max_n_testing))

  msg <- c(
    sprintf("%s is too short for future=%d and n_testing=%d under the expanding-window rule.", context, future, n_testing),
    sprintf("Effective transformed length L=%d (N=%d, future=%d).", L, N, future),
    sprintf("Minimum training window size: floor(L / n_testing) = floor(%d / %d) = %d.", L, n_testing, min_win),
    sprintf("Feasibility requires min_win >= min_train, but min_train=%d.", min_train),
    "",
    "Try one of the following:",
    sprintf("  * Reduce future (h): maximum future for n_testing=%d is about: %d", n_testing, max_future),
    sprintf("  * Reduce n_testing: maximum n_testing for future=%d is about: %d", future, max_n_testing)
  )

  stop(paste(msg, collapse = "\n"))
}



