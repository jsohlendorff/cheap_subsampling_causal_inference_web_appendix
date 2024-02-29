get_cheap_subsampling_ci <-
  function(est, boot_est, m_val, n_val, alpha) {
    s_val <- sqrt(mean((est - boot_est)^2))
    tq <- qt(1 - alpha / 2, df = length(boot_est))
    list(
      lower_b = est - tq * sqrt((m_val) / (n_val - m_val)) * s_val,
      upper_b = est + tq * sqrt((m_val) / (n_val - m_val)) * s_val
    )
  }

get_cheap_bootstrap_ci <-
  function(est, boot_est, n_val, alpha) {
    b_val <- length(boot_est)
    s_val <- sqrt(mean((est - boot_est)^2))
    tq <- qt(1 - alpha / 2, df = b_val)
    list(
      lower_b = est - tq * s_val,
      upper_b = est + tq * s_val
    )
  }
