rm(list = ls())

simulate_data <- function(n_group, n_pergroup, n_fixef, var_rand, var_resid) {
  b <- rnorm(n = n_group, mean = 0, sd = sqrt(var_rand))
  intercept <- ifelse(n_fixef > 0, runif(n = 1, min = -5, max = 5), 0)
  
  fix_effects <- intercept
  
  if (n_fixef > 1) {
    slopes <- runif(n = n_fixef - 1, min = -5, max = 5)
    X <- replicate(n_fixef - 1, runif(n = n_group * n_pergroup))
    fix_effects <- as.vector(fix_effects + X %*% slopes)
  }
  
  d <- data.frame(y = rnorm(n = n_group * n_pergroup,
                            mean = fix_effects + rep(b, each = n_pergroup),
                            sd = sqrt(var_resid)),
                   g = gl(n_group, n_pergroup))
  
  if (n_fixef > 1) {
    d <- cbind(d, as.data.frame(X))
  }
  
  d
}

test_method <- function(n_group, n_fixef, var_rand, var_resid, n_replicate = 100, N = 100) {

  n_pergroup <- round(N / n_group)
  N <- n_group * n_pergroup # exact N may differ from input N due to previous line
  
  test_lambda <- replicate(n_replicate, {
    d <- simulate_data(n_group =  n_group, n_pergroup = n_pergroup,
                       n_fixef = n_fixef, var_rand = var_rand, var_resid = var_resid)
    if (n_fixef > 1) {
      form <- formula(y ~ . -g + (1|g))
    } else if (n_fixef == 1) {
      form <- formula(y ~ 1 + (1|g))
    } else {
      form <- formula(y ~ 0 + (1|g))
    }
    fit_ML   <- spaMM::fitme(form, data = d, method = "ML")
    fit_REML <- spaMM::fitme(form, data = d, method = "REML")
    
    c(ML = get_ranPars(fit_ML)$lambda[[1]],
      REML = get_ranPars(fit_REML)$lambda[[1]])
  })
  
  c(bias = apply(test_lambda, 1, function(x) abs(mean(x) - var_rand)),
    var = apply(test_lambda, 1, var),
    n_group =  n_group,
    n_pergroup = n_pergroup,
    n_fixef = n_fixef,
    var_rand = var_rand,
    var_resid = var_resid,
    n_replicate = n_replicate,
    N = N)
}

library(tidyverse)
library(patchwork)

settings <- expand.grid(n_group = c(2, 20),
                        n_fixef = c(2, 20),
                        var_rand = 2,
                        var_resid = 1,
                        n_replicate = 30)

replicates <- rep(1:nrow(settings), each = 10)
settings <- settings[replicates, ]

set.seed(1)
test <- apply(settings, 1, function(x) do.call(test_method, as.list(x)))

as.data.frame(t(test)) %>%
  mutate(across(.cols = contains("_"), .fns = as.factor)) %>%
  ggplot() +
    aes(y = bias.REML, x = bias.ML,
        size = n_group, shape = n_fixef, colour = n_fixef) +
    geom_point() +
    geom_abline(slope = 1, linetype = "dashed") +
    scale_shape_manual(values = seq_len(length(unique(test["n_group", ])))) +
    theme_classic() +
    theme(legend.position = "none", panel.grid = element_blank()) -> p1

as.data.frame(t(test)) %>%
  mutate(across(.cols = contains("_"), .fns = as.factor)) %>%
  ggplot() +
  aes(y = var.REML, x = var.ML,
      size = n_group, shape = n_fixef, colour = n_fixef) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed") +
  scale_shape_manual(values = seq_len(length(unique(test["n_group", ])))) +
  theme_classic() +
  labs(size = "Nb of levels\n in random factor",
       shape = "Nb of fixed-effect\n parameters",
       colour = "Nb of fixed-effect\n parameters") +
  theme(panel.grid = element_blank()) -> p2

p1 + p2 +
  plot_annotation(title = expression("Bias and variance of "~widehat(lambda)~"in LMM fited by REML & ML"),
                  subtitle = expression("(N = 100,"~phi~"= 1, "~lambda~"= 3)")) 

ggsave("img/simu_MLvsREML.png", width = 7, height = 3)
