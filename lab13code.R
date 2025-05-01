library(pwr)
library(tidyverse)
library(patchwork)
library(effectsize)
library(e1071)
library(boot)
library(boot.pval)

finch.dat <- read_csv("zebrafinches.csv")

n <- nrow(finch.dat)
#######################################################################
#task 1
#######################################################################
####  A   ####

sum.farther.vals <- finch.dat |>
  summarise(
    mean = mean(further),
    sd = sd(further),
    skew = skewness(further)
  )
#T-test and t-stat
t.farther <- t.test(finch.dat$further, mu = 0, alternative = "less")

t.stat <- t.farther$statistic

fz <- dnorm(t.stat)


error <- ((sum.farther.vals$skew/sqrt(n))*((2*(t.stat)^2+1)/6)*fz)

####    B   ####

t.vals <- seq(-10, 10, length.out = 1000 )

FZ.vals <- dnorm(t.vals)

error.vals <- (sum.farther.vals$skew / sqrt(n)) * ((2 * t.vals^2 + 1) / 6) * FZ.vals

error.df <- tibble(
    t     = t.vals,
    error = error.vals
  )

edgeworth.graph <- ggplot(error.df, aes(x = t, y = error)) +
  geom_line(size = 1) +
  labs(
    title    = "Edgeworth Approximation Error vs t",
    subtitle = paste0("n = ", n, ", skewness = ", round(sum.farther.vals$skew, 3)),
    x        = "t",
    y        = "Edgeworth Error"
  ) +
  theme_minimal()

ggsave("edgeworth.pdf", plot = edgeworth.graph, height = 4, width= 7)

####    C   ####
alpha <- .05

theoretical.t.value  <- abs(qnorm(alpha))
fz                   <- dnorm(theoretical.t.value)


solved.sample.size <- (
  abs(sum.farther.vals$skew)
  / (6 * (0.10 * alpha) * (2 * theoretical.t.value^2 + 1) * fz)
)^2

#######################################################################
#######################################################################
#######################################################################
#######################################################################
#task 2

################################################################################
################################################################################

######
R <- 10000
# Helper to bootstrap t-statistics and shift to null
bootstrap_t_null <- function(x, R) {
  n  <- length(x)
  s  <- sd(x)
  
  # vector to store T-statistics
  tstats <- numeric(R)
  
  for (i in seq_len(R)) {
    samp    <- sample(x, size = n, replace = TRUE)
    xbar    <- mean(samp)
    # T = (x̄_r - 0) / (s / sqrt(n))
    tstats[i] <- (xbar - 0) / (s / sqrt(n))
  }
  # Shift so that mean(tstats) = 0 under H0
  tstats - mean(tstats)
}

#
resamples.null.closer <- bootstrap_t_null(finch.dat$closer, R)

resamples.null.further <- bootstrap_t_null(finch.dat$further, R)

resamples.null.diff    <- bootstrap_t_null(finch.dat$diff, R)


#########################
## b
#compute t stat
obs_t_closer <- (mean(finch.dat$closer)   - 0) / (sd(finch.dat$closer)   / sqrt(n))
obs_t_further <- (mean(finch.dat$further) - 0) / (sd(finch.dat$further) / sqrt(n))
obs_t_diff    <- (mean(finch.dat$diff) - 0) / (sd(finch.dat$diff) / sqrt(n))

#compute p-values for mean
p_boot_closer <- mean(abs(resamples.null.closer) >= abs(obs_t_closer))
p_boot_further <- mean(abs(resamples.null.further) >= abs(obs_t_further))
p_boot_diff    <- mean(abs(resamples.null.diff)    >= abs(obs_t_diff))

# 3. Compute standard two-sided t-test p-values
tt_closer <- t.test(finch.dat$closer)$p.value
tt_further <- t.test(finch.dat$further)$p.value
tt_diff    <- t.test(finch.dat$diff)$p.value

# 4. Summarize side-by-side
results <- data.frame(
  Test       = c("Closer", "Further", "Difference"),
  Obs_T      = c(obs_t_closer, obs_t_further, obs_t_diff),
  P_bootstrap = c(p_boot_closer, p_boot_further, p_boot_diff),
  P_t_test   = c(tt_closer, tt_further, tt_diff)
)



#########################
## c
df <- n - 1

# Empirical 5th percentiles
q5_closer <- quantile(resamples.null.closer, 0.05)
q5_further <- quantile(resamples.null.further, 0.05)
q5_diff    <- quantile(resamples.null.diff,    0.05)

# Theoretical t cutoff at α=0.05 (one-tailed)
t5_theoretical <- qt(0.05, df = df)

# Summarize
results_c <- data.frame(
  Test           = c("Closer", "Further", "Difference", "Theoretical"),
  `5th Percentile` = c(q5_closer, q5_further, q5_diff, t5_theoretical)
)

############################
##d

R <- 10000

# Helper to get percentile CI for the mean
bootstrap_ci_mean <- function(x, R) {
  means <- replicate(R, mean(sample(x, size = length(x), replace = TRUE)))
  quantile(means, c(0.025, 0.975))
}

# 1. Bootstrap percentile CIs
ci_boot_closer <- bootstrap_ci_mean(finch.dat$closer,   R)
ci_boot_further <- bootstrap_ci_mean(finch.dat$further, R)
ci_boot_diff    <- bootstrap_ci_mean(finch.dat$diff, R)

# 2. Classical t-test CIs
ci_t_closer   <- t.test(finch.dat$closer)$conf.int
ci_t_further  <- t.test(finch.dat$further)$conf.int
ci_t_diff     <- t.test(finch.dat$diff)$conf.int

# 3. Summarize comparison
ci_comparison <- data.frame(
  Test        = c("Closer", "Further", "Difference"),
  Boot_Lower  = c(ci_boot_closer[1], ci_boot_further[1], ci_boot_diff[1]),
  Boot_Upper  = c(ci_boot_closer[2], ci_boot_further[2], ci_boot_diff[2]),
  T_Lower     = c(ci_t_closer[1],    ci_t_further[1],    ci_t_diff[1]),
  T_Upper     = c(ci_t_closer[2],    ci_t_further[2],    ci_t_diff[2])
)



################################################################################
################################################################################
#By Hand 
################################################################################
# Bootstrap Test Confidence Interval Closer

R <- 100000
resamples <- tibble(xbars = rep(NA, R))
s <- sd(finch.dat$closer)  
n <- nrow(finch.dat)
for(i in 1:R){
  curr.resample <- sample(finch.dat$closer,
                          size = n,
                          replace = TRUE)
  resamples$xbars[i] <- (mean(curr.resample) - 0) / (s / sqrt(n))
}

ggdat <- tibble(x=seq(3, 13, length.out=1000))|>
  mutate(pdf = dnorm(x, mean(resamples$xbars), sd(resamples$xbars)))

ggplot()+
  geom_histogram(data=resamples, aes(x=xbars, y=after_stat(density)))+
  geom_line(data=ggdat, aes(x=x, y=pdf), color="red")

# Confidence Interval
quantile(resamples$xbars, c(0.025, 0.975))


# Plot the Bootstrap Hypothesis Test
# shift so H0 is true
mean(resamples$xbars)
(delta <- mean(resamples$xbars) - 0)


resamples <- resamples |>
  mutate(xbars.shifted = xbars - delta)

low <- 0 - delta
high <- 0 + delta

resamples |>
  summarize(mean = mean(xbars.shifted),
            p.low = mean(xbars.shifted <= low),
            p.high = mean(xbars.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()

# Plot the Bootstrap Test
ggdat.obs <- tibble(xbar=mean(finch.dat$closer),
                    y=0)|>
  mutate(mirror = mu0 - (xbar-mu0))

bootstrap.plot.closer <- ggplot() +
  # Plot Resampling Distribution
  stat_density(data = resamples, aes(x = xbars),            geom = "line", color = "lightgrey") +
  stat_density(data = resamples, aes(x = xbars.shifted),    geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  # Plot Observation
  geom_point(data = ggdat.obs, aes(x = xbar,   y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"),
             shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name   = bquote(bar(x)),
    limits = range(c(resamples$xbars, resamples$xbars.shifted)) * 1.1,
    breaks = pretty(c(resamples$xbars, resamples$xbars.shifted), n = 5)
  ) +
  ylab("Density") +
  ggtitle("Bootstrap Hypothesis Test")

################################################################################
################################################################################
# Bootstrap Test Confidence Interval Further 
R <- 100000
resamples <- tibble(xbars = rep(NA, R))
s <- sd(finch.dat$further)  # for "further", adjust for "further" and "diff"
n <- nrow(finch.dat)
for(i in 1:R){
  curr.resample <- sample(finch.dat$closer,
                          size = n,
                          replace = TRUE)
  resamples$xbars[i] <- (mean(curr.resample) - 0) / (s / sqrt(n))
}


# Confidence Interval
quantile(resamples$xbars, c(0.025, 0.975))


# Plot the Bootstrap Hypothesis Test
# shift so H0 is true
obs.mean <- mean(finch.dat$further)
(delta <- obs.mean - 0)

resamples <- resamples |>
  mutate(xbars.shifted = xbars - delta)

low  <- min(0 - delta, 0 + delta)
high <- max(0 - delta, 0 + delta)

resamples |>
  summarize(mean = mean(xbars.shifted),
            p.low = mean(xbars.shifted <= low),
            p.high = mean(xbars.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()


# Plot the Bootstrap Test
ggdat.obs <- tibble(xbar=mean(finch.dat$further),
                    y=0)|>
  mutate(mirror = mu0 - (xbar-mu0))

bootstrap.plot.further <- ggplot() +
  # Plot Resampling Distribution
  stat_density(data = resamples, aes(x = xbars),            geom = "line", color = "lightgrey") +
  stat_density(data = resamples, aes(x = xbars.shifted),    geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  # Plot Observation
  geom_point(data = ggdat.obs, aes(x = xbar,   y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"),
             shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name   = bquote(bar(x)),
    limits = range(c(resamples$xbars, resamples$xbars.shifted)) * 1.1,
    breaks = pretty(c(resamples$xbars, resamples$xbars.shifted), n = 5)
  ) +
  ylab("Density") +
  ggtitle("Bootstrap Hypothesis Test")

################################################################################
################################################################################
# Bootstrap Test Confidence Interval Difference
R <- 100000
resamples <- tibble(xbars = rep(NA, R))

for(i in 1:R){
  curr.resample <- sample(finch.dat$diff,
                          size = nrow(dat),
                          replace = T)
  
  resamples$xbars[i] <- mean(curr.resample)
}

# Confidence Interval
quantile(resamples$xbars, c(0.025, 0.975))


################################################################################
# Plot the Bootstrap Hypothesis Test
# shift so H0 is true
mean(resamples$xbars)
(delta <- mean(resamples$xbars) - 0)


resamples <- resamples |>
  mutate(xbars.shifted = xbars - delta)

low <- 0 - delta
high <- 0 + delta

resamples |>
  summarize(mean = mean(xbars.shifted),
            p.low = mean(xbars.shifted <= low),
            p.high = mean(xbars.shifted >= high))|>
  mutate(p = p.low + p.high) |>
  view()

################################################################################
# Plot the Bootstrap Test
ggdat.obs <- tibble(xbar=mean(finch.dat$diff),
                    y=0)|>
  mutate(mirror = mu0 - (xbar-mu0))

bootstrap.plot.diff <- ggplot() +
  # Plot Resampling Distribution
  stat_density(data = resamples, aes(x = xbars),            geom = "line", color = "lightgrey") +
  stat_density(data = resamples, aes(x = xbars.shifted),    geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  # Plot Observation
  geom_point(data = ggdat.obs, aes(x = xbar,   y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"),
             shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name   = bquote(bar(x)),
    limits = range(c(resamples$xbars, resamples$xbars.shifted)) * 1.1,
    breaks = pretty(c(resamples$xbars, resamples$xbars.shifted), n = 5)
  ) +
  ylab("Density") +
  ggtitle("Bootstrap Hypothesis Test")


################################################################################
################################################################################
################################################################################
#task 3
################################################################################
#Randomization test for Closer
R   <- 10000
mu0 <- 0

#— Build null distribution for "closer"
rand <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- finch.dat$closer - mu0

# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift * 
    sample(x      = c(-1, 1),
           size   = length(x.shift),
           replace = TRUE)
  
  rand$xbars[i] <- mean(curr.rand)
}

# shift back
rand <- rand |> mutate(xbars = xbars + mu0)

#— Compute two‐sided p‐value
(delta <- abs(mean(finch.dat$closer) - mu0))
(low   <- mu0 - delta)  # mirror
(high  <- mu0 + delta)  # observed x̄

p_value_closer <- 
  mean(rand$xbars <= low) +
  mean(rand$xbars >= high)

p_value_closer  # randomization p-value for closer

# Plot the Randomization Test for "closer"
ggdat.obs <- tibble(
  xbar   = mean(finch.dat$closer),
  y      = 0
) |> mutate(mirror = mu0 - (xbar - mu0))

randomization.plot.closer <- ggplot() +
  stat_density(data = rand, aes(x = xbars), geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  geom_point(data = ggdat.obs, aes(x = xbar,   y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"),
             shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name   = bquote(bar(x)),
    limits = range(rand$xbars) * 1.1,         # ±10% padding around your actual null
    breaks = pretty(rand$xbars, n = 5)        # 5 “pretty” tick marks
  ) +
  ylab("Density") +
  ggtitle("Randomization Hypothesis Test (closer)")


# Randomization Test Confidence Interval for "closer"
R             <- 1000
mu0.iterate   <- 0.0001
starting.point <- mean(finch.dat$closer)

#— Lower bound search
mu.lower <- starting.point
repeat {
  rand <- tibble(xbars = rep(NA, R))
  x.shift <- finch.dat$closer - mu.lower
  
  for(i in 1:R){
    curr.rand     <- x.shift * 
      sample(c(-1, 1), length(x.shift), replace = TRUE)
    rand$xbars[i] <- mean(curr.rand)
  }
  rand <- rand |> mutate(xbars = xbars + mu.lower)
  
  delta <- abs(mean(finch.dat$closer) - mu.lower)
  low   <- mu.lower - delta
  high  <- mu.lower + delta
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  if (p.val < 0.05) break
  mu.lower <- mu.lower - mu0.iterate
}

#— Upper bound search
mu.upper <- starting.point
repeat {
  rand <- tibble(xbars = rep(NA, R))
  x.shift <- finch.dat$closer - mu.upper
  
  for(i in 1:R){
    curr.rand     <- x.shift * 
      sample(c(-1, 1), length(x.shift), replace = TRUE)
    rand$xbars[i] <- mean(curr.rand)
  }
  rand <- rand |> mutate(xbars = xbars + mu.upper)
  
  delta <- abs(mean(finch.dat$closer) - mu.upper)
  low   <- mu.upper - delta
  high  <- mu.upper + delta
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  if (p.val < 0.05) break
  mu.upper <- mu.upper + mu0.iterate
}

ciR_closer <- c(mu.lower, mu.upper)
ciR_closer  # randomization 95% CI for closer







################################################################################
################################################################################
#Randomization test for Further
R   <- 10000
mu0 <- 0

#— Build null distribution for "further"
rand <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- finch.dat$further - mu0

# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift * 
    sample(x      = c(-1, 1),
           size   = length(x.shift),
           replace = TRUE)
  
  rand$xbars[i] <- mean(curr.rand)
}

# shift back
rand <- rand |> mutate(xbars = xbars + mu0)

#— Compute two‐sided p‐value
(delta <- abs(mean(finch.dat$further) - mu0))
(low   <- mu0 - delta)  # mirror
(high  <- mu0 + delta)  # observed x̄

p_value_further <- 
  mean(rand$xbars <= low) +
  mean(rand$xbars >= high)

p_value_further  # randomization p-value for further

# Plot the Randomization Test for "further"
ggdat.obs <- tibble(
  xbar   = mean(finch.dat$further),
  y      = 0
) |> mutate(mirror = mu0 - (xbar - mu0))

randomization.plot.further <- ggplot() +
  stat_density(data = rand, aes(x = xbars), geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  geom_point(data = ggdat.obs, aes(x = xbar,   y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"),
             shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name   = bquote(bar(x)),
    limits = range(rand$xbars) * 1.1,         # ±10% padding around your actual null
    breaks = pretty(rand$xbars, n = 5)        # 5 “pretty” tick marks
  ) +
  ylab("Density") +
  ggtitle("Randomization Hypothesis Test (Fruther)")


# Randomization Test Confidence Interval for "further"
R             <- 1000
mu0.iterate   <- 0.0001
starting.point <- mean(finch.dat$further)

#— Lower bound search
mu.lower <- starting.point
repeat {
  rand <- tibble(xbars = rep(NA, R))
  x.shift <- finch.dat$further - mu.lower
  
  for(i in 1:R){
    curr.rand     <- x.shift * 
      sample(c(-1, 1), length(x.shift), replace = TRUE)
    rand$xbars[i] <- mean(curr.rand)
  }
  rand <- rand |> mutate(xbars = xbars + mu.lower)
  
  delta <- abs(mean(finch.dat$further) - mu.lower)
  low   <- mu.lower - delta
  high  <- mu.lower + delta
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  if (p.val < 0.05) break
  mu.lower <- mu.lower - mu0.iterate
}

#— Upper bound search
mu.upper <- starting.point
repeat {
  rand <- tibble(xbars = rep(NA, R))
  x.shift <- finch.dat$further - mu.upper
  
  for(i in 1:R){
    curr.rand     <- x.shift * 
      sample(c(-1, 1), length(x.shift), replace = TRUE)
    rand$xbars[i] <- mean(curr.rand)
  }
  rand <- rand |> mutate(xbars = xbars + mu.upper)
  
  delta <- abs(mean(finch.dat$further) - mu.upper)
  low   <- mu.upper - delta
  high  <- mu.upper + delta
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  if (p.val < 0.05) break
  mu.upper <- mu.upper + mu0.iterate
}

ciR_further <- c(mu.lower, mu.upper)
ciR_further  # randomization 95% CI for further

################################################################################
################################################################################
#Randomization test for diff
R   <- 10000
mu0 <- 0

#— Build null distribution for "diff"
rand <- tibble(xbars = rep(NA, R))

# PREPROCESSING: shift the data to be mean 0 under H0
x.shift <- finch.dat$diff - mu0

# RANDOMIZE / SHUFFLE
for(i in 1:R){
  curr.rand <- x.shift * 
    sample(x      = c(-1, 1),
           size   = length(x.shift),
           replace = TRUE)
  
  rand$xbars[i] <- mean(curr.rand)
}

# shift back
rand <- rand |> mutate(xbars = xbars + mu0)

#— Compute two‐sided p‐value
(delta <- abs(mean(finch.dat$diff) - mu0))
(low   <- mu0 - delta)  # mirror
(high  <- mu0 + delta)  # observed x̄

p_value_diff <- 
  mean(rand$xbars <= low) +
  mean(rand$xbars >= high)

p_value_diff  # randomization p-value for diff

# Plot the Randomization Test for "diff"
ggdat.obs <- tibble(
  xbar   = mean(finch.dat$diff),
  y      = 0
) |> mutate(mirror = mu0 - (xbar - mu0))

randomization.plot.closer.diff <- ggplot() +
  stat_density(data = rand, aes(x = xbars), geom = "line", color = "black") +
  geom_hline(yintercept = 0) +
  geom_point(data = ggdat.obs, aes(x = xbar,   y = y, color = "Observation")) +
  geom_point(data = ggdat.obs, aes(x = mirror, y = y, color = "Mirror Observation"),
             shape = 19) +
  scale_color_manual("", values = c("black", "red")) +
  theme_bw() +
  scale_x_continuous(
    name   = bquote(bar(x)),
    limits = range(rand$xbars) * 1.1,         # ±10% padding around your actual null
    breaks = pretty(rand$xbars, n = 5)        # 5 “pretty” tick marks
  ) +
  ylab("Density") +
  ggtitle("Randomization Hypothesis Test (Fruther)")


# Randomization Test Confidence Interval for "diff"
R             <- 1000
mu0.iterate   <- 0.0001
starting.point <- mean(finch.dat$diff)

#— Lower bound search
mu.lower <- starting.point
repeat {
  rand <- tibble(xbars = rep(NA, R))
  x.shift <- finch.dat$diff - mu.lower
  
  for(i in 1:R){
    curr.rand     <- x.shift * 
      sample(c(-1, 1), length(x.shift), replace = TRUE)
    rand$xbars[i] <- mean(curr.rand)
  }
  rand <- rand |> mutate(xbars = xbars + mu.lower)
  
  delta <- abs(mean(finch.dat$diff) - mu.lower)
  low   <- mu.lower - delta
  high  <- mu.lower + delta
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  if (p.val < 0.05) break
  mu.lower <- mu.lower - mu0.iterate
}

#— Upper bound search
mu.upper <- starting.point
repeat {
  rand <- tibble(xbars = rep(NA, R))
  x.shift <- finch.dat$diff - mu.upper
  
  for(i in 1:R){
    curr.rand     <- x.shift * 
      sample(c(-1, 1), length(x.shift), replace = TRUE)
    rand$xbars[i] <- mean(curr.rand)
  }
  rand <- rand |> mutate(xbars = xbars + mu.upper)
  
  delta <- abs(mean(finch.dat$diff) - mu.upper)
  low   <- mu.upper - delta
  high  <- mu.upper + delta
  p.val <- mean(rand$xbars <= low) +
    mean(rand$xbars >= high)
  
  if (p.val < 0.05) break
  mu.upper <- mu.upper + mu0.iterate
}

ciR_diff <- c(mu.lower, mu.upper)
ciR_diff  # randomization 95% CI for diff



























