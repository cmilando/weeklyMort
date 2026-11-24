library(data.table)
library(ggplot2)

# df <- read.csv("heat_mortality_asofnov6.csv")
# setDT(df)
# 
# head(df)
# 
# 
# # just grab one
# phx <- subset(df, occurrence_county_code == "4013")
# saveRDS(phx, "phx.RDS")

phx <- setDT(readRDS("phx.RDS"))
phx[, dt_fmt := as.IDate(X_week_end_date)]

#
ggplot(phx) +
  geom_line(aes(x = dt_fmt, y = X_deaths))

ggplot(phx) +
  geom_line(aes(x = dt_fmt, y = mean_temp_C))

ggplot(phx) +
  geom_histogram(aes(x = mean_temp_C))
quantile(phx$mean_temp_C, probs = 0.95)

# make strata of county-year-month
phx[, strata := paste0(occurrence_county_code, 
                       ".", year(dt_fmt), ".",
                       month(dt_fmt))]
head(phx)


# make a baseline trend
phx[, yday := yday(dt_fmt)]

season_sine <- function(x) {
  A = 200 # amplitude
  B = 2*pi/365.25 # frequency in days
  C = 4.44 # phase shift
  D = 1000 # vertical shift
  A * sin(B*x + C) + D
}
phx[,sim_death := 10 * season_sine(yday)]
# phx[,sim_death := 10 * 1000]

ggplot(phx) +
  geom_line(aes(x = dt_fmt, y = sim_death))

# add noise based on sd
# applied in log space to make sure its positive
rand = 0.00
noise = runif(nrow(phx), min = -rand, max = rand)
phx[, sim_death := exp(log(sim_death) + noise)]

# round to integer
phx[, sim_death := round(sim_death)]
phx[, sim_death_baseline := sim_death]

ggplot(phx) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline))

# now introduce heat response
# doing this just with one variable not with crossbasis

phx[, gt40 := 1*(mean_temp_C > 40)]
phx

beta_temp = log(2)

phx[, sim_death := sim_death_baseline * exp(gt40 * beta_temp)]
phx

ggplot(subset(phx, year(dt_fmt) == 2018)) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline), col = 'grey') +
  geom_line(aes(x = dt_fmt, y = sim_death)) +
  geom_point(aes(x = dt_fmt, y = sim_death, color = factor(gt40)))

# 
library(gnm)

m1 <- gnm(sim_death ~ gt40, 
          data = phx,
          family = quasipoisson,
          eliminate = factor(strata))

m1

coef(m1)
exp(coef(m1))

#
ggplot(subset(phx, year(dt_fmt) == 2018)) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline), col = 'grey') +
  geom_line(aes(x = dt_fmt, y = sim_death)) +
  geom_point(aes(x = dt_fmt, y = sim_death, color = factor(gt40))) + 
  facet_wrap(~strata)


#
library(splines)
phx[, yr := year(dt_fmt)]
nyears = length(unique(phx$yr))
nyears
df_per_year = 8

m1.mod <- gnm(sim_death ~ gt40 + 
                ns(dt_fmt,df = df_per_year * nyears), 
              data = phx,
              family = quasipoisson)

m1.mod

# =============================================================================
# /////////////////////////////////////////////////////////////////////////////
# TRIAL #1: COVID EFFECT IS ADDITIVE
#
# log(death) = log(baseline) + beta_temp * gt40 + beta_covid * is_covid_yr
#
# /////////////////////////////////////////////////////////////////////////////
# =============================================================================

beta_covid = log(2)

phx[, is_covid_year := year(dt_fmt) %in% c(2020, 2021)]

phx[, sim_death_t1 := sim_death_baseline * 
      exp(gt40 * beta_temp + is_covid_year * beta_covid)
]

# plot
ggplot(subset(phx, year(dt_fmt) == 2020)) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline), col = 'grey') +
  geom_line(aes(x = dt_fmt, y = sim_death_t1)) +
  geom_point(aes(x = dt_fmt, y = sim_death_t1, color = factor(gt40)))

ggplot(phx) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline), col = 'grey') +
  geom_line(aes(x = dt_fmt, y = sim_death_t1)) +
  geom_point(aes(x = dt_fmt, y = sim_death_t1, color = factor(gt40)))

# 
library(gnm)

m1.t1 <- gnm(sim_death_t1 ~ gt40 + 
               ns(dt_fmt,df = df_per_year * nyears), 
             data = phx,
             family = quasipoisson)

coef(m1.t1 )[1:2]
exp(coef(m1.t1 )[1:2])s