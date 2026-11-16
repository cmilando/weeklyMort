
library(data.table)
library(ggplot2)

df <- read.csv("heat_mortality_asofnov6.csv")
setDT(df)

df[, fips := sprintf("%05s", occurrence_county_code)]
df$fips
# 
# head(df)
# 
# 
# # just grab one
phx <- subset(df, grepl("^25", fips))
unique(phx$fips)
# saveRDS(phx, "phx.RDS")

phx <- setDT(phx)
phx[, dt_fmt := as.IDate(X_week_end_date)]

#
ggplot(phx) +
  geom_line(aes(x = dt_fmt, y = X_deaths, color = fips))

ggplot(phx) +
  geom_line(aes(x = dt_fmt, y = mean_temp_C, color = fips))

ggplot(phx) +
  geom_histogram(aes(x = mean_temp_C, fill = fips))

quantile(phx$mean_temp_C, probs = 0.95)

# make strata of county-year-month
phx[, strata := paste0(fips,".", 
                       year(dt_fmt), ".",
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
hist(phx$mean_temp_C)
qX = quantile(phx$mean_temp_C, probs = .90)
phx[, grQX := 1*(mean_temp_C >= qX)]
phx
table(phx$grQX)

beta_temp = log(2)

phx[, sim_death := sim_death_baseline * exp(grQX * beta_temp)]
phx

ggplot(subset(phx, year(dt_fmt) == 2018)) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline), col = 'grey') +
  geom_line(aes(x = dt_fmt, y = sim_death)) +
  geom_point(aes(x = dt_fmt, y = sim_death, color = factor(grQX)))

# 
library(gnm)

m1 <- gnm(sim_death ~ grQX, 
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
  geom_point(aes(x = dt_fmt, y = sim_death, color = factor(grQX))) + 
  facet_wrap(~strata)


#
library(splines)
phx[, yr := year(dt_fmt)]
nyears = length(unique(phx$yr))
nyears
df_per_year = 8

m1.mod <- gnm(sim_death ~ grQX + 
                ns(dt_fmt,df = df_per_year * nyears), 
              data = phx,
              family = quasipoisson)

coef(m1.mod )[1:2]
exp(coef(m1.mod )[1:2])

# =============================================================================
# /////////////////////////////////////////////////////////////////////////////
# TRIAL #1: COVID EFFECT IS ADDITIVE
#
# log(death) = log(baseline) + beta_temp * grQX + beta_covid * is_covid_yr
#
# /////////////////////////////////////////////////////////////////////////////
# =============================================================================

beta_covid_add = log(2)

phx[, is_covid_year := year(dt_fmt) %in% c(2020, 2021)]

phx[, sim_death_t1 := sim_death_baseline * 
      exp(grQX * beta_temp + is_covid_year * beta_covid_add)
]

# plot
ggplot(subset(phx, year(dt_fmt) == 2020)) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline), col = 'grey') +
  geom_line(aes(x = dt_fmt, y = sim_death_t1)) +
  geom_point(aes(x = dt_fmt, y = sim_death_t1, color = factor(grQX)))

ggplot(phx) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline), col = 'grey') +
  geom_line(aes(x = dt_fmt, y = sim_death_t1)) +
  geom_point(aes(x = dt_fmt, y = sim_death_t1, color = factor(grQX)))

# 
library(gnm)

m1.t1 <- gnm(sim_death_t1 ~ grQX + 
            ns(dt_fmt,df = df_per_year * nyears), 
          data = phx,
          family = quasipoisson)

coef(m1.t1 )[1:2]
exp(coef(m1.t1 )[1:2])
 
# ----
# RESULT IS
# if covid effect is additive ONLY
# the time-series method adjusts properly without needing a COVID factor
# wherease the time-sratified by month-county-year
# does not
#
# -----


# =============================================================================
# /////////////////////////////////////////////////////////////////////////////
# TRIAL # 2: COVID EFFECT IS MULTIPLICATIVE
#
# log(death) = log(baseline) + 
#              beta_temp * grQX * (1 - is_covid_yr) + 
#              beta_covid * grQX * (is_covid_yr)
#
# /////////////////////////////////////////////////////////////////////////////
# =============================================================================

beta_covid_mult = log(3)

phx[, sim_death_t2 := sim_death_baseline * 
      exp(grQX * beta_temp +
            is_covid_year * (beta_covid_add + grQX * beta_covid_mult))
]

# plot
ggplot(subset(phx, year(dt_fmt) == 2020)) +
  geom_line(aes(x = dt_fmt, y = sim_death_baseline), col = 'grey') +
  geom_line(aes(x = dt_fmt, y = sim_death_t2)) +
  geom_point(aes(x = dt_fmt, y = sim_death_t2, color = factor(grQX)))

#
m1.t2 <- gnm(sim_death_t2 ~ grQX * is_covid_year + 
               ns(dt_fmt,df = df_per_year * nyears), 
             data = phx,
             family = quasipoisson)

coef(m1.t2 )[1:4]
exp(coef(m1.t2 )[1:4])

#
m1.t2b <- gnm(sim_death_t2 ~ grQX  + 
               ns(dt_fmt,df = df_per_year * nyears), 
             data = phx,
             family = quasipoisson)

coef(m1.t2b )[1:2]
exp(coef(m1.t2b )[1:2])

# ----
# RESULT IS
# if covid effect is additive AND multiplicative
# you need to add an interaction term with is_covid_year
#
# -----

# =============================================================================
# /////////////////////////////////////////////////////////////////////////////
# TRIAL # 3: THERE IS HOPE
#
# log(death) = log(baseline) + 
#              beta_temp * grQX * (1 - is_covid_yr) + 
#              beta_covid * grQX * (is_covid_yr)
#
# /////////////////////////////////////////////////////////////////////////////
# =============================================================================


m1.x1 <- gnm(X_deaths ~  grQX * is_covid_year +
                ns(dt_fmt,df = df_per_year * nyears), 
              data = phx,
              family = quasipoisson)

exp(confint(m1.x1))[1:4,]


m1.x1b <- gnm(X_deaths ~  grQX  +
               ns(dt_fmt,df = df_per_year * nyears), 
             data = phx,
             family = quasipoisson)

exp(confint(m1.x1b))[1:4,]

# =============================================================================
# /////////////////////////////////////////////////////////////////////////////
# TRIAL # 4: DLNM
#
# https://github.com/gasparrini/CTS-smallarea/blob/ec37a627611547e1ea81ce2df329dd072a63d349/04.intmod.R#L14
#
# log(death) = log(baseline) + 
#              beta_temp * grQX * (1 - is_covid_yr) + 
#              beta_covid * grQX * (is_covid_yr)
#
# /////////////////////////////////////////////////////////////////////////////
# =============================================================================

library(dlnm)

cb.1 <- crossbasis(phx$mean_temp_C,
                   argvar = list(fun = 'ns', knots = c(10, 20)),
                   arglag = list(fun = 'poly', degree = 4),
                   lag = 21,
                   group = phx$fips)

cb.1

#
# cb.1 
# 0.213397150 
# is_covid_yearTRUE 
# -0.093127869 
# cb.1:is_covid_yearTRUE 
# -0.241527604

# i think you just have to write out the vars by hand in order to use 
# them in crosspred
# -- is_covid_year_TRUE
# -- cb1 * is_covid_year_TRUE
is_covid_yearTRUE = 1*(phx$is_covid_year)
covid.int =  cb.1 * is_covid_yearTRUE
covid.int

m1.x1 <- gnm(X_deaths ~  cb.1 +
               ns(dt_fmt,df = df_per_year * nyears), 
             data = phx,
             family = quasipoisson)

m1.int1 <- update(m1.x1, .~. + covid.int)
m1.int1

# TEST SIGNIFICANCE OF INTERACTION
anova(m1.x1, m1.int1, test="Chisq")

coef(m1.x1)[1:5]

any(is.na(coef(m1.int1)))

cp.1 = crosspred(cb.1, m1.x1)
plot(cp.1, "overall")


cp.1 = crosspred(cb.1, m1.int1)
plot(cp.1, "overall")

exp(confint(m1.x1))[1:4,]


m1.x1b <- gnm(X_deaths ~  grQX  +
                ns(dt_fmt,df = df_per_year * nyears), 
              data = phx,
              family = quasipoisson)

exp(confint(m1.x1b))[1:4,]
