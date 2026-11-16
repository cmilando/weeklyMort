library(data.table)
library(ggplot2)
library(dlnm)
library(splines)
library(gnm)
library(patchwork)

df <- read.csv("heat_mortality_asofnov6.csv")
setDT(df)

# get rid of low pop
df <- subset(df, population_estimate > 100000)



# set X_deaths NA to 0
rr <- which(is.na(df$X_deaths))
df$X_deaths[rr] = 0

# make FIPS
df[, fips := sprintf("%05s", occurrence_county_code)]

# make strata
df[, dt_fmt := as.IDate(X_week_end_date)]
df[, stratum := paste0(fips,".", 
                      year(dt_fmt), ".",
                      month(dt_fmt))]

#
# just 2018 and 2019
df <- subset(df, year(dt_fmt) < 2020)

# 
df[, state_fips := substr(fips, 1, 2)]
unique_states <- unique(df$state_fips)
unique_states

# just grab one state
# ma = 25
# california = 06
# az = 04
plotlist <- vector("list", length(unique_states))

for(ii in 1:length(unique_states)) {
  
  # ii = 4
  
  cat(">", unique_states[ii], "\t")
  
  state_df <- subset(df, grepl(paste0("^", unique_states[ii]), fips))
  state_df <- setDT(state_df)
  unique(state_df$fips)
  
  # make some additional variables
  state_df[, log_pop := log(population_estimate)]
  state_df[, yday := yday(dt_fmt)]
  state_df[, is_covid_year := year(dt_fmt) %in% c(2020, 2021)]
  
  # make crossbasis
  qEXP <- quantile(state_df$max_temp_C, probs = c(10, 50, 90)/100)

  cb.1 <- crossbasis(
    x = state_df$max_temp_C,
    argvar = list(fun = 'ns', knots = qEXP),
    arglag = list(fun = 'poly', 2),
    lag = 4, # because its weekly
    group = state_df$fips) # because you are trying this in a 1-stage
  
  # make the interaction terms
  is_covid_yearTRUE = 1*(state_df$is_covid_year)
  covid.int =  cb.1 * is_covid_yearTRUE
  covid.int
  
  # get the spline DF
  nyears = length(unique(year(state_df$dt_fmt)))
  nyears
  df_per_year = 8
  
  # run the CP model
  state_df[,  keep:=sum(X_deaths)>0, by=stratum]
  m0.cp <- gnm(X_deaths ~  cb.1, 
               eliminate = factor(stratum),
               subset = keep,
               data = state_df,
               offset = log_pop,
               family = quasipoisson)
  
  cp.0 = crosspred(cb.1, m0.cp, cen = qEXP[2])
  
  
  # run the baseline time-series model
  m1.x1 <- gnm(X_deaths ~  cb.1 +
                 ns(dt_fmt, df = df_per_year * nyears), 
               data = state_df,
               offset = log_pop,
               family = quasipoisson)
  
  cp.1 = crosspred(cb.1, m1.x1, cen = qEXP[2])
  
  # run the interaction model
  m1.int1 <- update(m1.x1, .~. + is_covid_yearTRUE + covid.int)
  m1.int1
  
  cp.1.int = crosspred(cb.1, m1.int1, cen = qEXP[2])



  
  # Extract predictions
  dt0 <- data.table(
    temp = cp.0$predvar,
    RR = cp.0$allRRfit,
    lower = cp.0$allRRlow,
    upper = cp.0$allRRhigh,
    model = "Time-StratCP"
  )
  
  dt1 <- data.table(
    temp = cp.1$predvar,
    RR = cp.1$allRRfit,
    lower = cp.1$allRRlow,
    upper = cp.1$allRRhigh,
    model = "TS Baseline"
  )
  
  dt2 <- data.table(
    temp = cp.1.int$predvar,
    RR = cp.1.int$allRRfit,
    lower = cp.1.int$allRRlow,
    upper = cp.1.int$allRRhigh,
    model = "TS w Interaction"
  )
  
  dt <- rbind(dt0, dt1, dt2)
  col <- c("forestgreen","royalblue", "tomato3")
  
  # Plot
  pX <- ggplot(dt, aes(x = temp, y = RR, color = model, fill = model)) +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                alpha = 0.2, color = NA) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    scale_y_log10() +
    coord_cartesian(ylim = c(0.5, 10), 
                    xlim = c(min(state_df$max_temp_C),
                             max(state_df$max_temp_C))) +
    labs(
      x = expression(paste("Temperature (" * degree * "C)")),
      y = "RR",
      title = paste0(ii, " - ", state_df$stname[1]),
      color = NULL,
      fill = NULL
    ) +
    theme_classic() +
    theme(
      axis.text.y = element_text(angle = 0),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  plotlist[[ii]] <- pX

}

pall <- wrap_plots(plotlist)

ggsave(pall, filename = "p_all.png", width = 30, height = 30)

# # TEST SIGNIFICANCE OF INTERACTION
# anova(m1.x1, m1.int1, test="Chisq")

# 
