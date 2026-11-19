library(data.table)
library(ggplot2)
library(dlnm)
library(splines)
library(gnm)
library(patchwork)

df <- read.csv("heat_mortality_asofnov6.csv")
setDT(df)

# get rid of low pop
df <- subset(df, population_estimate >= 1e5)

# set X_deaths NA to 0
rr <- which(is.na(df$X_deaths))
df$X_deaths[rr] = 0

# make FIPS
df[, fips := sprintf("%05s", occurrence_county_code)]

# get NCA4 Region
get_region <- function(STATE_COUNTY_FIPS) {
  
  STATE = substr(STATE_COUNTY_FIPS, 1, 2)
  NE <- c('23', '33', '50', '25', '44', '09', '36', '42', 
          '34', '10', '24', '42', '54', '11' )
  MW <- c("27", "55", "26", "19", "29", "17", "18", "39")
  NW <- c("53", "41", "16")
  NGP <- c('30', '56', '38', '46', '31')
  SGP <- c('20', '40', '48')
  SW <- c('06', '32', '49', '04', '08', '35')
  SE <- c('51', '37', '21', '47', '05', '22', '28',
          '01', '13', '45', 12)
  
  if(STATE %in% NE) return('NE')
  if(STATE %in% MW) return('MW')
  if(STATE %in% NW) return('NW')
  if(STATE %in% NGP) return('NGP')
  if(STATE %in% SGP) return('SGP')
  if(STATE %in% SW) return('SW')
  if(STATE %in% SE) return('SE')
  return('Other')
  
}
df$NCAregion <- sapply(df$fips, get_region)

# make strata
df[, dt_fmt := as.IDate(X_week_end_date)]
df[, dt_int := as.integer(dt_fmt)]
df[, stratum := paste0(fips,".", year(dt_fmt), ".", month(dt_fmt))]
df[, state_fips := substr(fips, 1, 2)]

# just 2018 and 2019
# df <- subset(df, year(dt_fmt) < 2020)
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
  n_counties <- length(unique(state_df$fips))
  
  # make some additional variables
  state_df[, log_pop := log(population_estimate)]
  state_df[, yday := yday(dt_fmt)]
  state_df[, is_covid_year := year(dt_fmt) %in% c(2020, 2021)]
  
  # make the exposure var a percentile
  ref_years = 2018:2025
  ref_ecdf <- ecdf(subset(state_df, year(dt_fmt) %in% ref_years)$max_temp_C)
  state_df[, exp_percentile := ref_ecdf(max_temp_C)]
  
  # make crossbasis
  qEXP <- quantile(state_df$exp_percentile, probs = c(10, 75, 90)/100)
  qEXP
  
  cb.1 <- crossbasis(
    x = state_df$exp_percentile,
    argvar = list(fun = 'bs', knots = qEXP, degree = 2),
    arglag = list(fun = 'ns'),
    lag = 3, # because its weekly
    group = state_df$fips) # because you are trying this in a 1-stage
  
  # run the CP model
  
  # 
  
  state_df[,  strata_rat:= max(X_deaths, na.rm = T) / min(X_deaths, na.rm = T), by=stratum]
  state_df
  
  state_df[,  strata_rat:= var(X_deaths, na.rm = T), by=stratum]
  state_df
  
  ## 
  # 4 point regression line
  # standardize the slope coefficient / divided by the variance
  # 
  state_df[,  strata_rat:= coef(lm(X_deaths ~ dt_fmt))[2], by=stratum]
  state_df
  
  # rate rank

  state_df[,  rate_rank:= ecdf(strata_rat)(strata_rat), by=stratum]
  state_df
  
  state_df_agg <- state_df[, .(strata_rat = unique(strata_rat),
                               first_date = min(dt_fmt)), 
                               by = .(stratum, X_year, fips, month, population_estimate)]
  
  state_df_agg$strata_rat_rank = NA
  
  unique_fips <- unique(state_df$fips)
  for(ii in 1:length(unique_fips)) {
    fips_i = unique_fips[ii]
    rr <- which(state_df_agg$fips == fips_i)
    xecdf <- ecdf(state_df_agg$strata_rat[rr])
    state_df_agg$strata_rat_rank[rr] =  xecdf(state_df_agg$strata_rat[rr])
  }
  
  hist(state_df_agg$population_estimate)
  
  kk <- quantile(state_df_agg$population_estimate, 
                 probs = c(0:3)/3)
  kk
  
  state_df_agg[, population_cut := cut(population_estimate, breaks = kk,
                                       include.lowest = T,
                                       labels = c('low', 'med', 'high'))]
  levels(state_df_agg$population_cut)
  table(state_df_agg$population_cut)
  
  quantile(state_df_agg$strata_rat)
  
  ggplot(state_df_agg) +
    geom_tile(aes(x = factor(first_date), 
                  y = reorder(fips, population_estimate) , 
                  fill = strata_rat_rank)) + 
    scale_fill_viridis_c() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 5)) +
    facet_wrap(population_cut ~ .,nrow = 3, scales = 'free_y')
  
  hist(state_df_agg$strata_rat)
  quantile(state_df_agg$strata_rat)
  
  hist(subset(state_df_agg, fips == '06037')$strata_rat)
  
  ggplot(subset(state_df_agg, strata_rat_rank > 0.95)) +
    geom_tile(aes(x = factor(first_date), 
                  y = reorder(fips, population_estimate) , 
                  fill = strata_rat_rank)) + 
    scale_fill_viridis_c(option = 'magma') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0, size = 5)) +
    facet_wrap(population_cut~.,nrow = 3, scales = 'free_y')
  
  # confirming when the COVID wave happened
  xx <- subset(state_df, fips == '06019')
  state_df_agg
  xx <- state_df_agg[xx, on = .(fips, strata_rat)]
  head(xx)
  
  p1 = ggplot(xx) + 
    geom_line(aes(x = dt_fmt, y = X_deaths))
  
  p2 = ggplot(xx) + 
    geom_line(aes(x = dt_fmt, 
                   y = strata_rat_rank)) +
    geom_point(aes(x = dt_fmt, 
                  y = strata_rat_rank,
                  color = cut(strata_rat_rank, 
                              breaks = c(0, 0.95, 1))),
               show.legend = F)
  
  p3 = ggplot(xx) + 
    geom_line(aes(x = dt_fmt, 
                  y = strata_rat)) +
    geom_point(aes(x = dt_fmt, 
                   y = strata_rat))
  
  
  p1 / p2 / p3
  
  state_df[,  
           keep:=(sum(X_deaths) > 0 & strata_rat < 1.2), 
           by=stratum]
  
  m0.cp <- gnm(X_deaths ~  cb.1, 
               eliminate = factor(stratum),
               subset = keep,
               data = state_df,
               offset = log_pop,
               family = quasipoisson)
  
  cp.0 = crosspred(cb.1, m0.cp, cen = qEXP[2])
  
  plot(cp.0, "overall")
  
  
  
}


# 
