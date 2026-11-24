library(data.table); library(ggplot2); library(cityClimateHealth)

# Anomaly detection:
# D1 absolute: |within-month slope of deaths|
# D2 relative: |slope| / stratum mean
# D3 neighbour: |log(stratum mean / mean of the previous and next month)|
# "Anomalous" = above the county's own 95th percentile of each metric in reference (non-COVID) years.

library(data.table); library(dplyr)
df <- read.csv("heat_mortality_asofnov6.csv") %>%
  select(!X)
setDT(df)

# get rid of low pop
df <- subset(df, population_estimate >= 1e5)

# set X_deaths NA to 0
rr <- which(is.na(df$X_deaths))
df$X_deaths[rr] = 0

# make FIPS
df$occurrence_county_code <- sprintf("%05d",df$occurrence_county_code)
df$state_fips <- sprintf("%02d", df$state_fips)

#correcting date and outcome format
df[, X_week_end_date := as.Date(X_week_end_date)]
df[, X_deaths := as.integer(X_deaths)]
head(df)

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
df$NCAregion <- sapply(df$state_fips, get_region)

head(df)

# make temp percentiles by county
pct_by <- "occurrence_county_code"
pct_ref_years <- 2018:2026

df[, pct_year := year(as.IDate(X_week_end_date))]
df[, exp_percentile := {
  ref <- max_temp_C[pct_year %in% pct_ref_years & !is.na(max_temp_C)]
  out <- 100 * ecdf(ref)(max_temp_C)
  out[is.na(max_temp_C)] <- NA_real_
  out
}, by = pct_grp]
cat("Exposure percentile by", pct_by, "| reference years", paste(range(pct_ref_years), collapse = "-"),
    "| mean in reference years:", round(df[pct_year %in% pct_ref_years, mean(exp_percentile, na.rm = TRUE)], 1), "\n")

flag_q <- 0.95 #pctile for flagging
min_r2 <- 0.5 #for slope flagging
ref_years <- c(2018, 2019, 2023, 2024, 2025, 2026)
defs <- c("D1 absolute", "D2 relative", "D3 neighbour")
knots <- c(0.10, 0.75, 0.90)

# Exposure matrix and outcome columns

exposure_columns_grp <- list(date = "X_week_end_date", exposure = "exp_percentile",
                         geo_unit = "occurrence_county_code", geo_unit_grp = "stname")
outcome_columns_grp  <- list(date = "X_week_end_date", outcome = "X_deaths",
                         geo_unit = "occurrence_county_code", geo_unit_grp = "stname")

exp_data_grp <- as.data.table(df)[, .(X_week_end_date, exp_percentile, occurrence_county_code, stname)]

out_data_grp <- as.data.table(df)[, .(X_week_end_date, X_deaths, occurrence_county_code, stname)]

exposure_matrix_grp <- make_exposure_matrix(
  exp_data_grp, exposure_columns_grp, dt_by = "week",
  time_subset = list(month = 1:12, year = c(2018:2026)),
  grp_level = TRUE, maxlag = 4, keep_unit_exposures = TRUE)


# Flag strata

wk <- as.data.table(df)[year(as.IDate(X_week_end_date)) %in% 2018:2026,
                              .(geo = occurrence_county_code, date = as.IDate(X_week_end_date),
                                deaths = as.numeric(X_deaths))]
wk <- wk[, .(deaths = sum(deaths)), by = .(geo, date)]
if (anyNA(wk$deaths)) wk[is.na(deaths), deaths := 5]       # suppressed weeks as 5 for the metrics only
wk[, `:=`(year = year(date), month = month(date))]
setorder(wk, geo, date)
wk[, t := seq_len(.N) - 1, by = .(geo, year, month)]

st <- wk[, {
  n  <- .N; m <- mean(deaths)
  sl <- if (n >= 3 && var(t) > 0) cov(t, deaths) / var(t) else NA_real_
  r2 <- if (n >= 3 && var(deaths) > 0) cor(t, deaths)^2 else 0
  .(mean_d = m, slope = sl, r2 = r2)
}, by = .(geo, year, month)]
st[, ym := year * 12L + month]
st <- merge(st, st[, .(geo, ym = ym + 1L, prev = mean_d)], by = c("geo", "ym"), all.x = TRUE)
st <- merge(st, st[, .(geo, ym = ym - 1L, nxt = mean_d)],  by = c("geo", "ym"), all.x = TRUE)
st[, nb := rowMeans(cbind(prev, nxt), na.rm = TRUE)]
st[, `:=`(m1 = abs(slope),
          m2 = fifelse(mean_d > 0, abs(slope) / mean_d, NA_real_),
          m3 = fifelse(mean_d > 0 & nb > 0, abs(log(mean_d / nb)), NA_real_))]
thr <- st[year %in% ref_years, .(t1 = quantile(m1, flag_q, na.rm = TRUE),
                                 t2 = quantile(m2[r2 >= min_r2], flag_q, na.rm = TRUE),
                                 t3 = quantile(m3, flag_q, na.rm = TRUE)), by = geo]
st <- merge(st, thr, by = "geo")
st[, `:=`(`D1 absolute`  = !is.na(m1) & m1 > t1,
          `D2 relative`  = !is.na(m2) & r2 >= min_r2 & !is.na(t2) & m2 > t2,
          `D3 neighbour` = !is.na(m3) & m3 > t3)]

cat("\nShare of county-months flagged (%), by year\n")
print(dcast(melt(st, id.vars = "year", measure.vars = defs, variable.name = "def")[
  , .(pct = round(100 * mean(value), 1)), by = .(year, def)], year ~ def, value.var = "pct"))

# remake outcome data with each anomaly def weeks removed

od <- as.data.table(out_data_grp)
od[, `:=`(year = year(X_week_end_date), month = month(X_week_end_date))]
od <- merge(od, st[, .(occurrence_county_code = geo, year, month,
                       D1 = `D1 absolute`, D2 = `D2 relative`, D3 = `D3 neighbour`)],
            by = c("occurrence_county_code", "year", "month"), all.x = TRUE)
od[is.na(D1), D1 := FALSE]; od[is.na(D2), D2 := FALSE]; od[is.na(D3), D3 := FALSE]

out_data_D1 <- od[D1 == FALSE, .(X_week_end_date, X_deaths, occurrence_county_code, stname)]
out_data_D2 <- od[D2 == FALSE, .(X_week_end_date, X_deaths, occurrence_county_code, stname)]
out_data_D3 <- od[D3 == FALSE, .(X_week_end_date, X_deaths, occurrence_county_code, stname)]

# make the outcome tables
outcomes_tbl_all <- make_outcome_table(
  out_data_grp, outcome_columns_grp, dt_by = "week",
  time_subset = list(month = 1:12, year = c(2018:2026)),
  grp_level = TRUE, keep_unit_outcomes = TRUE)
outcomes_tbl_D1 <- make_outcome_table(
  out_data_D1, outcome_columns_grp, dt_by = "week",
  time_subset = list(month = 1:12, year = c(2018:2026)),
  grp_level = TRUE, keep_unit_outcomes = TRUE)
outcomes_tbl_D2 <- make_outcome_table(
  out_data_D2, outcome_columns_grp, dt_by = "week",
  time_subset = list(month = 1:12, year = c(2018:2026)),
  grp_level = TRUE, keep_unit_outcomes = TRUE)
outcomes_tbl_D3 <- make_outcome_table(
  out_data_D3, outcome_columns_grp, dt_by = "week",
  time_subset = list(month = 1:12, year = c(2018:2026)),
  grp_level = TRUE, keep_unit_outcomes = TRUE)

#run the models
model_all <- condPois_2stage(exposure_matrix_grp, outcomes_tbl_all,
                             argvar = list(fun = "ns", pct = knots, degree = 2),
                             arglag = list(fun = "poly"), maxlag = 3, verbose = 1, truncate = 0)
                             # arglag = list(fun = "ns"), maxlag = 3, verbose = 1, truncate = 0)
model_D1 <- condPois_2stage(exposure_matrix_grp, outcomes_tbl_D1,
                            argvar = list(fun = "ns", pct = knots, degree = 2),
                            arglag = list(fun = "poly"), maxlag = 3, verbose = 1, truncate = 0)
model_D2 <- condPois_2stage(exposure_matrix_grp, outcomes_tbl_D2,
                            argvar = list(fun = "ns", pct = knots, degree = 2),
                            arglag = list(fun = "poly"), maxlag = 3, verbose = 1, truncate = 0)
model_D3 <- condPois_2stage(exposure_matrix_grp, outcomes_tbl_D3,
                            argvar = list(fun = "ns", pct = knots, degree = 2),
                            arglag = list(fun = "poly"), maxlag = 3, verbose = 1, truncate = 0)

models <- list("All weeks" = model_all, "Absolute" = model_D1,
               "Relative" = model_D2, "Neighbor" = model_D3)


# Plot ERFs
cols_m <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"), names(models))
theme_erf <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "#e6e5e0", linewidth = 0.3),
        axis.text = element_text(colour = "#52514e"), axis.title = element_text(colour = "#52514e"),
        strip.text = element_text(face = "bold", hjust = 0), legend.position = "top",
        legend.title = element_blank(), plot.caption = element_text(colour = "#52514e", hjust = 0),
        plot.background = element_rect(fill = "#fcfcfb", colour = NA))

# pooled curve
pooled <- rbindlist(lapply(names(models), \(nm)
  data.table(model = nm, as.data.table(models[[nm]]$`_`$grp_plt$data))))
pooled[, model := factor(model, levels = names(models))]
p_pooled <- ggplot(pooled, aes(x, RR, colour = model, fill = model)) +
  geom_hline(yintercept = 1, colour = "#898781", linewidth = 0.3) +
  geom_ribbon(aes(ymin = RRlb, ymax = RRub), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.8) +
  scale_colour_manual(values = cols_m) + scale_fill_manual(values = cols_m) +
  scale_y_continuous(transform = "log", breaks = c(0.8, 0.9, 1, 1.1, 1.25, 1.5, 2)) +
  labs(title = "", x = "Weekly temperature percentile",
       y = "RR") +
  theme_erf
print(p_pooled)


# state curves
st_rr <- rbindlist(lapply(names(models), \(nm) data.table(model = nm, as.data.table(getRR(models[[nm]])))))
xcol <- exposure_columns_grp$exposure
setnames(st_rr, xcol, "x", skip_absent = TRUE)
gcol <- intersect(c("stname", outcome_columns_grp$geo_unit_grp), names(st_rr))[1]
if ("stage" %in% names(st_rr)) st_rr <- st_rr[stage == "stage2"]
st_rr[, model := factor(model, levels = names(models))]
ylim_st <- c(max(0.7, min(st_rr$RR)), min(2, max(st_rr$RR)))
p_states <- ggplot(st_rr, aes(x, RR, colour = model)) +
  geom_hline(yintercept = 1, colour = "#898781", linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  scale_colour_manual(values = cols_m) +
  scale_y_continuous(transform = "log", breaks = c(0.8, 0.9, 1, 1.03, 1.05, 1.07, 1.1, 1.25, 1.5)) +
  # coord_cartesian(ylim = ylim_st) +
  facet_wrap(as.formula(paste("~", gcol)), ncol = 7, scales = "free_y") +
  labs(title = "State ERFs", x = "Weekly temperature percentile",
       y = "RR") +
  theme_erf + theme(axis.text = element_text(size = 7), strip.text = element_text(size = 8))
print(p_states)

#mmts across different anom definitions
mmt <- rbindlist(lapply(names(models), \(nm) data.table(model = nm,
  unit = names(models[[nm]]$`_`$out), mmt = sapply(models[[nm]]$`_`$out, `[[`, "cen"))))
cat("\nState MMT percentile under each model\n")
print(dcast(mmt, unit ~ model, value.var = "mmt"))

