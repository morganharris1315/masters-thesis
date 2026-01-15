# -------------------------------------------------------------------------
# 03a-blockmaxmia-modeldata.R
# -------------------------------------------------------------------------
# Dec 2025/ Jan 2025
# 
# -------------------------------------------------------------------------

#install.packages('ismev')
#install.packages('fExtremes')
#install.packages('extRemes')

library(ismev)
library(fExtremes)

getwd()
model_data_dir <- glue("{base_raw_dir}/model_data")
milford_CD <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/wah_ens_daily_rain_MilfordSound_CD.csv")

milford_CD$RX1day <-apply(milford_CD, 1, max, na.rm=TRUE)
milford_CD <- milford_CD %>%
  mutate(Year = row_number())
print(milford_CD)

ggplot(data = milford_CD, aes(x = Year, y = RX1day)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Milford Sound Current Day",
    x = "Year",
    y = "Rainfall (mm) on RX1day")


# Defining thresholds -----------------------------------------------------
minRX1day <-min(milford_CD$RX1day)
print(minRX1day)
meanRX1day <-mean(milford_CD$RX1day)
print(meanRX1day)
maxRX1day <-max(milford_CD$RX1day)
print(maxRX1day)
midRX1day<- (minRX1day + meanRX1day) / 2
print(midRX1day)


# Plotting thresholds -----------------------------------------------------
ggplot(data = milford_CD, aes(x = Year, y = RX1day)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Milford Sound Current Day",
    x = "Year",
    y = "Rainfall (mm) on RX1day") +
  geom_hline(yintercept = minRX1day, color = "red", linetype = "dashed", size = 1)+
  geom_hline(yintercept = meanRX1day, color = "red", linetype = "dashed", size = 1)+
  geom_hline(yintercept = midRX1day, color = "red", linetype = "dashed", size = 1)+
  geom_hline(yintercept = maxRX1day, color = "red", linetype = "dashed", size = 1)


# Proportion of days above threshold --------------------------------------
# select only daily rainfall columns
milford_CD_daily <- milford_CD %>%
  select(-RX1day, -Year)

milford_CD_daily_vec <- as.numeric(as.matrix(milford_CD_daily))

n_daily_exceed_min  <- sum(milford_CD_daily_vec > minRX1day,  na.rm = TRUE)
n_daily_exceed_mean <- sum(milford_CD_daily_vec > meanRX1day, na.rm = TRUE)
n_daily_exceed_mid  <- sum(milford_CD_daily_vec > midRX1day,  na.rm = TRUE)
n_daily_exceed_max  <- sum(milford_CD_daily_vec > maxRX1day,  na.rm = TRUE)

# total number of valid daily observations
n_days <- sum(!is.na(milford_CD_daily_vec))

prop_exceed_min  <- n_daily_exceed_min  / n_days
prop_exceed_mean <- n_daily_exceed_mean / n_days
prop_exceed_mid  <- n_daily_exceed_mid  / n_days
prop_exceed_max  <- n_daily_exceed_max  / n_days

print(prop_exceed_min)
print(prop_exceed_mid)
print(prop_exceed_mean)
print(prop_exceed_max)

prop_exceed_onein360 <- 1 / 360
onein360_quantile <- 1 - prop_exceed_onein360

threshold_onein360 <- quantile(milford_CD_daily_vec, probs = onein360_quantile, na.rm = TRUE, type = 7)
print(threshold_onein360)

# common y-axis limits for comparability -------------------------------
y_max_common <- max(
  milford_CD$RX1day,
  milford_FP$RX1day,
  na.rm = TRUE
)
y_limits_common <- c(0, y_max_common * 1.05)  # 5% headroom

# re-plotting thresholds -----------------------------------------------------
# legend labels for current day -------------------------------------------
labels_CD <- c(
  "Min RX1day" = paste0(
    "Min RX1day (",
    round(minRX1day, 1), " mm, ",
    round(prop_exceed_min * 100, 2), "%)"
  ),
  "Mid RX1day" = paste0(
    "Mid RX1day (",
    round(midRX1day, 1), " mm, ",
    round(prop_exceed_mid * 100, 2), "%)"
  ),
  "Mean RX1day" = paste0(
    "Mean RX1day (",
    round(meanRX1day, 1), " mm, ",
    round(prop_exceed_mean * 100, 2), "%)"
  ),
  "1 in 360 day threshold" = paste0(
    "1 in 360 day (",
    round(threshold_onein360, 1), " mm, ",
    round(prop_exceed_onein360 * 100, 2), "%)"
  ),
  "Max RX1day" = paste0(
    "Max RX1day (",
    round(maxRX1day, 1), " mm, ",
    round(prop_exceed_max * 100, 2), "%)"
  )
)

p_CD <- ggplot(data = milford_CD, aes(x = Year, y = RX1day)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Milford Sound Current Day",
    x = "Year",
    y = "Rainfall (mm) on RX1day",
    colour = "Threshold"
  ) +
  geom_hline(aes(yintercept = minRX1day, colour = "Min RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = meanRX1day, colour = "Mean RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = midRX1day, colour = "Mid RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = maxRX1day, colour = "Max RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = threshold_onein360, colour = "1 in 360 day threshold"),
             linetype = "solid", size = 1.1) +
scale_y_continuous(limits = y_limits_common) +
scale_colour_manual(
  breaks = c(
    "Max RX1day",
    "Mean RX1day",
    "1 in 360 day threshold",
    "Mid RX1day",
    "Min RX1day"
  ),
  labels = labels_CD,
  values = c(
    "Min RX1day"             = "#56B4E9",  
    "Mid RX1day"             = "#0072B2", 
    "1 in 360 day threshold" = "#E69F00",  
    "Mean RX1day"            = "#009E73",  
    "Max RX1day"             = "#D55E00"   
  )
)


# Future projections ------------------------------------------------------
milford_FP <- read.csv("C:/Users/morga/OneDrive - The University of Waikato/Masters Thesis/Thesis/Historic Compound Events/model_data/wah_ens_daily_rain_MilfordSound_3deg.csv")

# Calculate annual maxima (RX1day) ----------------------------------------
milford_FP$RX1day <- apply(milford_FP, 1, max, na.rm=TRUE)
milford_FP <- milford_FP %>%
  mutate(Year = row_number())
print(milford_FP)

# Plot RX1day for future projections --------------------------------------
ggplot(data = milford_FP, aes(x = Year, y = RX1day)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Milford Sound Future Projections",
    x = "Year",
    y = "Rainfall (mm) on RX1day"
  )

minRX1day_FP <- min(milford_FP$RX1day)
meanRX1day_FP <- mean(milford_FP$RX1day)
maxRX1day_FP <- max(milford_FP$RX1day)
midRX1day_FP <- (minRX1day_FP + meanRX1day_FP) / 2

print(minRX1day_FP)
print(meanRX1day_FP)
print(maxRX1day_FP)
print(midRX1day_FP)

milford_FP_daily <- milford_FP %>%
  select(-RX1day, -Year)

milford_FP_daily_vec <- as.numeric(as.matrix(milford_FP_daily))

n_days_FP <- sum(!is.na(milford_FP_daily_vec))

prop_exceed_min_FP  <- sum(milford_FP_daily_vec > minRX1day_FP,  na.rm = TRUE) / n_days_FP
prop_exceed_mid_FP  <- sum(milford_FP_daily_vec > midRX1day_FP,  na.rm = TRUE) / n_days_FP
prop_exceed_mean_FP <- sum(milford_FP_daily_vec > meanRX1day_FP, na.rm = TRUE) / n_days_FP
prop_exceed_max_FP  <- sum(milford_FP_daily_vec > maxRX1day_FP,  na.rm = TRUE) / n_days_FP

print(prop_exceed_min_FP)
print(prop_exceed_mid_FP)
print(prop_exceed_mean_FP)
print(prop_exceed_max_FP)

prop_exceed_onein360_FP <- 1 / 360
onein360_quantile_FP <- 1 - prop_exceed_onein360_FP
threshold_onein360_FP <- quantile(milford_FP_daily_vec, probs = onein360_quantile_FP, na.rm = TRUE, type = 7)
print(threshold_onein360_FP)

p_FP <- ggplot(data = milford_FP, aes(x = Year, y = RX1day)) +
  geom_line() +
  theme_minimal() +
  labs(
    title = "Milford Sound Future Projections",
    x = "Year",
    y = "Rainfall (mm) on RX1day",
    colour = "Threshold"
  ) +
  geom_hline(aes(yintercept = minRX1day_FP, colour = "Min RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = meanRX1day_FP, colour = "Mean RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = midRX1day_FP, colour = "Mid RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = maxRX1day_FP, colour = "Max RX1day"),
             linetype = "dashed", size = 1.1) +
  geom_hline(aes(yintercept = threshold_onein360_FP, colour = "1 in 360 day threshold"),
             linetype = "solid", size = 1.1) +
  scale_y_continuous(limits = y_limits_common) +
  scale_colour_manual(
  breaks = c(
    "Max RX1day",
    "Mean RX1day",
    "1 in 360 day threshold",
    "Mid RX1day",
    "Min RX1day"
  ),
  labels = labels_FP,
  values = c(
    "Min RX1day"             = "#56B4E9",  
    "Mid RX1day"             = "#0072B2", 
    "1 in 360 day threshold" = "#E69F00",  
    "Mean RX1day"            = "#009E73",  
    "Max RX1day"             = "#D55E00"   
  )
)

p_CD / p_FP +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


# histogram version one ---------------------------------------------------------------

# Helper function 
count_exceedances_per_year <- function(df, threshold) {
  df_daily <- df %>% select(-RX1day, -Year)
  apply(df_daily, 1, function(x) sum(x > threshold, na.rm = TRUE))
}

# Calculate exceedances
cd_exceed_mid  <- count_exceedances_per_year(milford_CD, midRX1day)
cd_exceed_360  <- count_exceedances_per_year(milford_CD, threshold_onein360)
cd_exceed_mean <- count_exceedances_per_year(milford_CD, meanRX1day)

fp_exceed_mid  <- count_exceedances_per_year(milford_FP, midRX1day_FP)
fp_exceed_360  <- count_exceedances_per_year(milford_FP, threshold_onein360_FP)
fp_exceed_mean <- count_exceedances_per_year(milford_FP, meanRX1day_FP)

# Data frames 
df_cd_mid  <- data.frame(days = cd_exceed_mid)
df_cd_360  <- data.frame(days = cd_exceed_360)
df_cd_mean <- data.frame(days = cd_exceed_mean)

df_fp_mid  <- data.frame(days = fp_exceed_mid)
df_fp_360  <- data.frame(days = fp_exceed_360)
df_fp_mean <- data.frame(days = fp_exceed_mean)

# Common axis limits
x_lim <- c(-0.5, 12.5)
y_lim <- c(0, 0.6)

# Current Day histograms 
p_cd_mid <- ggplot(df_cd_mid, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#0072B2"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim,
                     expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> Mid RX1day",
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )

p_cd_360 <- ggplot(df_cd_360, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#E69F00"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim,
                     expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> 1-in-360 day threshold",
    x = "Number of exceedance days per year",
    y = NULL
  )

p_cd_mean <- ggplot(df_cd_mean, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#009E73"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim,
                     expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = "> Mean RX1day",
    x = "Number of exceedance days per year",
    y = NULL
  )

# Future Projections histograms 
p_fp_mid <- ggplot(df_fp_mid, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#0072B2"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim,
                     expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = NULL,
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )

p_fp_360 <- ggplot(df_fp_360, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#E69F00"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim,
                     expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = NULL,
    x = "Number of exceedance days per year",
    y = NULL
  )

p_fp_mean <- ggplot(df_fp_mean, aes(x = days)) +
  geom_histogram(
    aes(y = after_stat(count / sum(count))),
    binwidth = 1, boundary = 0, closed = "left",
    fill = "#009E73"
  ) +
  scale_x_continuous(limits = x_lim, breaks = 0:12) +
  scale_y_continuous(limits = y_lim,
                     expand = expansion(mult = c(0, 0.02))) +
  theme_minimal() +
  labs(
    title = NULL,
    x = "Number of exceedance days per year",
    y = NULL
  )

# Final layout 
(p_cd_mid | p_cd_360 | p_cd_mean) /
  (p_fp_mid | p_fp_360 | p_fp_mean)


# histogram version two ---------------------------------------------------------------

count_exceedances_per_year <- function(df, threshold) {
  df_daily <- df %>% select(-RX1day, -Year)
  apply(df_daily, 1, function(x) sum(x > threshold, na.rm = TRUE))
}

# Build tidy exceedance dataset 
hist_df <- bind_rows(
  # Current Day
  data.frame(days = count_exceedances_per_year(milford_CD, midRX1day),
             Threshold = "> Mid RX1day",
             Period = "Current Day"),
  data.frame(days = count_exceedances_per_year(milford_CD, threshold_onein360),
             Threshold = "> 1-in-360 day threshold",
             Period = "Current Day"),
  data.frame(days = count_exceedances_per_year(milford_CD, meanRX1day),
             Threshold = "> Mean RX1day",
             Period = "Current Day"),
  
  # Future Projections
  data.frame(days = count_exceedances_per_year(milford_FP, midRX1day_FP),
             Threshold = "> Mid RX1day",
             Period = "Future Projection"),
  data.frame(days = count_exceedances_per_year(milford_FP, threshold_onein360_FP),
             Threshold = "> 1-in-360 day threshold",
             Period = "Future Projection"),
  data.frame(days = count_exceedances_per_year(milford_FP, meanRX1day_FP),
             Threshold = "> Mean RX1day",
             Period = "Future Projection")
)

# Factor ordering 
hist_df$Threshold <- factor(
  hist_df$Threshold,
  levels = c("> Mid RX1day",
             "> 1-in-360 day threshold",
             "> Mean RX1day")
)

# Convert counts to proportions per panel 
hist_df_prop <- hist_df %>%
  count(Period, Threshold, days) %>%
  group_by(Period, Threshold) %>%
  mutate(prop_years = n / sum(n)) %>%
  ungroup()

# Final plot 
ggplot(hist_df_prop, aes(x = days, y = prop_years, fill = Threshold)) +
  geom_col(width = 0.9) +
  facet_grid(
    Period ~ Threshold,
    switch = "y"
  ) +
  scale_x_continuous(
    limits = c(-0.5, 12.5),
    breaks = 0:12
  ) +
  scale_y_continuous(
    limits = c(0, 0.6),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    values = c(
      "> Mid RX1day" = "#0072B2",
      "> 1-in-360 day threshold" = "#E69F00",
      "> Mean RX1day" = "#009E73"
    ),
    guide = "none"
  ) +
  theme_minimal() +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90),
    panel.spacing = unit(1, "lines")
  ) +
  labs(
    x = "Number of exceedance days per year",
    y = "Proportion of years"
  )


# Year comparisions --------------------------------------------------------


# MidRX1day ---------------------------------------------------------------

# Identifying muted and high-exceedance years (Current Day, mid RX1day)
cd_exceed_mid <- count_exceedances_per_year(milford_CD, midRX1day)

exceedance_mid_df_cd <- milford_CD %>%
  select(Year) %>%
  mutate(
    exceedances = cd_exceed_mid,
    )

# Inspecting years - Current Day 
exceedance_mid_df_cd %>%
  filter(exceedances==0)
# There are 305 years where the mid threshold was not exceed. 
# Example year 18

exceedance_mid_df_cd %>%
  arrange(desc(exceedances))
# Two years have 10 exceedances (the highest number of exceedances above the mid threshold. 
# These are years 1464 and 1617

# Identifying muted and high-exceedance years (Future Projections, mid RX1day)
fp_exceed_mid <- count_exceedances_per_year(milford_FP, midRX1day_FP)

exceedance_mid_df_fp <- milford_FP %>%
  select(Year) %>%
  mutate(
    exceedances = fp_exceed_mid,
  )

# Inspecting years - Future Projections
exceedance_mid_df_fp %>%
  filter(exceedances==0)
#there are 226 years where the mid threshold was not exceed. 
#example year 12. 

exceedance_mid_df_fp %>%
  arrange(desc(exceedances))
# 12 times is the largest number of exceedances this happens in the year 1633

#Filtering years
#Muted curent day
milford_CD_year18_muted <- milford_CD %>%
  filter(Year == 18)
milford_CD_year18_muted

milford_CD_year18_muted <- milford_CD_year18_muted %>%
  select(starts_with("X")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "day",
    values_to = "rainfall_mm"
  ) %>%
  mutate(
    day = as.integer(sub("X", "", day))  
  )
print(milford_CD_year18_muted)

#High exccedance current day
milford_CD_year1464_highexceed <- milford_CD %>%
  filter(Year == 1464)
milford_CD_year1464_highexceed

milford_CD_year1464_highexceed <- milford_CD_year1464_highexceed %>%
  select(starts_with("X")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "day",
    values_to = "rainfall_mm"
  ) %>%
  mutate(
    day = as.integer(sub("X", "", day))  
  )
print(milford_CD_year1464_highexceed)

#Plotting time series
plotmilford_CD_year1464_highexceed <- ggplot(milford_CD_year1464_highexceed, aes(x=day, y=rainfall_mm)) +
  geom_line() +
  labs(
    title = "High Exceedance Current Day (Year 1464)",
    x = "Day of Year",
    y = "Daily Rainfall (mm)"
  )+ 
  scale_y_continuous(limits = c(0, 260), breaks = seq(50, 250, by = 50))+
  theme_minimal() + 
  geom_hline(aes(yintercept = midRX1day), col="#0072B2", linetype = "dashed", size = 1.1) 

plotmilford_CD_year18_muted <- ggplot(milford_CD_year18_muted, aes(x=day, y=rainfall_mm)) +
  geom_line() + 
  labs(
    title = "Muted Current Day (Year 18)",
    x = "Day of Year",
    y = ""
  )+ 
  scale_y_continuous(limits = c(0, 260), breaks = seq(50, 250, by = 50))+
  theme_minimal() + 
  geom_hline(aes(yintercept = midRX1day), col="#0072B2", linetype = "dashed", size = 1.1) 

plotmilford_CD_year1464_highexceed + plotmilford_CD_year18_muted

#Filtering years
n_days <- ncol(milford_FP) - 2  # Number of columns to rename (everything except last 2)
new_names <- c(as.character(1:n_days), "RX1day", "Year")
colnames(milford_FP) <- new_names
head(milford_FP)

#Muted future projection
milford_FP_year12_muted <- milford_FP %>%
  filter(Year == 12) %>%
  pivot_longer(
    cols = `1`:(ncol(milford_FP)-2), 
    names_to = "day",
    values_to = "rainfall_mm"
  ) %>%
  mutate(day = as.integer(day))%>%
  select(day, rainfall_mm)

print(milford_FP_year12_muted)

# High exceedance future projection
milford_FP_year1633_highexceed <- milford_FP %>%
  filter(Year == 1633) %>%
  pivot_longer(
    cols = `1`:(ncol(milford_FP)-2),  # all day columns
    names_to = "day",
    values_to = "rainfall_mm"
  ) %>%
  mutate(day = as.integer(day))%>%
  select(day, rainfall_mm)

print(milford_FP_year1633_highexceed)


#Plotting time series
plotmilford_FP_year1633_highexceed <- ggplot(milford_FP_year1633_highexceed, aes(x=day, y=rainfall_mm)) +
  geom_line() + 
  labs(
    title = "High Exceedance Future Projection (Year 1633)",
    x = "Day of Year",
    y = "Rainfal (mm)"
  )+ 
  scale_y_continuous(limits = c(0, 260), breaks = seq(50, 250, by = 50))+
  theme_minimal() + 
  geom_hline(aes(yintercept = midRX1day_FP), col="#0072B2", linetype = "dashed", size = 1.1) 

plotmilford_FP_year12_muted <- ggplot(milford_FP_year12_muted, aes(x=day, y=rainfall_mm)) +
  geom_line() + 
  labs(
    title = "Muted Future Projection (Year 12)",
    x = "Day of Year",
    y = ""
  )+ 
  scale_y_continuous(limits = c(0, 260), breaks = seq(50, 250, by = 50))+
  theme_minimal() + 
  geom_hline(aes(yintercept = midRX1day_FP), col="#0072B2", linetype = "dashed", size = 1.1) 

plotmilford_FP_year1633_highexceed + plotmilford_FP_year12_muted

(plotmilford_CD_year1464_highexceed + plotmilford_CD_year18_muted) / (plotmilford_FP_year1633_highexceed + plotmilford_FP_year12_muted)



# 1 in 360 ----------------------------------------------------------------

# Identifying muted and high-exceedance years (Current Day, 1 in 360 days)
cd_exceed_onein360 <- count_exceedances_per_year(milford_CD, threshold_onein360)

exceedance_onein360_df_cd <- milford_CD %>%
  select(Year) %>%
  mutate(
    exceedances = cd_exceed_onein360,
  )

# Inspecting years - Current Day 
exceedance_onein360_df_cd %>%
  filter(exceedances==0)
#there are 1227 years where the onein360 threshold was not exceed. 
#Example year 5

exceedance_onein360_df_cd %>%
  arrange(desc(exceedances))
# Two years have 6 exceedances above the onein360 threshold. These are years 2112 and 2167

# Identifying muted and high-exceedance years (Future Projections, mid RX1day)
fp_exceed_onein360 <- count_exceedances_per_year(milford_FP, threshold_onein360_FP)

exceedance_onein360_df_fp <- milford_FP %>%
  select(Year) %>%
  mutate(
    exceedances = fp_exceed_onein360,
  )

# Inspect years - Future Projections
exceedance_onein360_df_fp %>%
  filter(exceedances==0)
#there are 986 years where the onein360 threshold was not exceed. 
#example year 12. 

exceedance_onein360_df_fp %>%
  arrange(desc(exceedances))
# 6 times is the largest number of exceedances this happens in the year 429

#Filtering years
#Muted curent day
milford_CD_year5_muted <- milford_CD %>%
  filter(Year == 5)
milford_CD_year5_muted

milford_CD_year5_muted <- milford_CD_year5_muted %>%
  select(starts_with("X")) %>%
  pivot_longer(
    cols = everything(),
    names_to = "day",
    values_to = "rainfall_mm"
  ) %>%
  mutate(
    day = as.integer(sub("X", "", day))  
  )
print(milford_CD_year5_muted)

#High exccedance current day
milford_CD_year2167_highexceed <- milford_CD %>%
  filter(Year == 2167)
milford_CD_year2167_highexceed

milford_CD_year2167_highexceed <- milford_CD_year2167_highexceed %>%
  select(starts_with("X")) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "day",
    values_to = "rainfall_mm"
  ) %>%
  mutate(
    day = as.integer(sub("X", "", day))  
  )
print(milford_CD_year2167_highexceed)

#Plotting time series
plotmilford_CD_year2167_highexceed <- ggplot(milford_CD_year2167_highexceed, aes(x=day, y=rainfall_mm)) +
  geom_line() +
  labs(
    title = "High Exceedance Current Day (Year 2167)",
    x = "Day of Year",
    y = "Daily Rainfall (mm)"
  )+
  scale_y_continuous(
    limits = c(0, 260),
    breaks = seq(50, 250, by = 50)
  )+
  theme_minimal() + 
  geom_hline(aes(yintercept = threshold_onein360), col="#E69F00", linetype = "dashed", size = 1.1) 
  
  plotmilford_CD_year5_muted <- ggplot(milford_CD_year5_muted, aes(x=day, y=rainfall_mm)) +
  geom_line() + 
  labs(
    title = "Muted Current Day (Year 5)",
    x = "Day of Year",
    y = ""
  )+
    scale_y_continuous(
      limits = c(0, 260),
      breaks = seq(50, 250, by = 50)
    )+
  theme_minimal() + 
  geom_hline(aes(yintercept = threshold_onein360), col="#E69F00", linetype = "dashed", size = 1.1) 
  
  plotmilford_CD_year2167_highexceed + plotmilford_CD_year5_muted

#Muted future projection
milford_FP_year12_muted <- milford_FP %>%
  filter(Year == 12) %>%
  pivot_longer(
    cols = `1`:(ncol(milford_FP)-2), 
    names_to = "day",
    values_to = "rainfall_mm"
  ) %>%
  mutate(day = as.integer(day))%>%
  select(day, rainfall_mm)

print(milford_FP_year12_muted)

# High exceedance future projection
milford_FP_year429_highexceed <- milford_FP %>%
  filter(Year ==429) %>%
  pivot_longer(
    cols = `1`:(ncol(milford_FP)-2),  # all day columns
    names_to = "day",
    values_to = "rainfall_mm"
  ) %>%
  mutate(day = as.integer(day))%>%
  select(day, rainfall_mm)

print(milford_FP_year429_highexceed)


#Plotting time series
plotmilford_FP_year429_highexceed <- ggplot(milford_FP_year429_highexceed, aes(x=day, y=rainfall_mm)) +
  geom_line() + 
  labs(
    title = "High Exceedance Future Projection (Year 429)",
    x = "Day of Year",
    y = "Rainfal (mm)"
  )+
  scale_y_continuous(
    limits = c(0, 260),
    breaks = seq(50, 250, by = 50)
  )+
  theme_minimal() + 
  geom_hline(aes(yintercept = threshold_onein360_FP), col="#E69F00", linetype = "dashed", size = 1.1) 
  
  plotmilford_FP_year12_muted_ <- ggplot(milford_FP_year12_muted, aes(x=day, y=rainfall_mm)) +
  geom_line() + 
  labs(
    title = "Muted Future Projection (Year 12)",
    x = "Day of Year",
    y = ""
  )+
    scale_y_continuous(
      limits = c(0, 260),
      breaks = seq(50, 250, by = 50)
    )+
  theme_minimal() + 
  geom_hline(aes(yintercept = threshold_onein360_FP), col= "#E69F00", linetype = "dashed", size = 1.1) 
  
  plotmilford_FP_year429_highexceed + plotmilford_FP_year12_muted_

(plotmilford_CD_year2167_highexceed + plotmilford_CD_year5_muted) / (plotmilford_FP_year429_highexceed + plotmilford_FP_year12_muted_)



# Correlation between exceedances and RX1day ------------------------------
  comparison_cd_360 <- milford_CD %>%
    mutate(
      exceedances = cd_exceed_onein360
    ) %>%
    select(Year, RX1day, exceedances)
  
  # Scatter + smooth
  ggplot(comparison_cd_360, aes(x = exceedances, y = RX1day)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, colour = "black") +
    theme_minimal() +
    labs(
      x = "Number of 1-in-360 exceedances per year",
      y = "RX1day (mm)",
      title = "Relationship between 1-in-360 exceedance frequency and RX1day (Current Day)"
    )
  
  # Spearman correlation
  cor(
    comparison_cd_360$exceedances,
    comparison_cd_360$RX1day,
    method = "spearman",
    use = "complete.obs"
  )
  
  
  comparison_cd_360 <- comparison_cd_360 %>%
    mutate(
      exceed_class = case_when(
        exceedances >= 4            ~ "High (≥4)",
        exceedances >= 1 & exceedances <= 3 ~ "Low (1–3)",
        TRUE                         ~ NA_character_
      )
    ) %>%
    filter(!is.na(exceed_class))
  
  ggplot(comparison_cd_360, aes(x = exceed_class, y = RX1day, fill = exceed_class)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.4) +
    theme_minimal() +
    labs(
      x = "Exceedance category",
      y = "RX1day (mm)",
      title = "RX1day for low vs high exceedance years (1-in-360-day threshold)"
    ) +
    guides(fill = "none")
  
  
  rx1_threshold <- quantile(comparison_cd_360$RX1day, 0.95)
  
  comparison_cd_360 %>%
    mutate(
      high_rx1 = RX1day >= rx1_threshold
    ) %>%
    group_by(exceed_class) %>%
    summarise(
      total_years = n(),
      high_rx1_years = sum(high_rx1),
      prop_high_rx1 = high_rx1_years / total_years
    )
  
  
  # Define high RX1day and high exceedances
  rx1_threshold <- quantile(comparison_cd_360$RX1day, 0.95)  # top 5% RX1day
  comparison_cd_360 <- comparison_cd_360 %>%
    mutate(
      high_rx1 = RX1day >= rx1_threshold,
      high_exceed = exceedances >= 4  # your "High (≥4)" category
    )
  
  # Build 2x2 contingency table
  table_cond <- table(
    High_RX1day = comparison_cd_360$high_rx1,
    High_Exceedances = comparison_cd_360$high_exceed
  )
  print(table_cond)
  
  # P(High RX1day | High exceedances)
  p_high_rx1_given_high_exceed <- table_cond["TRUE", "TRUE"] / sum(table_cond[, "TRUE"])
  print(p_high_rx1_given_high_exceed)
  
  # P(High exceedances | High RX1day)
  p_high_exceed_given_high_rx1 <- table_cond["TRUE", "TRUE"] / sum(table_cond["TRUE", ])
  print(p_high_exceed_given_high_rx1)
  


