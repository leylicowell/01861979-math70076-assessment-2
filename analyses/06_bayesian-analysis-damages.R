# in this script, we implement a bayesian regression model for economical damages
# using our merged data set

library(data.table) 
library(ggplot2) 
library(ggsci) 
library(hexbin)  
library(bayesplot) 
library(kableExtra)  
library(cmdstanr)  
library(webshot2) 
library(here)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(stringr)
library(geosphere)
library(corrr)
library(forcats)


# set colour scheme for Stan
bayesplot::color_scheme_set("brewer-RdYlBu")

#==============================================================================
# we load our data set
#==============================================================================

data <- read.csv(here("data", "derived", "merged-data.csv"))

#==============================================================================
# we make a pairwise correlation plot of wind at time of landfall, 
# overall curvature of hurricane path, duration of hurricane and adjusted total 
# damages
#==============================================================================

# reformat date and time to calculate duration
data <- data %>%
  mutate(TIME = sprintf("%04s", as.integer(TIME)),
         DATE_TIME = as.POSIXct(paste(YEAR, 
                                      MONTH, 
                                      DAY, 
                                      str_sub(TIME, 1, 2), 
                                      str_sub(TIME, 3, 4)),
                                format = "%Y %m %d %H %M"))


storm_predictors <- data %>%
  arrange(STORM_ID, DATE_TIME) %>%
  group_by(STORM_ID, NAME)%>%
  mutate(
    NEXT_LAT = lead(LAT),
    NEXT_LON = lead(LON),
    BEARING = bearing(cbind(LON, LAT), cbind(NEXT_LON, NEXT_LAT)),
    # change bearings from -180 to 180 to 0 to 360 degrees
    BEARING = BEARING + 180,
    DIFF_BEARING = abs(lead(BEARING)-BEARING)) %>%
  summarise(
    total_bearing_change = sum(DIFF_BEARING, na.rm = TRUE),
    avg_landfall_wind = mean(WIND[RECORD_ID =="L"]),
    total_damage = mean(ADJ_TOTAL_DAMAGE), # damage is the same for all rows for each storm
    duration_h = as.numeric(difftime(max(DATE_TIME), min(DATE_TIME), units = "hours")),
    .groups = "drop")

# extract relevant columns for pairwise correlation plot
corr_matrix <- storm_predictors %>%
  select("avg_landfall_wind", 
         "total_bearing_change", 
         "duration_h", 
         "total_damage")%>%
  correlate(diagonal = 1) %>%
  shave(upper = FALSE)

corr_matrix<- corr_matrix %>%
  pivot_longer(cols = -term,
               names_to = "colname",
               values_to = "corr") %>%
  mutate(rowname = fct_inorder(term),
         colname = fct_inorder(colname),
         label = if_else(is.na(corr), "", sprintf("%1.2f", corr)))

p1 <- ggplot(corr_matrix, aes(rowname, fct_rev(colname),
                 fill = corr)) +
  geom_tile() +
  geom_text(aes(
    label = label,
    color = abs(corr) < .75
  )) +
  coord_fixed(expand = FALSE) +
  scale_color_manual(
    values = c("white", "black"),
    guide = "none"
  ) +
  scale_fill_distiller(
    palette = "PuOr", na.value = "white",
    direction = 1, limits = c(-1, 1),
    name = "Pearson\nCorrelation:"
  ) +
  labs(x = NULL, y = NULL) +
  theme(panel.border = element_rect(color = NA, fill = NA),
        legend.position = c(.85, .8))

ggsave(here("outputs", "bayesian-analysis-damages", "correlation-plot.pdf"), 
       plot =p1,
       height = 7, 
       width = 7, 
       dpi=600)



    