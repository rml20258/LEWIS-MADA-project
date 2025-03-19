###############################
# analysis script
#
#this script loads the processed, cleaned data, does a simple analysis
#and saves the results to the results folder

#load needed packages. make sure they are installed.
library(ggplot2) #for plotting
library(broom) #for cleaning up output from lm()
library(here) #for data loading/saving
library(tidyverse)
library(patchwork)
library(dplyr)
library(tidyr)
library(gt)
#path to data
#note the use of the here() package and not absolute paths
data_location <- here::here("data","processed-data","project_data.rds")

#load data. 
project_data <- readRDS(data_location)


######################################
#Data fitting/statistical analysis
######################################






#Violin plots of the outcome by each categorical variable
demo1 <- project_data %>%
  ggplot(aes(coverage_diff, WHO_REGION)) + geom_violin() + labs(x = "Difference in vaccination coverage", y = "WHO region")

demo2 <- project_data %>%
  ggplot(aes(coverage_diff, WB_INCOME)) + geom_violin() + labs(x = "Difference in vaccination coverage", y = "Country income level")

demo3 <- project_data %>%
  ggplot(aes(coverage_diff, total_pop)) + geom_point() + geom_smooth() + labs(x = "Difference in vaccination coverage", y = "Total population (in thousands)")

demo4 <- project_data %>%
  ggplot(aes(coverage_diff, pop_density)) + geom_point() + geom_smooth() + labs(x = "Difference in vaccination coverage", y = "Population density (persons per square km)")

#Building 4 panel figure of country demographics
bivar_demo_1 <- (demo1 | demo2) / (demo3 | demo4)
bivar_demo_1_loc <- here("results","figures","bivar_demo_1.png")
# Save the plot created as a PNG
ggsave(bivar_demo_1_loc, plot = bivar_demo_1, width = 8, height = 6, dpi = 300)


demo5 <- project_data %>%
  ggplot(aes(coverage_diff, median_age)) + geom_point() + geom_smooth() + labs(x = "Difference in vaccination coverage", y = "Median age of population")

demo6 <- project_data %>%
  ggplot(aes(coverage_diff, fertility_rate)) + geom_point() + geom_smooth() + labs(x = "Difference in vaccination coverage", y = "Fertility rate (live births per woman)")

demo7 <- project_data %>%
  ggplot(aes(coverage_diff, death_rate)) + geom_point() + geom_smooth() + labs(x = "Difference in vaccination coverage", y = "Death rate (deaths per 1,000 population)")

demo8 <- project_data %>%
  ggplot(aes(coverage_diff, life_expect)) + geom_point() + geom_smooth() + labs(x = "Difference in vaccination coverage", y = "Life expectancy (in years)")

demo9 <- project_data %>%
  ggplot(aes(coverage_diff, migration)) + geom_point() + geom_smooth() + labs(x = "Difference in vaccination coverage", y = "Migration rate (per 1000 population)")

#Building 5 panel figure of the rest of the demographics
bivar_demo_2 <- (demo5 | demo6 | demo7) / (demo8 | demo9)
bivar_demo_2_loc <- here("results","figures","bivar_demo_2.png")
# Save the plot created as a PNG
ggsave(bivar_demo_2_loc, plot = bivar_demo_2, width = 8, height = 6, dpi = 300)

#Building 9 panel demographic figure
demo_fig_loc <- here("results","figures","demo_fig.png")
demo_fig <- (demo1 | demo2 | demo3) / (demo4 | demo5 | demo6) / (demo7 | demo8 | demo9) 
ggsave(demo_fig_loc, plot = demo_fig, width = 8, height = 10, dpi = 300)

#Vaccination variables
vac1 <- project_data %>%
  ggplot(aes(coverage_diff, HPV_NATIONAL_SCHEDULE)) + geom_violin() + labs(x = "Difference in vaccination coverage", y = "National vaccine schedule")

vac2 <- project_data %>%
  ggplot(aes(coverage_diff, HPV_PRIM_DELIV_STRATEGY)) + geom_violin() + labs(x = "Difference in vaccination coverage", y = "Primary vaccination delivery strategy")

vac3 <- project_data %>%
  ggplot(aes(coverage_diff, HPV_SEX)) + geom_violin() + labs(x = "Difference in vaccination coverage", y = "Genders vaccinated")

vac4 <- project_data %>%
  ggplot(aes(coverage_diff, as.factor(doses_rec))) + geom_violin() + labs(x = "Difference in vaccination coverage", y = "Number of recommended doses")

vac5 <- project_data %>%
  ggplot(aes(coverage_diff, HPV_YEAR_INTRODUCTION)) + geom_point() + geom_smooth() + labs(x = "Difference in vaccination coverage", y = "Year of vaccine introduction")

#Building 5 panel figure of vaccination program variables
bivar_vac <- (vac1 | vac2 | vac3) / (vac4 | vac5)
bivar_vac_loc <- here("results","figures","bivar_vac.png")
# Save the plot created as a PNG
ggsave(bivar_vac_loc, plot = bivar_vac, width = 8, height = 6, dpi = 300)

############################
#### First model fit
# fit linear model using coverage difference as outcome, COVID case rate as predictor

lmfit1 <- lm(coverage_diff ~ Cumulative_case_rate, project_data)  

# place results from fit into a data frame with the tidy function
lmtable1 <- broom::tidy(lmfit1)

#look at fit results
print(lmtable1)
fit_stat_1 <- glance(lmfit1)

# save fit results table  
table_file1 = here("results", "tables", "resulttable1.rds")
saveRDS(lmtable1, file = table_file1)

############################
#### Second and third model fit
# creating centered COVID case variables and seconrd order term
project_data <- project_data %>%
  mutate(covid_case_center = Cumulative_case_rate - mean(Cumulative_case_rate),
         covid_case_center_sq = covid_case_center * covid_case_center)

#Fit model with centered variable
lmfit2 <- lm(coverage_diff ~ covid_case_center, project_data)  
fit_stat_2 <- glance(lmfit2)
#Confirmed the model matches the uncentered model above

#Fit model with first and second order terms for COVID cases
lmfit3 <- lm(coverage_diff ~ covid_case_center + covid_case_center_sq, project_data)  
fit_stat_3 <- glance(lmfit3)

# place results from fit into a data frame with the tidy function
lmtable2 <- broom::tidy(lmfit2)
lmtable3 <- broom::tidy(lmfit3)

#look at fit results
print(lmtable2)
print(lmtable3)
#Model with a quadratic term for case count didn't help with model fit, so I'll stick with the first order model



# Next, fit full model and models with one additional predictor (national HPV schedule removed due to only 2 countries not having a national schedule)
# For each model, I run the model, convert it to the tidymodels framework, get select fit statistics, then 
# create a row with the beta, se, and p-value for the covid case count variable, plus the model r squared, AIC, and BIC for model comparison
lmfitfull <- lm(coverage_diff ~ Cumulative_case_rate + HPV_YEAR_INTRODUCTION +
                  delivery_strat + HPV_SEX + doses_rec + WHO_REGION + WB_INCOME + total_pop +  
                  total_pop + pop_density + median_age + pop_growth_rate + fertility_rate + death_rate +
                  life_expect + migration, project_data)  
lmtablefull <- broom::tidy(lmfitfull)
fit_stat_full <- glance(lmfitfull) %>% select(r.squared, AIC, BIC)
print(lmtablefull, n = 25)

full <- lmtablefull %>% 
  mutate(Adjustment = "Full model") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_full)

# HPV vax year of introduction 
lmfityear <- lm(coverage_diff ~ Cumulative_case_rate + HPV_YEAR_INTRODUCTION , project_data)  
lmtableyear <- broom::tidy(lmfityear)
fit_stat_year <- glance(lmfityear) %>% select(r.squared, AIC, BIC)
year <- lmtableyear %>% 
  mutate(Adjustment = "Year of vaccine introduction") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_year)

# Delivery strategy 
lmfitdeliv <- lm(coverage_diff ~ Cumulative_case_rate  + delivery_strat, project_data)  
lmtabledeliv <- broom::tidy(lmfitdeliv)
fit_stat_deliv <- glance(lmfitdeliv) %>% select(r.squared, AIC, BIC)
deliv <- lmtabledeliv %>% 
  mutate(Adjustment = "Primary delivery method") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_deliv)

# Gender nuetral 
lmfitgen <- lm(coverage_diff ~ Cumulative_case_rate + HPV_SEX, project_data)  
lmtablegen <- broom::tidy(lmfitgen)
fit_stat_gen <- glance(lmfitgen) %>% select(r.squared, AIC, BIC)
gen <- lmtablegen %>% 
  mutate(Adjustment = "Genders vaccinated") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_gen)

# Number of recommended doses
lmfitdoses <- lm(coverage_diff ~ Cumulative_case_rate + doses_rec, project_data)  
lmtabledoses <- broom::tidy(lmfitdoses)
fit_stat_doses <- glance(lmfitdoses) %>% select(r.squared, AIC, BIC)
doses <- lmtabledoses %>% 
  mutate(Adjustment = "Number of recommended doses") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_doses)

# WHO region 
lmfitwho <- lm(coverage_diff ~ Cumulative_case_rate + WHO_REGION, project_data)  
lmtablewho <- broom::tidy(lmfitwho)
fit_stat_who <- glance(lmfitwho) %>% select(r.squared, AIC, BIC)
who <- lmtablewho %>% 
  mutate(Adjustment = "WHO region") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_who)

# Income level
lmfitincome <- lm(coverage_diff ~ Cumulative_case_rate + WB_INCOME, project_data)  
lmtableincome <- broom::tidy(lmfitincome)
fit_stat_income <- glance(lmfitincome) %>% select(r.squared, AIC, BIC)
income <- lmtableincome %>% 
  mutate(Adjustment = "Country income level") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_income)

# total population 
lmfitpop <- lm(coverage_diff ~ Cumulative_case_rate + total_pop, project_data)  
lmtablepop <- broom::tidy(lmfitpop)
fit_stat_pop <- glance(lmfitpop) %>% select(r.squared, AIC, BIC)
pop <- lmtablepop %>% 
  mutate(Adjustment = "Total population") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_pop)

# population density 
lmfitpopdens <- lm(coverage_diff ~ Cumulative_case_rate + pop_density, project_data)  
lmtablepopdens <- broom::tidy(lmfitpopdens)
fit_stat_popdens <- glance(lmfitpopdens) %>% select(r.squared, AIC, BIC)
popdens <- lmtablepopdens %>% 
  mutate(Adjustment = "Population density") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_popdens)

# median age
lmfitage <- lm(coverage_diff ~ Cumulative_case_rate + median_age, project_data)  
lmtableage <- broom::tidy(lmfitage)
fit_stat_age <- glance(lmfitage) %>% select(r.squared, AIC, BIC)
age <- lmtableage %>% 
  mutate(Adjustment = "Median age") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_age)

# population growth rate 
lmfitpopgrow <- lm(coverage_diff ~ Cumulative_case_rate + pop_growth_rate, project_data)  
lmtablepopgrow <- broom::tidy(lmfitpopgrow)
fit_stat_popgrow <- glance(lmfitpopgrow) %>% select(r.squared, AIC, BIC)
popgrow <- lmtablepopgrow %>% 
  mutate(Adjustment = "Population growth rate") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_popgrow)

# fertility rate 
lmfitfertility <- lm(coverage_diff ~ Cumulative_case_rate + fertility_rate, project_data)  
lmtablefertility <- broom::tidy(lmfitfertility)
fit_stat_fertility <- glance(lmfitfertility) %>% select(r.squared, AIC, BIC)
fertility <- lmtablefertility %>% 
  mutate(Adjustment = "Fertility rate") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_fertility)

# death rate 
lmfitdeath <- lm(coverage_diff ~ Cumulative_case_rate + death_rate, project_data)  
lmtabledeath <- broom::tidy(lmfitdeath)
fit_stat_death <- glance(lmfitdeath) %>% select(r.squared, AIC, BIC)
death <- lmtabledeath %>% 
  mutate(Adjustment = "Death rate") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_death)

# life expectancy 
lmfitlifeexp <- lm(coverage_diff ~ Cumulative_case_rate + life_expect, project_data)  
lmtablelifeexp <- broom::tidy(lmfitlifeexp)
fit_stat_lifeexp <- glance(lmfitlifeexp) %>% select(r.squared, AIC, BIC)
lifeexp <- lmtablelifeexp %>% 
  mutate(Adjustment = "Life expectancy") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_lifeexp)

# migration rate 
lmfitmig <- lm(coverage_diff ~ Cumulative_case_rate + migration, project_data)  
lmtablemig <- broom::tidy(lmfitmig)
fit_stat_mig <- glance(lmfitmig) %>% select(r.squared, AIC, BIC)
mig <- lmtablemig %>% 
  mutate(Adjustment = "Migration rate") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_mig)

# unadjusted
lmfitunadj <- lm(coverage_diff ~ Cumulative_case_rate , project_data)  
lmtableunadj <- broom::tidy(lmfitunadj)
fit_stat_unadj <- glance(lmfitunadj) %>% select(r.squared, AIC, BIC)
unadj <- lmtableunadj %>% 
  mutate(Adjustment = "Unadjusted") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_unadj)

table_model_adj <- as_tibble(rbind(full, year, deliv, gen, doses, who, income, pop, popdens, age, popgrow, fertility, death, lifeexp, mig, unadj)) %>%
  mutate(estimate = round(estimate, digits = 6), 
         std.error = round(std.error, digits = 5), 
         p.value = num(p.value, digits = 3), 
         r.squared = num(r.squared, digits = 3), 
         AIC = num(AIC, digits = 1), 
         BIC = num(BIC, digits = 1)) %>%
  rename(beta = estimate)
  

# Find the minimum AIC and BIC values
min_AIC <- min(table_model_adj$AIC, na.rm = TRUE)
min_BIC <- min(table_model_adj$BIC, na.rm = TRUE)
min_rsq <- min(table_model_adj$r.squared, na.rm = TRUE)


# Create a GT table with conditional formatting
table_model_adj_format <- table_model_adj %>%
  gt() %>%
  data_color(
    columns = AIC,
    fn = scales::col_numeric(palette = c("lightgreen", "white"), domain = c(min_AIC, max(table_model_adj$AIC)))
  ) %>%
  data_color(
    columns = BIC,
    fn = scales::col_numeric(palette = c("lightblue", "white"), domain = c(min_BIC, max(table_model_adj$BIC)))
  ) %>%
  data_color(
    columns = r.squared,
    fn = scales::col_numeric(palette = c("white", "lightpink"), domain = c(min_rsq, max(table_model_adj$r.squared)))
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "green"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(columns = AIC, rows = AIC == min_AIC)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "skyblue"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(columns = BIC, rows = BIC == min_BIC)
  ) %>%
  tab_style(
    style = list(
      cell_fill(color = "pink"),
      cell_text(weight = "bold")
    ),
    locations = cells_body(columns = r.squared, rows = r.squared == max(table_model_adj$r.squared))
  ) 

#Saving formatted table for all model comparisons
gtsave(table_model_adj_format, filename = here("results","tables","table_model_adj_format.png"))

#Converting full table to GT table to save as an image
lmtablefull_gt <- lmtablefull %>%
  gt()
#Saving formatted table for all model comparisons
gtsave(lmtablefull_gt, filename = here("results","tables","lmtablefull_gt.png"))


# save fit results table  
#table_file2 = here("results", "tables", "resulttable2.rds")
#saveRDS(lmtable2, file = table_file2)


# Leaving the code below here in case I want to add a table of mean coverage diff and correlations by demographics

# 
# # Define categorical and continuous variables
# categorical_vars <- c("WHO_REGION", "WB_INCOME", "HPV_NATIONAL_SCHEDULE", 
#                       "HPV_YEAR_INTRODUCTION", "HPV_PRIM_DELIV_STRATEGY", 
#                       "HPV_SEX", "doses_rec")
# 
# continuous_vars <- c("Cumulative_case_rate", "Cumulative_death_rate", 
#                      "total_pop", "pop_density", "median_age", 
#                      "pop_growth_rate", "fertility_rate", "death_rate", 
#                      "life_expect", "migration")
# 
# # Assuming `project_data` contains all variables
# # Function to compute Mean (SD) and ANOVA p-value for categorical variables
# categorical_summary <- categorical_vars %>%
#   lapply(function(var) {
#     # Run ANOVA
#     anova_result <- aov(coverage_diff ~ get(var), data = project_data)
#     p_value <- summary(anova_result)[[1]][["Pr(>F)"]][1]  # Extract p-value
#     
#     # Compute Mean (SD) for each category
#     project_data %>%
#       group_by(across(all_of(var))) %>%
#       summarise(
#         statistic = paste0(round(mean(coverage_diff, na.rm = TRUE), 2), 
#                            " (", round(sd(coverage_diff, na.rm = TRUE), 2), ")"),
#         .groups = "drop"
#       ) %>%
#       rename(category = 1) %>%  # Rename the first column to 'category'
#       mutate(variable = var, p_value = round(p_value, 4))
#   }) %>%
#   bind_rows() %>%
#   select(variable, category, statistic, p_value)
# 
# # Function to compute correlation and p-value for continuous variables
# continuous_summary <- continuous_vars %>%
#   lapply(function(var) {
#     cor_test <- cor.test(project_data[[var]], project_data$coverage_diff, use = "pairwise.complete.obs")
#     
#     tibble(
#       variable = var,
#       category = NA,  # No category for continuous variables
#       statistic = round(cor_test$estimate, 2),
#       p_value = round(cor_test$p.value, 4)
#     )
#   }) %>% 
#     mutate(category = as.character(category)) %>%
#   bind_rows()
# str(categorical_summary)
# # Combine results into a single table
# final_table <- bind_rows(categorical_summary, continuous_summary) %>%
#   arrange(variable, category)
# 
# # Print the table
# print(final_table)