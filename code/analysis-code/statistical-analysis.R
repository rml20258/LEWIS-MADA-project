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
library(yardstick)

library(tidymodels)
library(here)
library(gtsummary)
library(caret)
library(pROC)
library(quarto)
library(GGally)
library(tune)
library(workflowsets)
library(parsnip)

library(glmnet) #LASSO
library(ranger) #Random forest plots

#Classification tree packages
library(rpart)
library(rpart.plot)


#path to data
#note the use of the here() package and not absolute paths
data_location <- here::here("data","processed-data","project_data.rds")

#load data. 
project_data <- readRDS(data_location)

#Restricting project_data to complete rows 
cc_project_data <- project_data %>%
  drop_na()

######################################
#Data fitting/statistical analysis
######################################



############################
#### First model fit
# fit linear model using coverage difference as outcome, COVID case rate as predictor

lmfit1 <- lm(coverage_diff ~ Cumulative_case_rate, cc_project_data)  

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
cc_project_data <- cc_project_data %>%
  mutate(covid_case_center = Cumulative_case_rate - mean(Cumulative_case_rate),
         covid_case_center_sq = covid_case_center * covid_case_center)

#Fit model with centered variable
lmfit2 <- lm(coverage_diff ~ covid_case_center, cc_project_data)  
fit_stat_2 <- glance(lmfit2)
#Confirmed the model matches the uncentered model above

#Fit model with first and second order terms for COVID cases
lmfit3 <- lm(coverage_diff ~ covid_case_center + covid_case_center_sq, cc_project_data)  
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
                  life_expect + migration, cc_project_data)  
lmtablefull <- broom::tidy(lmfitfull)
fit_stat_full <- glance(lmfitfull) %>% select(r.squared, AIC, BIC)
print(lmtablefull, n = 25)

full <- lmtablefull %>% 
  mutate(Adjustment = "Full model") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_full)

# HPV vax year of introduction 
lmfityear <- lm(coverage_diff ~ Cumulative_case_rate + HPV_YEAR_INTRODUCTION , cc_project_data)  
lmtableyear <- broom::tidy(lmfityear)
fit_stat_year <- glance(lmfityear) %>% select(r.squared, AIC, BIC)
year <- lmtableyear %>% 
  mutate(Adjustment = "Year of vaccine introduction") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_year)

# Delivery strategy 
lmfitdeliv <- lm(coverage_diff ~ Cumulative_case_rate  + delivery_strat, cc_project_data)  
lmtabledeliv <- broom::tidy(lmfitdeliv)
fit_stat_deliv <- glance(lmfitdeliv) %>% select(r.squared, AIC, BIC)
deliv <- lmtabledeliv %>% 
  mutate(Adjustment = "Primary delivery method") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_deliv)

# Gender nuetral 
lmfitgen <- lm(coverage_diff ~ Cumulative_case_rate + HPV_SEX, cc_project_data)  
lmtablegen <- broom::tidy(lmfitgen)
fit_stat_gen <- glance(lmfitgen) %>% select(r.squared, AIC, BIC)
gen <- lmtablegen %>% 
  mutate(Adjustment = "Genders vaccinated") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_gen)

# Number of recommended doses
lmfitdoses <- lm(coverage_diff ~ Cumulative_case_rate + doses_rec, cc_project_data)  
lmtabledoses <- broom::tidy(lmfitdoses)
fit_stat_doses <- glance(lmfitdoses) %>% select(r.squared, AIC, BIC)
doses <- lmtabledoses %>% 
  mutate(Adjustment = "Number of recommended doses") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_doses)

# WHO region 
lmfitwho <- lm(coverage_diff ~ Cumulative_case_rate + WHO_REGION, cc_project_data)  
lmtablewho <- broom::tidy(lmfitwho)
fit_stat_who <- glance(lmfitwho) %>% select(r.squared, AIC, BIC)
who <- lmtablewho %>% 
  mutate(Adjustment = "WHO region") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_who)

# Income level
lmfitincome <- lm(coverage_diff ~ Cumulative_case_rate + WB_INCOME, cc_project_data)  
lmtableincome <- broom::tidy(lmfitincome)
fit_stat_income <- glance(lmfitincome) %>% select(r.squared, AIC, BIC)
income <- lmtableincome %>% 
  mutate(Adjustment = "Country income level") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_income)

# total population 
lmfitpop <- lm(coverage_diff ~ Cumulative_case_rate + total_pop, cc_project_data)  
lmtablepop <- broom::tidy(lmfitpop)
fit_stat_pop <- glance(lmfitpop) %>% select(r.squared, AIC, BIC)
pop <- lmtablepop %>% 
  mutate(Adjustment = "Total population") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_pop)

# population density 
lmfitpopdens <- lm(coverage_diff ~ Cumulative_case_rate + pop_density, cc_project_data)  
lmtablepopdens <- broom::tidy(lmfitpopdens)
fit_stat_popdens <- glance(lmfitpopdens) %>% select(r.squared, AIC, BIC)
popdens <- lmtablepopdens %>% 
  mutate(Adjustment = "Population density") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_popdens)

# median age
lmfitage <- lm(coverage_diff ~ Cumulative_case_rate + median_age, cc_project_data)  
lmtableage <- broom::tidy(lmfitage)
fit_stat_age <- glance(lmfitage) %>% select(r.squared, AIC, BIC)
age <- lmtableage %>% 
  mutate(Adjustment = "Median age") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_age)
#Saving results of median age adjusted model to call in the supplementary material
lmtableage_loc = here("results", "tables", "lmtableage.rds") #Saving table as rds
saveRDS(lmtableage, file = lmtableage_loc)

# population growth rate 
lmfitpopgrow <- lm(coverage_diff ~ Cumulative_case_rate + pop_growth_rate, cc_project_data)  
lmtablepopgrow <- broom::tidy(lmfitpopgrow)
fit_stat_popgrow <- glance(lmfitpopgrow) %>% select(r.squared, AIC, BIC)
popgrow <- lmtablepopgrow %>% 
  mutate(Adjustment = "Population growth rate") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_popgrow)

# fertility rate 
lmfitfertility <- lm(coverage_diff ~ Cumulative_case_rate + fertility_rate, cc_project_data)  
lmtablefertility <- broom::tidy(lmfitfertility)
fit_stat_fertility <- glance(lmfitfertility) %>% select(r.squared, AIC, BIC)
fertility <- lmtablefertility %>% 
  mutate(Adjustment = "Fertility rate") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_fertility)

# death rate 
lmfitdeath <- lm(coverage_diff ~ Cumulative_case_rate + death_rate, cc_project_data)  
lmtabledeath <- broom::tidy(lmfitdeath)
fit_stat_death <- glance(lmfitdeath) %>% select(r.squared, AIC, BIC)
death <- lmtabledeath %>% 
  mutate(Adjustment = "Death rate") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_death)

# life expectancy 
lmfitlifeexp <- lm(coverage_diff ~ Cumulative_case_rate + life_expect, cc_project_data)  
lmtablelifeexp <- broom::tidy(lmfitlifeexp)
fit_stat_lifeexp <- glance(lmfitlifeexp) %>% select(r.squared, AIC, BIC)
lifeexp <- lmtablelifeexp %>% 
  mutate(Adjustment = "Life expectancy") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_lifeexp)

# migration rate 
lmfitmig <- lm(coverage_diff ~ Cumulative_case_rate + migration, cc_project_data)  
lmtablemig <- broom::tidy(lmfitmig)
fit_stat_mig <- glance(lmfitmig) %>% select(r.squared, AIC, BIC)
mig <- lmtablemig %>% 
  mutate(Adjustment = "Migration rate") %>%
  filter(term == "Cumulative_case_rate") %>%
  select(Adjustment, estimate, std.error, p.value) %>%
  cbind(fit_stat_mig)

# unadjusted
lmfitunadj <- lm(coverage_diff ~ Cumulative_case_rate , cc_project_data)  
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
  )  %>%
  cols_label(
    beta = "Beta for COVID case burden",  # Rename columns
    std.error = "Standard error", 
    p.value = "p-value",
    r.squared = "R-squared"
  )

#Saving formatted table for all model comparisons
gtsave(table_model_adj_format, filename = here("results","tables","table_model_adj_format.png"))

#Saving table so the numbers can be referenced in text
table_model_adj_format_loc = here("results", "tables", "table_model_adj_format.rds") #Saving table as rds
saveRDS(table_model_adj_format, file = table_model_adj_format_loc)

######### Model fit plots ###########
# Compute predictions on training data
lmfitfull_pred <- predict(lmfitfull, cc_project_data) %>% bind_cols(cc_project_data)
# Compute RMSE for both models
lmfitfull_rmse <- rmse(lmfitfull_pred, truth = coverage_diff, estimate = ...1)

# Create the scatter plot of observed vs predicted
full_pred_vs_obs <- ggplot(lmfitfull_pred, aes(y = coverage_diff, x = ...1)) +
  geom_point(color = "blue", alpha = 0.6) +  # Plot points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal fit line
  labs(title = "Observed vs. Predicted Values",
       x = "Observed",
       y = "Predicted") +
  theme_minimal()

# Save the predicted vs observed plot as a PNG
full_pred_vs_obs_loc <- here("results","figures","full_pred_vs_obs.png")
ggsave(full_pred_vs_obs_loc, plot = full_pred_vs_obs, width = 6, height = 3, dpi = 300)


# Create the scatter plot of residuals vs predicted
full_pred_vs_res <- ggplot(lmfitfull_pred, aes(x = ...1, y = lmfitfull$residuals)) +
  geom_point() +  # Scatter plot 
  geom_abline(slope = 0, intercept = 0, linetype = "solid", color = "black") +  # line at 0
  theme_minimal() +  # Use a minimal theme
  labs(title = "Predicted values vs. residuals",
       x = "Predicted Values",
       y = "Residuals")

# Save the 2 panel figure of fully adjusted linear model as a PNG
full_lm_diagnostics <- (full_pred_vs_obs|full_pred_vs_res)
full_lm_diagnostics_loc <- here("results","figures","full_lm_diagnostics.png")
ggsave(full_lm_diagnostics_loc, plot = full_lm_diagnostics, width = 8, height = 3, dpi = 300)

############ LASSO regression ################

set.seed(12681)

# Define a grid of penalty values (lambda) from 1E-5 to 1E2 on a log scale
lambda_project_grid <- 10^seq(log10(1E-5), log10(1E2), length.out = 50)

# Define LASSO model with a tunable penalty
lasso_project_model <- linear_reg(penalty = tune()) %>% 
  set_engine("glmnet") %>% 
  set_mode("regression")

# Create workflow
lasso_project_wf <- workflow() %>% 
  add_model(lasso_project_model) %>% 
  add_formula(coverage_diff ~ Cumulative_case_rate + HPV_YEAR_INTRODUCTION +
                delivery_strat + HPV_SEX + doses_rec + WHO_REGION + WB_INCOME + total_pop +  
                total_pop + pop_density + median_age + pop_growth_rate + fertility_rate + death_rate +
                life_expect + migration)

# Create resampling object using apparent()
resamples_project <- apparent(cc_project_data)

# Perform tuning with tune_grid()
lasso_tune_project_results <- tune_grid(
  lasso_project_wf,
  resamples = resamples_project,
  grid = tibble(penalty = lambda_project_grid),
  metrics = metric_set(yardstick::rsq),
  control = control_grid(save_pred = TRUE)  # Ensure predictions are stored
)

# Collect tuning results
#lasso_results_project <- collect_metrics(lasso_tune_project_results)
rsq <- lasso_tune_project_results$.metrics[[1]]


#Graphing rsq by penalty
lasso_rsq_plot <- ggplot(rsq, aes(x = penalty, y = .estimate, color = .config)) +
  geom_point() +
  scale_x_log10() +  # Use log scale for penalty
  labs(title = "LASSO Tuning - R-squared vs. Penalty",
       x = "Penalty (log scale)",
       y = "R-squared") +
  theme_minimal() + theme(legend.position = "none")


# Create the scatter plot of observed vs predicted
lasso_pred_obs_plot <- ggplot(lasso_tune_project_results[[5]][[1]], aes(x = coverage_diff, y = .pred, color = .config)) +
  geom_point( alpha = 0.6) +  # Plot points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # Ideal fit line
  labs(title = "Observed vs. Predicted Values - LASSO",
       x = "Observed",
       y = "Predicted") +
  theme_minimal() + theme(legend.position = "none")


#Building 2 panel figure of LASSO diagnostics
lasso_fit_plot <- (lasso_rsq_plot | lasso_pred_obs_plot) 
lasso_fit_plot_loc <- here("results","figures","lasso_fit_plot.png")
# Save the plot created as a PNG
ggsave(lasso_fit_plot_loc, plot = lasso_fit_plot, width = 8, height = 3, dpi = 300)

### I want to compare the betas for the LASSO model with tuning parameter = 1e-2 to the linear regression model
# Finalize the workflow with the chosen lambda
final_wf <- finalize_workflow(lasso_project_wf, tibble(penalty = 1e-2))

# Fit the model using the full dataset
final_fit <- fit(final_wf, data = cc_project_data)

# Extract the fitted LASSO model
lasso_fit <- extract_fit_parsnip(final_fit)

# Extracting the beta coefficients from the selected model
lasso_coefs <- tidy(lasso_fit)
print(lasso_coefs, n = 23)


######### Creating table with linear regression and LASSO betas #########


#Converting full table to GT table to save as an image
lmtablefull_gt <- lmtablefull %>%
  mutate(lasso = lasso_coefs$estimate) %>%
  gt() %>%
  fmt_number(
    columns = everything(),  # Apply to all numeric columns
    n_sigfig = 3             # Adjust decimal places as needed
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),  # Make headers bold
    locations = cells_column_labels()
  ) %>%
  cols_label(
    term = "Coefficient",  # Rename columns
    estimate = "Beta estimate", 
    std.error = "Standard error",
    statistic = "t-statistic", 
    p.value = "p-value", 
    lasso = "LASSO beta estimate"
  ) %>%
  tab_spanner(
    label = md("**Linear regression model**"),  # denotes columns from the linear regression model
    columns = c("estimate", "std.error", "statistic", "p.value")  # Columns to span
  ) %>% 
  text_case_match(
    "Cumulative_case_rate" ~ "Cumulative COVID cases",
    "HPV_YEAR_INTRODUCTION" ~ "Year of HPV vaccine introduction",
    "delivery_stratMixed" ~ "Mixed delivery strategy",
    "delivery_stratSchool-based" ~ "School-based delivery strategy",
    "HPV_SEXFemale" ~ "Female only vaccination",
    "doses_rec" ~ "Number of recommended doses",
    "WHO_REGIONAMR" ~ "Region of the Americas",
    "WHO_REGIONEMR" ~ "Eastern Mediterranean region",
    "WHO_REGIONEUR" ~ "European region",
    "WHO_REGIONSEAR" ~ "South-East Asian region",
    "WHO_REGIONWPR" ~ "Western Pacific region",
    "WB_INCOMELow income" ~ "Low income level",
    "WB_INCOMELower middle income" ~ "Lower middle income level",
    "WB_INCOMEUpper middle income" ~ "Upper middle income level",
    "total_pop" ~ "Total population",
    "pop_density" ~ "Population density",
    "median_age" ~ "Median age",
    "pop_growth_rate" ~ "Population growth rate",
    "fertility_rate" ~ "Total fertility rate",
    "death_rate" ~ "Crude death rate",
    "life_expect" ~ "Life expectancy at birth",
    "migration" ~ "Net migration rate"
  )

#Saving formatted table for all model comparisons
gtsave(lmtablefull_gt, filename = here("results","tables","lmtablefull_gt.png"))

#Creating and saving results as table so they can be referenced in the manuscript
lmtablefull_gt_ref <- lmtablefull %>%
  mutate(lasso = lasso_coefs$estimate) #Creating table to be referenced

lmtablefull_gt_ref_loc = here("results", "tables", "lmtablefull_gt_ref.rds") #Saving table as rds
saveRDS(lmtablefull_gt_ref, file = lmtablefull_gt_ref_loc)

############# Classification tree #############


#### Classification tree
# Renaming some variables so that the tree has better formatting
cc_project_data_tree <- cc_project_data %>%
  rename(`Median age` = median_age,
         `Cumulative COVID cases` = Cumulative_case_rate,
         `Poplation density` = pop_density,
         `Fertility rate` = fertility_rate,
         `Migration rate` = migration,
         `Number of recommended doses` = doses_rec)

# Fit a regression tree model
reg_tree <- rpart(coverage_diff ~ total_pop + `Poplation density` + `Median age` + pop_growth_rate + 
                    `Fertility rate` + death_rate + life_expect + `Migration rate` + `Cumulative COVID cases` + 
                    HPV_YEAR_INTRODUCTION + WHO_REGION + WB_INCOME + HPV_NATIONAL_SCHEDULE + 
                    HPV_YEAR_INTRODUCTION + HPV_PRIM_DELIV_STRATEGY + HPV_SEX + `Number of recommended doses`, 
                  data = cc_project_data_tree, 
                  method = "anova")  # anova method for regression trees



# Visualize and save the tree
png(here("results", "figures", "regression_tree.png"), width = 800, height = 600)  # Adjust size as needed

# Plot the tree
tree <- rpart.plot(reg_tree, type = 2, extra = 101, fallen.leaves = TRUE,
           main = "Regression Tree for Coverage Difference")

# Close the graphics device
dev.off()

#Saving formatted table for all model comparisons
#ggsave(here("results","tables","tree.png"), plot = tree, width = 6, height = 3, dpi = 300)

# Print a summary of the tree
summary(reg_tree)

printcp(reg_tree)


# Calculate R-squared for the tree
# Get predicted values from the model
predictions <- predict(reg_tree, newdata = cc_project_data_tree)

# Actual values (target variable)
actual_values <- cc_project_data_tree$coverage_diff

# Calculate the residual sum of squares (RSS)
rss <- sum((actual_values - predictions)^2)

# Calculate the total sum of squares (TSS)
tss <- sum((actual_values - mean(actual_values))^2)

# Calculate R-squared
tree_r_squared <- round(1 - (rss / tss),3)

# Print R-squared
print(paste("R-squared: ", round(tree_r_squared, 4)))

########## Table of R-squared values from each model ###########
r_sq_lm <- glance(lmfitfull) %>% select(r.squared) %>% as.numeric() %>% round(3) #R-squared from full lienar regression model
r_sq_medianage <- glance(lmfitage) %>% select(r.squared) %>% as.numeric() %>% round(3) #R-squared from median age linear regression model
r_sq_lasso <- rsq %>%
  dplyr::filter(penalty == 0.1) %>%
  dplyr::select(.estimate) %>% 
  as.numeric() %>% 
  round(3)

# R-squared from tree is above (tree_r_squared)
r_sq_table <- tibble(Model = c("Full multivariable linear regression","Median age adjusted linear regression", "LASSO regression", "Regression tree"),
                     `R-squared` = c(r_sq_lm, r_sq_medianage, r_sq_lasso,tree_r_squared))

r_sq_table_loc = here("results", "tables", "r_sq_table.rds")
saveRDS(r_sq_table, file = r_sq_table_loc)

