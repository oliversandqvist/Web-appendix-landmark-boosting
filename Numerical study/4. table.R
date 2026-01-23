library(dplyr)
library(tidyr)
library(stringr)
library(xtable)

##Computational details
sessionInfo()
parallel::detectCores()
Sys.info()
system("lscpu | grep 'Model name'")

##Training time table
extract_times <- function(model_list, scenario_name, method_type) {
  
  do.call(rbind, lapply(names(model_list), function(nm) {
    fit <- model_list[[nm]]
    parts <- str_split(nm, "_")[[1]]
    model_algo <- ifelse(grepl("xgb", parts[1]), "Boosted Trees", "Cox")
    
    n_val <- as.numeric(parts[2])
    q_val <- if(method_type == "Landmark") as.numeric(parts[3]) else NA
    
    time_val <- as.numeric(fit$model$trainTime, units = "secs") 
    if(is.null(time_val)) time_val <- as.numeric(fit$trainTime, units="secs")
    
    data.frame(
      Scenario = scenario_name,
      Model = model_algo,
      Method = method_type,
      n = n_val,
      Q = q_val,
      Time = time_val,
      stringsAsFactors = FALSE
    )
  }))
}

df_a <- bind_rows(
  extract_times(coxNaive_linear, "Scenario 1 (Linear)", "Naive"),
  extract_times(coxStack_linear, "Scenario 1 (Linear)", "Landmark"),
  extract_times(xgbNaive_linear, "Scenario 1 (Linear)", "Naive"),
  extract_times(xgbStack_linear, "Scenario 1 (Linear)", "Landmark")
)

df_b <- bind_rows(
  extract_times(coxNaive_nonlinear, "Scenario 2 (Nonlinear)", "Naive"),
  extract_times(coxStack_nonlinear, "Scenario 2 (Nonlinear)", "Landmark"),
  extract_times(xgbNaive_nonlinear, "Scenario 2 (Nonlinear)", "Naive"),
  extract_times(xgbStack_nonlinear, "Scenario 2 (Nonlinear)", "Landmark")
)

df_c <- bind_rows(
  extract_times(coxNaive_large, "Scenario 3 (High-Dim)", "Naive"),
  extract_times(coxStack_large, "Scenario 3 (High-Dim)", "Landmark"),
  extract_times(xgbNaive_large, "Scenario 3 (High-Dim)", "Naive"),
  extract_times(xgbStack_large, "Scenario 3 (High-Dim)", "Landmark")
)

all_times <- bind_rows(df_a, df_b, df_c)

table_data <- all_times %>%
  filter(Method == "Naive" | (Method == "Landmark" & Q == 10)) %>%
  mutate(MethodLabel = paste(Model, ifelse(Method == "Naive", "(Naive)", "(Q=10)"))) %>%
  select(Scenario, n, MethodLabel, Time) %>%
  pivot_wider(names_from = MethodLabel, values_from = Time) %>%
  arrange(Scenario, n)

table_data <- table_data %>%
  mutate(across(where(is.numeric) & !c(n), ~ sprintf("%.2f", .x))) %>%
  mutate(n = scales::comma(n))

print_latex_table <- function(df) {
  xt <- xtable(df, caption = "Training times (in seconds) for the Naive and Landmark (Q=10) estimators across varying sample sizes $n$.")
  
  print(xt, 
        include.rownames = FALSE, 
        booktabs = TRUE, 
        sanitize.text.function = function(x) x, 
        floating = FALSE) 
}

print_latex_table(table_data)