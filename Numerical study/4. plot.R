setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library("dplyr")
library("xgboost")
library("survival")
library("ggplot2")
library("tidyr")
library("stringr")
library("gridExtra") 

options(scipen = 999)

dataSurvPred_linear <- readRDS("Results/dataSurvPred_linear.rds")
dataSurvPred_nonlinear <- readRDS("Results/dataSurvPred_nonlinear.rds")
dataSurvPred_large <- readRDS("Results/dataSurvPred_large.rds")

dataSurvPred_large[is.nan(as.matrix(dataSurvPred_large))] <- 0 #NaN occur as 0/0

calc_metrics <- function(df, scenario_label) {
  truth <- df$surv
  
  pred_cols <- names(df)[grep("^(cox|xgb)", names(df))]
  
  results_list <- lapply(pred_cols, function(col_name) {
    preds <- df[[col_name]]
    rmse_val <- sqrt(mean((preds - truth)^2, na.rm = TRUE))
    mape_val <- mean(abs(preds - truth)/truth)
    
    parts <- str_split(col_name, "_")[[1]]
    model_type <- ifelse(grepl("xgb", parts[1]), "Boosted Trees", "Cox")
    method <- ifelse(grepl("Naive", parts[1]), "Naive", "Landmark")
    
    n_val <- as.numeric(parts[2])
    q_val <- if (method == "Landmark") as.numeric(parts[3]) else NA
    
    data.frame(
      Scenario = scenario_label,
      Model = model_type,
      Method = method,
      n = n_val,
      Q = q_val,
      RMSE = rmse_val,
      MAPE = mape_val,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(results_list)
}

df_linear <- calc_metrics(dataSurvPred_linear, "Scenario 1 (Linear)")
df_nonlinear <- calc_metrics(dataSurvPred_nonlinear, "Scenario 2 (Nonlinear)")
df_large <- calc_metrics(dataSurvPred_large, "Scenario 3 (High-Dim)")

all_res <- bind_rows(df_linear, df_nonlinear, df_large)
all_res$Scenario <- factor(all_res$Scenario, levels = c("Scenario 1 (Linear)", "Scenario 2 (Nonlinear)", "Scenario 3 (High-Dim)"))
all_res$Model <- factor(all_res$Model, levels = c("Boosted Trees", "Cox"))
n_max <- max(all_res$n)

# Metrics vs Sample Size (n) [Fixed Q=10]
data_row1 <- all_res %>%
  filter((Method == "Landmark" & Q == 10) | Method == "Naive") %>%
  mutate(LegendLabel = paste(Method, Model)) 

data_row2_lines <- all_res %>%
  filter(n == n_max*0.001, Method == "Landmark")

data_row2_hline <- all_res %>%
  filter(n == n_max, Method == "Naive") %>%
  group_by(Scenario, Model) %>%
  summarise(RMSE = mean(RMSE), MAPE = mean(MAPE), .groups = "drop") 


# Plotting RMSE

p1_RMSE <- ggplot(data_row1, aes(x = n, y = RMSE, color = Model, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = unique(sort(data_row1$n)), labels = scales::comma) +
  scale_linetype_manual(values = c("Landmark" = "solid", "Naive" = "dotted")) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "RMSE", x = "Sample Size (n) [Log Scale]", title = "RMSE for varying sample size (Q = 10)") +
  theme(legend.position = "top", plot.title = element_text(face = "bold"), legend.key.width = unit(1.5, "cm"))

p2_RMSE <- ggplot(data_row2_lines, aes(x = Q, y = RMSE, color = Model)) +
  geom_hline(data = data_row2_hline, aes(yintercept = RMSE, color = Model), 
             linetype = "dotted", linewidth = 0.7) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(1, 2, 5, 10)) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "RMSE", x = "Number of Landmarks (Q)", 
       title = paste0("RMSE for varying Q (n = ", scales::comma(n_max*0.001), ")")) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"), legend.key.width = unit(1.5, "cm"))

plot_RMSE <- grid.arrange(p1_RMSE, p2_RMSE, nrow = 2, heights = c(1, 1))
ggsave("Plots/RMSE.pdf", plot = plot_RMSE, width = 10, height = 6, units = "in")



# Plotting MAPE

p1_MAPE <- ggplot(data_row1, aes(x = n, y = MAPE, color = Model, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = unique(sort(data_row1$n)), labels = scales::comma) +
  scale_linetype_manual(values = c("Landmark" = "solid", "Naive" = "dotted")) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "MAPE", x = "Sample Size (n) [Log Scale]", title = "MAPE for varying sample size (Q = 10)") +
  theme(legend.position = "top", plot.title = element_text(face = "bold"), legend.key.width = unit(1.5, "cm"))

p2_MAPE <- ggplot(data_row2_lines, aes(x = Q, y = MAPE, color = Model)) +
  geom_hline(data = data_row2_hline, aes(yintercept = MAPE, color = Model), 
             linetype = "dotted", linewidth = 0.7) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = c(1, 2, 5, 10)) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "MAPE", x = "Number of Landmarks (Q)", 
       title = paste0("MAPE for varying Q (n = ", scales::comma(n_max*0.001), ")")) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"), legend.key.width = unit(1.5, "cm"))

plot_MAPE <- grid.arrange(p1_MAPE, p2_MAPE, nrow = 2, heights = c(1, 1))
ggsave("Plots/MAPE.pdf", plot = plot_MAPE, width = 10, height = 6, units = "in")



# Plotting RMSE for varying Q and remaining n

data_row1_lines <- all_res %>%
  filter(n == n_max*0.01, Method == "Landmark")

data_row1_hline <- all_res %>%
  filter(n == n_max*0.01, Method == "Naive") %>%
  group_by(Scenario, Model) %>%
  summarise(RMSE = mean(RMSE), .groups = "drop") 

p1_RMSE <- ggplot(data_row1_lines, aes(x = Q, y = RMSE, color = Model)) +
  geom_hline(data = data_row1_hline, aes(yintercept = RMSE, color = Model, linetype = "Naive"), 
             linewidth = 0.7) +
  geom_line(aes(linetype = "Landmark"), linewidth = 0.8) +
  geom_point(size = 2) +
  
  scale_x_continuous(breaks = c(1, 2, 5, 10)) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  
  scale_linetype_manual(name = "Method", values = c("Landmark" = "solid", "Naive" = "dotted")) +
  
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "RMSE", x = NULL, 
       title = paste0("RMSE for varying Q (n = ", scales::comma(n_max*0.01), ")")) +
  theme(legend.position = "top", 
        plot.title = element_text(face = "bold"), 
        legend.key.width = unit(1.5, "cm"))


data_row2_lines <- all_res %>%
  filter(n == n_max*0.1, Method == "Landmark")

data_row2_hline <- all_res %>%
  filter(n == n_max*0.1, Method == "Naive") %>%
  group_by(Scenario, Model) %>%
  summarise(RMSE = mean(RMSE), .groups = "drop") 

p2_RMSE <- ggplot(data_row2_lines, aes(x = Q, y = RMSE, color = Model)) +
  geom_hline(data = data_row2_hline, aes(yintercept = RMSE, color = Model, linetype = "Naive"), 
             linewidth = 0.7) +
  geom_line(aes(linetype = "Landmark"), linewidth = 0.8) +
  geom_point(size = 2) +
  
  scale_x_continuous(breaks = c(1, 2, 5, 10)) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  scale_linetype_manual(name = "Method", values = c("Landmark" = "solid", "Naive" = "dotted")) +
  
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "RMSE", x = NULL, 
       title = paste0("RMSE for varying Q (n = ", scales::comma(n_max*0.1), ")")) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))


data_row3_lines <- all_res %>%
  filter(n == n_max, Method == "Landmark")

data_row3_hline <- all_res %>%
  filter(n == n_max, Method == "Naive") %>%
  group_by(Scenario, Model) %>%
  summarise(RMSE = mean(RMSE), .groups = "drop") 

p3_RMSE <- ggplot(data_row3_lines, aes(x = Q, y = RMSE, color = Model)) +
  geom_hline(data = data_row3_hline, aes(yintercept = RMSE, color = Model, linetype = "Naive"), 
             linewidth = 0.7) +
  geom_line(aes(linetype = "Landmark"), linewidth = 0.8) +
  geom_point(size = 2) +
  
  scale_x_continuous(breaks = c(1, 2, 5, 10)) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  scale_linetype_manual(name = "Method", values = c("Landmark" = "solid", "Naive" = "dotted")) +
  
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "RMSE", x = "Number of Landmarks (Q)", 
       title = paste0("RMSE for varying Q (n = ", scales::comma(n_max), ")")) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

plot_RMSE <- grid.arrange(p1_RMSE, p2_RMSE, p3_RMSE, nrow = 3, heights = c(1.25, 1, 1.15))
ggsave("Plots/RMSE_Q.pdf", plot = plot_RMSE, width = 10, height = 8, units = "in")


# Plotting RMSE for varying n and remaining Q

data_row1 <- all_res %>%
  filter((Method == "Landmark" & Q == 1) | Method == "Naive") %>%
  mutate(LegendLabel = paste(Method, Model)) 

p1_RMSE <- ggplot(data_row1, aes(x = n, y = RMSE, color = Model, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = unique(sort(data_row1$n)), labels = scales::comma) +
  scale_linetype_manual(values = c("Landmark" = "solid", "Naive" = "dotted")) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "RMSE", x = NULL, title = "RMSE for varying sample size (Q = 1)") +
  theme(legend.position = "top", plot.title = element_text(face = "bold"), legend.key.width = unit(1.5, "cm"))

data_row2 <- all_res %>%
  filter((Method == "Landmark" & Q == 2) | Method == "Naive") %>%
  mutate(LegendLabel = paste(Method, Model)) 

p2_RMSE <- ggplot(data_row2, aes(x = n, y = RMSE, color = Model, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = unique(sort(data_row2$n)), labels = scales::comma) +
  scale_linetype_manual(values = c("Landmark" = "solid", "Naive" = "dotted")) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "RMSE", x = NULL, title = "RMSE for varying sample size (Q = 2)") +
  theme(legend.position = "none", plot.title = element_text(face = "bold"), legend.key.width = unit(1.5, "cm"))


data_row3 <- all_res %>%
  filter((Method == "Landmark" & Q == 5) | Method == "Naive") %>%
  mutate(LegendLabel = paste(Method, Model)) 

p3_RMSE <- ggplot(data_row3, aes(x = n, y = RMSE, color = Model, linetype = Method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_log10(breaks = unique(sort(data_row3$n)), labels = scales::comma) +
  scale_linetype_manual(values = c("Landmark" = "solid", "Naive" = "dotted")) +
  scale_color_manual(values = c("Boosted Trees" = "#000000", "Cox" = "#D55E00")) +
  facet_wrap(~Scenario, scales = "free_y", nrow = 1) +
  theme_bw() +
  labs(y = "RMSE", x = "Sample Size (n) [Log Scale]", title = "RMSE for varying sample size (Q = 5)") +
  theme(legend.position = "none", plot.title = element_text(face = "bold"), legend.key.width = unit(1.5, "cm"))

plot_RMSE <- grid.arrange(p1_RMSE, p2_RMSE, p3_RMSE, nrow = 3, heights = c(1.25, 1, 1.15))
ggsave("Plots/RMSE_n.pdf", plot = plot_RMSE, width = 10, height = 6, units = "in")
