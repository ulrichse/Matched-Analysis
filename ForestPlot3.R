
## Forest plot for multiple lags (not cumulative)

library(grid)
library(forestploter)
library(data.table)
library(tidyr)

forest_data <- read.csv("Data/Forest_Plot_Results.csv")

forest_data <- forest_data %>%
  filter(Outcome=="SMI")%>%
  filter(Model=="GLM")%>%
  dplyr::select(Subgroup, Type, lag0, lag1, lag2, lag3, lag4, lag5, lag6, lag7)

forest_data <- pivot_wider(forest_data, names_from = Type, values_from = c("lag0", "lag1", "lag2", "lag3", "lag4", "lag5", "lag6", "lag7"))

dt <- forest_data

dt$`    lag0` <- paste(rep(" ", 20), collapse = " ")
dt$` lag1` <- paste(rep(" ", 20), collapse = " ")
dt$`    lag2` <- paste(rep(" ", 20), collapse = " ")
dt$` lag3` <- paste(rep(" ", 20), collapse = " ")
dt$`    lag4` <- paste(rep(" ", 20), collapse = " ")
dt$` lag5` <- paste(rep(" ", 20), collapse = " ")
dt$`    lag6` <- paste(rep(" ", 20), collapse = " ")
dt$` lag7` <- paste(rep(" ", 20), collapse = " ")

dt$'lag0 95% CI' <- paste(sprintf("%.1f (%.1f, %.1f)", dt$lag0_matRRfit, dt$lag0_matRRlow, dt$lag0_matRRhigh))
dt$'lag1 95% CI' <- paste(sprintf("%.1f (%.1f, %.1f)", dt$lag1_matRRfit, dt$lag1_matRRlow, dt$lag1_matRRhigh))
dt$'lag2 95% CI' <- paste(sprintf("%.1f (%.1f, %.1f)", dt$lag2_matRRfit, dt$lag2_matRRlow, dt$lag2_matRRhigh))
dt$'lag3 95% CI' <- paste(sprintf("%.1f (%.1f, %.1f)", dt$lag3_matRRfit, dt$lag3_matRRlow, dt$lag3_matRRhigh))
dt$'lag4 95% CI' <- paste(sprintf("%.1f (%.1f, %.1f)", dt$lag4_matRRfit, dt$lag4_matRRlow, dt$lag4_matRRhigh))
dt$'lag5 95% CI' <- paste(sprintf("%.1f (%.1f, %.1f)", dt$lag5_matRRfit, dt$lag5_matRRlow, dt$lag5_matRRhigh))
dt$'lag6 95% CI' <- paste(sprintf("%.1f (%.1f, %.1f)", dt$lag6_matRRfit, dt$lag6_matRRlow, dt$lag6_matRRhigh))
dt$'lag7 95% CI' <- paste(sprintf("%.1f (%.1f, %.1f)", dt$lag7_matRRfit, dt$lag7_matRRlow, dt$lag7_matRRhigh))

est = list(dt$lag0_matRRfit,
           dt$lag1_matRRfit,
           dt$lag2_matRRfit,
           dt$lag3_matRRfit,
           dt$lag4_matRRfit,
           dt$lag5_matRRfit,
           dt$lag6_matRRfit,
           dt$lag7_matRRfit)
lower = list(dt$lag0_matRRlow,
             dt$lag1_matRRlow,
             dt$lag2_matRRlow,
             dt$lag3_matRRlow,
             dt$lag4_matRRlow,
             dt$lag5_matRRlow,
             dt$lag6_matRRlow,
             dt$lag7_matRRlow)
upper = list(dt$lag0_matRRhigh,
             dt$lag1_matRRhigh,
             dt$lag2_matRRhigh,
             dt$lag3_matRRhigh,
             dt$lag4_matRRhigh,
             dt$lag5_matRRhigh,
             dt$lag6_matRRhigh,
             dt$lag7_matRRhigh)

#dt <- dt %>%
 # dplyr::select(-c(2:25))

dt <- dt %>%
  dplyr::select(Subgroup, '    lag0', 'lag0 95% CI', ' lag1', 'lag1 95% CI', '    lag2', 'lag2 95% CI', ' lag3', 'lag3 95% CI',
                '    lag4', 'lag4 95% CI', ' lag5', 'lag5 95% CI', '    lag6', 'lag6 95% CI', ' lag7', 'lag7 95% CI')

tm <- forest_theme(base_size = 10,
                   refline_lty = "solid",
                   ci_pch = c(15, 18),
                   ci_col = c("#377eb8", "#4daf4a"),
                   footnote_gp = gpar(col = "blue"),
                   legend_name = "Group",
                   legend_value = c("Trt 1", "Trt 2"),
                   vertline_lty = c("dashed", "dotted"),
                   vertline_col = c("#d6604d", "#bababa"),
                   # Table cell padding, width 4 and heights 3
                   core = list(padding = unit(c(4, 3), "mm")))

p <- forest(dt,
            est = est,
            lower = lower, 
            upper = upper,
            ci_column = c(2, 4, 6, 8, 10, 12, 14, 16),
            ref_line = 1,
            xlim = c(0.5, 1.5),
            ticks_at = c(0.5, 1, 1.5),
            title = "SMI",
            theme = tm)

plot(p)  # Replace with your plot code

