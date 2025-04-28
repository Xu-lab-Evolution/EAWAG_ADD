# --------------------------------------------------------------------------
# R Script for Fitting ODE Models to Pond Ecosystem Data (Year 2021)
# SINGLE GLOBAL FIT with DAM Model (v19b)
# - Fits Full (DAM) model using ALL pond data combined
# - Includes time-dependent mortality ONLY for Duckweed (D)
# - Uses simpler N:P effect formulation (optimum + sensitivity)
# - Includes shading effect of Duckweed (D) on Temperature for A & M
# - Includes pH effect ONLY for Daphnia (M) growth
# - Includes density-dependent mortality for Daphnia (M)
# - Uses Quadratic temperature response for growth/consumption
# - Temp Opt/Range & Nutrient Optima parameters guided by literature
# - Includes temperature-dependent mortality for D/M (higher at low temps)
# - Rescales DuckweedCoverage if > 100% (Note: scaling handled differently now)
# - Performs Residual Analysis for diagnostics
# - Adds species time series plots facetted by Pond AND Species
# - Adds observed interaction plots (including N:P ratio and pH)
# --------------------------------------------------------------------------

# --- Load Required Libraries ---
library(deSolve)
library(dplyr)
library(nloptr)
library(ggplot2)
library(tidyr)
library(purrr)
library(zoo)
library(patchwork)

# --- Specify Current Time and Location ---
current_time <- Sys.time()
current_timezone <- Sys.timezone()
message(paste("Script run initiated at:", format(current_time, "%Y-%m-%d %H:%M:%S %Z"), "in timezone:", current_timezone))

# --- Specify Input Data File ---
# Assuming 'environment.txt' is accessible
data_file <- if (file.exists("environment.txt")) {
  "environment.txt"
} else if (file.exists("data/environment.txt")) {
  "data/environment.txt"
} else {
  stop("Could not find environment.txt in current directory or ./data/. Please specify the correct path.")
}

# --- Define Atomic Weights ---
# Used for N:P molar ratio calculation 
MW_NH4 <- 18.04   # ammonium
MW_PO4 <- 94.97   # phosphate


# -------------------------------
# Read and Prepare the Data
# -------------------------------
message("Reading data from: ", data_file)
df_raw <- tryCatch({
  read.table(data_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
             comment.char = "", check.names = FALSE, na.strings = c("NA", "NaN", ""))
}, error = function(e) {
  stop("Error reading data file: ", data_file, "\n", e$message,
       "\nPlease ensure the path is correct and the file is a tab-separated text file with a header.")
})

message("Initial data rows:")
print(head(df_raw))

message("Preparing data...")

has_treatment_column <- "Treatment" %in% colnames(df_raw)
# Assuming 'ammonium' is ug-N/L and 'phosphate' is ug-P/L
key_cols_numeric <- c("AphidDensity", "DuckweedCoverage", "ammonium", "phosphate", "chlA", "dmagna", "temperature", "pH")
key_cols_base <- c("Time", "Pond", "Year")
if (has_treatment_column) key_cols_base <- c(key_cols_base, "Treatment")
key_cols <- c(key_cols_base, key_cols_numeric)

df_filtered_year <- df_raw %>%
  filter(Year == 2021) %>%
  select(any_of(key_cols)) %>%
  mutate(across(all_of(key_cols_numeric), ~suppressWarnings(as.numeric(as.character(.)))))

# Note: DuckweedCoverage kept potentially > 100 for now, used in ODE, compared to raw obs in SSE.
df_processed_temp <- df_filtered_year

# --- Calculate N:P Molar Ratio from ug/L ---
message("Calculating N:P molar ratio from ug/L ")
df_processed_temp <- df_processed_temp %>%
  mutate(
    Phosphate_safe = pmax(phosphate, 1e-6), # Avoid division by zero (using a tiny value instead of 1e-9)
    Ammonium_umol = ammonium / MW_NH4,       # Convert ug-N/L -> umol-N/L
    Phosphate_umol = Phosphate_safe / MW_PO4, # Convert ug-P/L -> umol-P/L
    NP_ratio = Ammonium_umol / Phosphate_umol # Molar ratio N/P
  ) %>%
  select(-Phosphate_safe, -Ammonium_umol, -Phosphate_umol) # Remove intermediate cols

if("NP_ratio" %in% colnames(df_processed_temp)) {
  key_cols_numeric <- c(key_cols_numeric, "NP_ratio")
} else {
  warning("NP_ratio column could not be calculated.")
}

# --- Final Processing Steps ---
df_processed_intermediate <- df_processed_temp %>%
  filter(if_all(all_of(c("Time", key_cols_numeric)), ~!is.na(.))) %>%
  mutate(across(all_of(setdiff(key_cols_numeric, "DuckweedCoverage")), ~pmax(0, .))) %>%
  mutate(DuckweedCoverage = pmax(0, DuckweedCoverage))


# --- Handle potential duplicate Time entries per Pond/Treatment by averaging ---
grouping_vars <- c("Time", "Pond")
if (has_treatment_column && "Treatment" %in% colnames(df_processed_intermediate)) {
  grouping_vars <- c(grouping_vars, "Treatment")
}
df_processed_avg <- df_processed_intermediate %>%
  group_by(across(all_of(grouping_vars))) %>%
  summarise(across(all_of(key_cols_numeric), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  filter(if_all(all_of(c("Time", key_cols_numeric)), ~!is.na(.)))

# --- Continue processing with averaged data ---
grouping_vars_final <- "Pond"
if (has_treatment_column && "Treatment" %in% colnames(df_processed_avg)) {
  grouping_vars_final <- c("Treatment", "Pond")
}
df_processed <- df_processed_avg %>%
  group_by(across(all_of(grouping_vars_final))) %>%
  filter(n() >= 3) %>% ungroup() %>%
  arrange(across(all_of(grouping_vars_final)), Time)

if(nrow(df_processed) == 0) { stop("No valid data remaining after filtering and processing.") }
message("Data preparation complete (duplicates averaged, NP molar ratio calculated from ug/L, pH included). Rows remaining: ", nrow(df_processed))
print(head(df_processed))
print(summary(df_processed$NP_ratio)) # Check calculated NP_ratio range
has_treatment_column_final <- "Treatment" %in% colnames(df_processed)

# -------------------------------
# Helper Functions
# -------------------------------
first_record <- function(x) { val <- x[!is.na(x)][1]; if (length(val) == 0 || is.na(val)) { warning("No valid first record. Returning 0."); return(0) }; return(max(0, val)) }
make_interp_func <- function(Time_vec, val_vec, fallback_val, var_name="variable") { valid_idx <- !is.na(Time_vec) & !is.na(val_vec); n_valid <- sum(valid_idx); if (n_valid < 2) { mean_val <- if (n_valid == 1) val_vec[valid_idx] else mean(val_vec[valid_idx], na.rm = TRUE); const_val <- if(is.finite(mean_val)) mean_val else fallback_val; warning(paste0("< 2 valid points for ", var_name ,". Using constant: ", round(const_val, 3))); return(function(t) const_val) }; order_idx <- order(Time_vec[valid_idx]); time_ordered <- Time_vec[valid_idx][order_idx]; val_ordered <- val_vec[valid_idx][order_idx]; approxfun(time_ordered, val_ordered, method = "linear", rule = 2, f = 0) }
calculate_temp_effect <- function(T_val, T_min, T_opt, T_max) { if (is.na(T_val) || is.na(T_min) || is.na(T_opt) || is.na(T_max)) return(0); if (T_opt < T_min - 1e-6 || T_opt > T_max + 1e-6 || T_min > T_max + 1e-6) return(0); if (T_val < T_min - 1e-6 || T_val > T_max + 1e-6) { return(0) } else if (abs(T_val - T_opt) < 1e-6) { return(1.0) } else if (T_val >= T_min - 1e-6 && T_val < T_opt) { denominator <- T_opt - T_min; if (denominator < 1e-6) return(1.0); effect <- 1 - ((T_val - T_opt) / denominator)^2 } else { denominator <- T_max - T_opt; if (denominator < 1e-6) return(1.0); effect <- 1 - ((T_val - T_opt) / denominator)^2 }; return(max(0, effect)) }
calculate_pH_effect <- function(pH_val, pH_min, pH_opt, pH_max) { if (is.na(pH_val) || is.na(pH_min) || is.na(pH_opt) || is.na(pH_max)) return(0); if (pH_opt < pH_min - 1e-6 || pH_opt > pH_max + 1e-6 || pH_min > pH_max + 1e-6) return(0); if (pH_val < pH_min - 1e-6 || pH_val > pH_max + 1e-6) { return(0) } else if (abs(pH_val - pH_opt) < 1e-6) { return(1.0) } else if (pH_val >= pH_min - 1e-6 && pH_val < pH_opt) { denominator <- pH_opt - pH_min; if (denominator < 1e-6) return(1.0); effect <- 1 - ((pH_val - pH_opt) / denominator)^2 } else { denominator <- pH_max - pH_opt; if (denominator < 1e-6) return(1.0); effect <- 1 - ((pH_val - pH_opt) / denominator)^2 }; return(max(0, effect)) }
calculate_np_effect_simple <- function(NP_ratio_val, optNP, sensitivity) { if (is.na(NP_ratio_val) || is.na(optNP) || is.na(sensitivity) || sensitivity <= 1e-6) return(0); effect <- 1 / (1 + abs(NP_ratio_val - optNP) / sensitivity); return(max(0, effect)) }

# -------------------------------
# Function to Calculate Initial Parameter Guesses (Literature Guided Temps/Nutrients)
# -------------------------------
calculate_initial_guesses <- function(data_subset) {
  calc_growth_rate <- function(value, Time) { value <- pmax(value, 1e-6); valid_idx <- which(!is.na(value) & !is.na(Time)); if(length(valid_idx) < 2) return(0.1); ordered_idx <- order(Time[valid_idx]); value_ord <- value[valid_idx][ordered_idx]; Time_ord <- Time[valid_idx][ordered_idx]; log_diffs <- diff(log(value_ord)); time_diffs <- diff(Time_ord); valid_rate_idx <- which(time_diffs > 1e-6); if(length(valid_rate_idx) == 0) return(0.1); rates <- log_diffs[valid_rate_idx] / time_diffs[valid_rate_idx]; mean_rate <- mean(rates[is.finite(rates) & rates > 0], na.rm = TRUE); if (is.na(mean_rate) || is.nan(mean_rate) || mean_rate <= 0) return(0.1); return(min(max(mean_rate, 0.01), 2.0)) }
  r_D0_est <- calc_growth_rate(data_subset$DuckweedCoverage, data_subset$Time)
  r_A0_est <- calc_growth_rate(data_subset$chlA, data_subset$Time)
  c_M0_est <- 0.2
  
  # Literature-informed temperature parameters
  T_min_D_est <- 6; T_opt_D_est <- 27; T_max_D_est <- 36;
  T_min_A_est <- 5; T_opt_A_est <- 24; T_max_A_est <- 35;
  T_min_M_est <- 7; T_opt_M_est <- 22; T_max_M_est <- 30;
  T_crit_D_est = max(T_min_D_est, 10); T_crit_M_est = max(T_min_M_est, 10);
  
  # Literature-informed N:P Optima (Molar Ratio)
  optNP_D_est <- 16
  optNP_A_est <- 16
  
  # Simpler N:P sensitivity parameter guesses
  sensitivity_D_est <- 10
  sensitivity_A_est <- 10
  
  # pH parameters for Daphnia
  pH_min_obs <- quantile(data_subset$pH, 0.05, na.rm = TRUE, type=8); pH_max_obs <- quantile(data_subset$pH, 0.95, na.rm = TRUE, type=8); pH_mean_obs <- mean(data_subset$pH, na.rm = TRUE)
  if(is.na(pH_min_obs)) pH_min_obs=6.5; if(is.na(pH_max_obs)) pH_max_obs=8.5; if(is.na(pH_mean_obs)) pH_mean_obs=7.5;
  pH_min_M_est <- max(5.0, pH_min_obs - 0.5); pH_opt_M_est <- max(pH_min_M_est + 0.2, min(pH_max_obs + 0.2, pH_mean_obs)); pH_max_M_est <- max(pH_opt_M_est + 0.2, pH_max_obs + 0.5)
  
  # Other parameter estimates
  m_D_est <- 0.05; m_A_est <- 0.05; m_M_est <- 0.1;
  m_D_lowT_est = 0.01; m_M_lowT_est = 0.02; m_D_t_est <- 0.001
  m_M_density_est <- 0.01
  e_M_est <- 0.4; beta_DD_est <- 1/100; beta_AA_est <- 0.01; beta_DA_est <- 0.005; beta_AD_est <- 0.005; gamma_X_est <- 0.01
  k_T_shade_est <- 0.1
  
  param_names <- c("r_D0", "T_min_D", "T_opt_D", "T_max_D", "m_D",
                   "r_A0", "T_min_A", "T_opt_A", "T_max_A", "m_A",
                   "c_M0", "T_min_M", "T_opt_M", "T_max_M", "e_M", "m_M",
                   "beta_DD", "beta_DA", "beta_AA", "beta_AD",
                   "optNP_D", "sensitivity_D", "optNP_A", "sensitivity_A",
                   "gamma_X",
                   "m_D_lowT", "T_crit_D", "m_M_lowT", "T_crit_M", "m_D_t",
                   "pH_min_M", "pH_opt_M", "pH_max_M",
                   "m_M_density",
                   "k_T_shade")
  
  theta_init <- setNames(c(r_D0_est, T_min_D_est, T_opt_D_est, T_max_D_est, m_D_est,
                           r_A0_est, T_min_A_est, T_opt_A_est, T_max_A_est, m_A_est,
                           c_M0_est, T_min_M_est, T_opt_M_est, T_max_M_est, e_M_est, m_M_est,
                           beta_DD_est, beta_DA_est, beta_AA_est, beta_AD_est,
                           optNP_D_est, sensitivity_D_est, optNP_A_est, sensitivity_A_est,
                           gamma_X_est,
                           m_D_lowT_est, T_crit_D_est, m_M_lowT_est, T_crit_M_est, m_D_t_est,
                           pH_min_M_est, pH_opt_M_est, pH_max_M_est,
                           m_M_density_est,
                           k_T_shade_est),
                         param_names)
  
  # Ensure basic validity checks
  theta_init <- pmax(theta_init, 1e-7)
  if("e_M" %in% names(theta_init)) theta_init["e_M"] <- max(1e-6, min(1.0, theta_init["e_M"]))
  for (sp in c("D", "A", "M")) {
    tmin_name = paste0("T_min_", sp); topt_name = paste0("T_opt_", sp); tmax_name = paste0("T_max_", sp)
    if(all(c(tmin_name, topt_name, tmax_name) %in% names(theta_init))) {
      tmin = theta_init[tmin_name]; topt = theta_init[topt_name]; tmax = theta_init[tmax_name]
      theta_init[topt_name] <- max(tmin + 1e-6, topt); theta_init[tmax_name] <- max(theta_init[topt_name] + 1e-6, tmax)
    }
  }
  if("T_crit_D" %in% names(theta_init)) theta_init["T_crit_D"] <- max(theta_init["T_min_D"], theta_init["T_crit_D"])
  if("T_crit_M" %in% names(theta_init)) theta_init["T_crit_M"] <- max(theta_init["T_min_M"], theta_init["T_crit_M"])
  if(all(c("pH_min_M", "pH_opt_M", "pH_max_M") %in% names(theta_init))) {
    phmin = theta_init["pH_min_M"]; phopt = theta_init["pH_opt_M"]; phmax = theta_init["pH_max_M"]
    theta_init["pH_opt_M"] <- max(phmin + 1e-6, phopt); theta_init["pH_max_M"] <- max(theta_init["pH_opt_M"] + 1e-6, phmax)
  }
  if("sensitivity_D" %in% names(theta_init)) theta_init["sensitivity_D"] <- max(1e-6, theta_init["sensitivity_D"])
  if("sensitivity_A" %in% names(theta_init)) theta_init["sensitivity_A"] <- max(1e-6, theta_init["sensitivity_A"])
  if("k_T_shade" %in% names(theta_init)) theta_init["k_T_shade"] <- max(0, min(1.0, theta_init["k_T_shade"]))
  
  return(list(theta_init = theta_init, param_names = param_names))
}

# -------------------------------
# ODE Model Function
# -------------------------------
model_func_v19b <- function(t, state, parameters, temp_func, aphid_func, np_ratio_func, pH_func) {
  params_list <- as.list(parameters)
  with(as.list(c(state, params_list)), {
    # State variables
    D <- max(state["D"], 1e-9); A <- max(state["A"], 1e-9); M <- max(state["M"], 0)
    
    # Environmental drivers
    T_val_surface <- temp_func(t); if (is.na(T_val_surface)) T_val_surface <- T_opt_D
    Aphid_val <- max(0, aphid_func(t)); if (is.na(Aphid_val)) Aphid_val <- 0
    NP_ratio_val <- np_ratio_func(t); if (is.na(NP_ratio_val)) NP_ratio_val <- optNP_D # Molar ratio from interpolation
    pH_val <- pH_func(t); if (is.na(pH_val)) pH_val <- pH_opt_M
    
    # Calculate adjusted temperature due to Duckweed shading
    D_effect_shade <- min(D, 100) / 100 # Assumes D units are % coverage
    T_val_adj = T_val_surface * (1 - k_T_shade * D_effect_shade)
    T_val_D = T_val_surface; T_val_A = T_val_adj; T_val_M = T_val_adj
    
    # Temperature effects
    temp_effect_D <- calculate_temp_effect(T_val_D, T_min_D, T_opt_D, T_max_D)
    temp_effect_A <- calculate_temp_effect(T_val_A, T_min_A, T_opt_A, T_max_A)
    temp_effect_M <- calculate_temp_effect(T_val_M, T_min_M, T_opt_M, T_max_M)
    
    # Nutrient limitation (Simplified N:P ratio effect - using molar ratio)
    np_ratio_effect_D = calculate_np_effect_simple(NP_ratio_val, optNP_D, sensitivity_D)
    np_ratio_effect_A = calculate_np_effect_simple(NP_ratio_val, optNP_A, sensitivity_A)
    
    # pH effect
    pH_effect_M = calculate_pH_effect(pH_val, pH_min_M, pH_opt_M, pH_max_M)
    
    # Base rates
    r_D_base <- r_D0 * temp_effect_D * np_ratio_effect_D
    r_A_base <- r_A0 * temp_effect_A * np_ratio_effect_A
    c_M_base <- c_M0 * temp_effect_M
    
    # Effective mortality rates
    m_D_lowT_val <- if(!is.na(m_D_lowT) && m_D_lowT > 0) m_D_lowT else 0; T_crit_D_val <- if(!is.na(T_crit_D)) T_crit_D else T_min_D
    m_M_lowT_val <- if(!is.na(m_M_lowT) && m_M_lowT > 0) m_M_lowT else 0; T_crit_M_val <- if(!is.na(T_crit_M)) T_crit_M else T_min_M
    m_D_t_val <- if(!is.na(m_D_t) && m_D_t > 0) m_D_t else 0
    m_D_eff <- m_D + m_D_lowT_val * max(0, T_crit_D_val - T_val_D) + m_D_t_val * t
    m_A_eff <- m_A
    m_M_eff_base <- m_M + m_M_lowT_val * max(0, T_crit_M_val - T_val_M)
    m_D_eff <- max(1e-6, m_D_eff); m_A_eff <- max(1e-6, m_A_eff); m_M_eff_base <- max(1e-6, m_M_eff_base)
    
    # Carrying capacities
    K_D <- if(beta_DD > 1e-9) 1 / beta_DD else 1e9; K_A <- if(beta_AA > 1e-9) 1 / beta_AA else 1e9
    
    # Differential Equations
    dD <- r_D_base * D * max(0, (1 - D / K_D - beta_DA * A / K_D)) - gamma_X * Aphid_val * D - m_D_eff * D
    dA <- r_A_base * A * max(0, (1 - A / K_A - beta_AD * D / K_A)) - c_M_base * M * A - m_A_eff * A
    dM <- e_M * c_M_base * M * A * pH_effect_M - (m_M_eff_base + m_M_density * M) * M
    
    list(c(dD, dA, dM))
  })
}

# -------------------------------
# Objective Functions (Calculate Unweighted SSE + Penalties)
# -------------------------------
# (No changes needed here)
calculate_sse_core <- function(parameters, data_df, model_func_ode, state_vars, obs_map, param_penalty_func) {
  penalty <- param_penalty_func(parameters); if (penalty > 1e10) return(penalty); sse_total <- 0
  ponds_in_subset <- unique(data_df$Pond)
  for(pond_id in ponds_in_subset) {
    group <- data_df %>% filter(Pond == pond_id) %>% arrange(Time); if(nrow(group) < 3) next
    state0_values <- sapply(state_vars, function(sv) first_record(group[[ obs_map[sv] ]]))
    state0 <- setNames(pmax(state0_values, c(D=1e-6, A=1e-6, M=0)[state_vars]), state_vars)
    
    temp_func <- make_interp_func(group$Time, group$temperature, fallback_val = parameters["T_opt_D"], "Temp")
    aphid_func <- if ("gamma_X" %in% names(parameters)) make_interp_func(group$Time, group$AphidDensity, fallback_val = 0, "Aphid") else function(t) 0
    np_ratio_func <- make_interp_func(group$Time, group$NP_ratio, fallback_val = parameters["optNP_D"], "NP_Ratio") # Interpolates molar ratio
    pH_func <- make_interp_func(group$Time, group$pH, fallback_val = parameters["pH_opt_M"], "pH")
    
    Times_sim <- sort(unique(group$Time)); if(length(Times_sim) < 2) next
    out <- tryCatch({ ode(y = state0, times = Times_sim, func = model_func_ode, parms = parameters,
                          temp_func = temp_func, aphid_func = aphid_func, np_ratio_func = np_ratio_func, pH_func = pH_func,
                          method = "lsoda", rtol=1e-5, atol=1e-5, maxsteps = 10000)
    }, error = function(e) { message("ODE Error Pond ", pond_id, ": ", e$message); NULL })
    if(is.null(out) || nrow(out) != length(Times_sim) || any(!is.finite(out[,-1]))) { sse_total <- sse_total + 1e12; next }
    pred <- as.data.frame(out); colnames(pred)[1] <- "Time"; pred <- pred %>% mutate(across(all_of(state_vars), ~pmax(0, .)))
    if(any(pred[, state_vars] > 1e7)) { sse_total <- sse_total + 1e12; next }
    
    obs_group <- group %>% select(Time, all_of(unname(obs_map[state_vars])));
    pred$Time <- as.numeric(pred$Time); obs_group$Time <- as.numeric(obs_group$Time)
    merged <- merge(pred, obs_group, by = "Time", all.x = TRUE)
    
    sse_pond <- 0; for(sv in state_vars) { obs_col <- obs_map[sv]; pred_col <- sv; sse_pond <- sse_pond + sum((merged[[pred_col]] - merged[[obs_col]])^2, na.rm = TRUE) }
    sse_total <- sse_total + sse_pond
  }
  final_sse <- sse_total + penalty; if (!is.finite(final_sse)) { warning("Non-finite total SSE calculated."); return(1e12 + penalty) }; return(final_sse)
}

# Penalty function
check_common_penalties <- function(params) {
  penalty <- 0
  for (sp in c("D", "A", "M")) { tmin_name = paste0("T_min_", sp); topt_name = paste0("T_opt_", sp); tmax_name = paste0("T_max_", sp)
  if(all(c(tmin_name, topt_name, tmax_name) %in% names(params))) { T_min_p <- params[tmin_name]; T_opt_p <- params[topt_name]; T_max_p <- params[tmax_name]; if (T_opt_p < T_min_p + 1e-3 || T_opt_p > T_max_p - 1e-3 || T_min_p > T_max_p - 1e-3) { penalty <- penalty + ( (max(0, T_min_p - T_opt_p)) + (max(0, T_opt_p - T_max_p)) + (max(0, T_min_p - T_max_p)) ) * 1e6 } } }
  if(all(c("T_crit_D", "T_min_D") %in% names(params))) { if (params["T_crit_D"] < params["T_min_D"] - 1e-3) penalty <- penalty + (params["T_min_D"] - params["T_crit_D"]) * 1e6 }
  if(all(c("T_crit_M", "T_min_M") %in% names(params))) { if (params["T_crit_M"] < params["T_min_M"] - 1e-3) penalty <- penalty + (params["T_min_M"] - params["T_crit_M"]) * 1e6 }
  phmin_name = "pH_min_M"; phopt_name = "pH_opt_M"; phmax_name = "pH_max_M"
  if(all(c(phmin_name, phopt_name, phmax_name) %in% names(params))) { pH_min_p <- params[phmin_name]; pH_opt_p <- params[phopt_name]; pH_max_p <- params[phmax_name]; if (pH_opt_p < pH_min_p + 1e-3 || pH_opt_p > pH_max_p - 1e-3 || pH_min_p > pH_max_p - 1e-3) { penalty <- penalty + ( (max(0, pH_min_p - pH_opt_p)) + (max(0, pH_opt_p - pH_max_p)) + (max(0, pH_min_p - pH_max_p)) ) * 1e6 } }
  eff_params = c("e_M"); for(eff_p in eff_params) { if(eff_p %in% names(params)) { if (params[eff_p] <= 1e-6 || params[eff_p] > 1.0 + 1e-3) penalty <- penalty + abs(params[eff_p] - 0.5) * 1e6 } }
  # Ensure k_T_shade is checked against bounds (0-1)
  if("k_T_shade" %in% names(params)) { if(params["k_T_shade"] < -1e-3 || params["k_T_shade"] > 1.0 + 1e-3) penalty <- penalty + abs(params["k_T_shade"] - 0.5) * 1e7 } # Penalize outside 0-1
  return(penalty)
}

# Specific Objective Function Wrapper
objective_sse_v19b <- function(log_theta, data_df, param_names_list) {
  parameters <- setNames(exp(log_theta), param_names_list)
  if (any(!is.finite(parameters)) || any(is.na(parameters))) return(1e12)
  calculate_sse_core(parameters, data_df, model_func_v19b, # Use new ODE func name
                     state_vars = c("D", "A", "M"),
                     obs_map = c(D = "DuckweedCoverage", A = "chlA", M = "dmagna"),
                     param_penalty_func = check_common_penalties)
}

# -------------------------------
# Optimization Runner Function
# -------------------------------
run_optimization <- function(model_type, obj_func, model_func_ode, relevant_params, data_df, log_bounds_global, param_names_global) {
  message("--- Running Optimization for Model: ", model_type, " ---")
  param_init_list_full <- calculate_initial_guesses(data_df)
  theta_init_full <- param_init_list_full$theta_init
  param_names_subset <- relevant_params
  theta_init_subset <- theta_init_full[param_names_subset]
  log_theta_init_subset <- log(pmax(1e-7, theta_init_subset))
  
  if (!all(param_names_subset %in% param_names_global)) {
    missing_params <- setdiff(param_names_subset, param_names_global)
    stop("Parameters missing from global list for model ", model_type, ": ", paste(missing_params, collapse=", "))
  }
  global_param_indices <- setNames(1:length(param_names_global), param_names_global)
  subset_indices <- global_param_indices[param_names_subset]
  log_lower_bounds_subset <- log_bounds_global$lower[subset_indices]
  log_upper_bounds_subset <- log_bounds_global$upper[subset_indices]
  
  if(length(log_theta_init_subset) != length(param_names_subset) || length(log_lower_bounds_subset) != length(param_names_subset) || length(log_upper_bounds_subset) != length(param_names_subset)) {
    stop("Parameter/bound length mismatch for model ", model_type, "\n",
         "Init: ", length(log_theta_init_subset), " Names: ", length(param_names_subset),
         " Lower: ", length(log_lower_bounds_subset), " Upper: ", length(log_upper_bounds_subset))
  }
  if(any(is.na(log_theta_init_subset)) || any(is.na(log_lower_bounds_subset)) || any(is.na(log_upper_bounds_subset))) {
    na_params <- names(log_theta_init_subset)[is.na(log_theta_init_subset)]
    na_lower <- names(log_lower_bounds_subset)[is.na(log_lower_bounds_subset)]
    na_upper <- names(log_upper_bounds_subset)[is.na(log_upper_bounds_subset)]
    stop("NA values found in initial parameters or bounds for model ", model_type,
         "\n NAs in init: ", paste(na_params, collapse=", "),
         "\n NAs in lower: ", paste(na_lower, collapse=", "),
         "\n NAs in upper: ", paste(na_upper, collapse=", "))
  }
  
  opt_res <- tryCatch({
    nloptr( x0 = log_theta_init_subset, eval_f = obj_func, lb = log_lower_bounds_subset, ub = log_upper_bounds_subset,
            opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 15000, xtol_rel = 1e-5, print_level=0),
            data_df = data_df, param_names_list = param_names_subset )
  }, error = function(e) { message("Optimization error for ", model_type, ": ", e$message); NULL })
  
  if (is.null(opt_res) || !(opt_res$status %in% c(1, 2, 3, 4, 5))) {
    status_msg <- if (!is.null(opt_res)) opt_res$status else "Error"; message("Optimization may not have converged successfully for model: ", model_type, " (Status: ", status_msg, ")")
    success_flag = FALSE; if (!is.null(opt_res) && opt_res$status == 5) { message("Maxeval reached. Check objective value: ", opt_res$objective) }
    return(list(success = success_flag, results = opt_res))
  } else {
    message("Optimization successful for model: ", model_type, " (Status: ", opt_res$status, ", Objective: ", round(opt_res$objective, 4), ")")
    theta_opt <- setNames(exp(opt_res$solution), param_names_subset)
    final_penalty = check_common_penalties(theta_opt)
    if(final_penalty > 1e-3) { warning("Optimized parameters for ", model_type, " might violate constraints (Penalty: ", final_penalty, ")") }
    # Ensure k_T_shade is within [0, 1] after exp transform (bounds should handle this mostly)
    if ("k_T_shade" %in% names(theta_opt)) {
      theta_opt["k_T_shade"] <- max(0, min(1.0, theta_opt["k_T_shade"]))
    }
    return(list(success = TRUE, results = opt_res, params_opt = theta_opt, param_names = param_names_subset, model_func = model_func_ode))
  }
}

# -------------------------------
# Simulation and Residual Calculation Function
# -------------------------------
simulate_and_get_residuals <- function(fit_results, data_df, state_vars) {
  if (!fit_results$success) { model_name_msg <- if(!is.null(fit_results$model_name)) fit_results$model_name else "failed fit"; message("Skipping simulation/residuals for ", model_name_msg); return(NULL) }
  params <- fit_results$params_opt; model_func_ode <- fit_results$model_func; param_names_list <- fit_results$param_names
  all_preds <- list(); all_residuals <- list(); obs_map <- c(D = "DuckweedCoverage", A = "chlA", M = "dmagna")
  has_treatment_col_local <- "Treatment" %in% colnames(data_df)
  for(pond_id in unique(data_df$Pond)) {
    group <- data_df %>% filter(Pond == pond_id) %>% arrange(Time); if(nrow(group) < 3) next
    state0_values <- sapply(state_vars, function(sv) first_record(group[[ obs_map[sv] ]])); state0 <- setNames(pmax(state0_values, c(D=1e-6, A=1e-6, M=0)[state_vars]), state_vars)
    temp_func <- make_interp_func(group$Time, group$temperature, fallback_val = params["T_opt_D"], "Temp")
    aphid_func <- if ("gamma_X" %in% names(params)) make_interp_func(group$Time, group$AphidDensity, fallback_val = 0, "Aphid") else function(t) 0
    np_ratio_func <- make_interp_func(group$Time, group$NP_ratio, fallback_val = params["optNP_D"], "NP_Ratio")
    pH_func <- make_interp_func(group$Time, group$pH, fallback_val = params["pH_opt_M"], "pH")
    Times_obs <- sort(unique(group$Time)); if(length(Times_obs) < 2) next
    times_dense <- seq(min(Times_obs), max(Times_obs), length.out = 100); times_sim_dense <- sort(unique(c(Times_obs, times_dense)))
    out_dense <- tryCatch({ ode(y = state0, times = times_sim_dense, func = model_func_ode, parms = params, temp_func = temp_func, aphid_func = aphid_func, np_ratio_func = np_ratio_func, pH_func = pH_func, method = "lsoda", rtol=1e-6, atol=1e-6) }, error = function(e) {warning("Dense Sim Error Pond ", pond_id, ": ", e$message); NULL})
    if (is.null(out_dense) || any(!is.finite(out_dense[,-1]))) {
      warning("Dense simulation failed for pond ", pond_id)
      out_obs <- tryCatch({ ode(y = state0, times = Times_obs, func = model_func_ode, parms = params, temp_func = temp_func, aphid_func = aphid_func, np_ratio_func = np_ratio_func, pH_func = pH_func, method = "lsoda", rtol=1e-6, atol=1e-6) }, error = function(e) NULL)
      if (is.null(out_obs) || nrow(out_obs) != length(Times_obs) || any(!is.finite(out_obs[,-1]))) { warning("Residual calculation failed for pond ", pond_id); next }
      pred_df_resid <- as.data.frame(out_obs); colnames(pred_df_resid)[1] <- "Time"; pred_df_plot <- pred_df_resid
    } else { pred_df_plot <- as.data.frame(out_dense); colnames(pred_df_plot)[1] <- "Time"; pred_df_resid <- pred_df_plot %>% filter(Time %in% Times_obs) }
    pred_df_resid$Pond <- pond_id; if(has_treatment_col_local) pred_df_resid$Treatment <- group$Treatment[1]
    pred_df_plot <- pred_df_plot %>% mutate(across(all_of(state_vars), ~pmax(0, .))); pred_df_plot$Pond <- pond_id; if(has_treatment_col_local) pred_df_plot$Treatment <- group$Treatment[1]; all_preds[[length(all_preds) + 1]] <- pred_df_plot
    pred_df_resid <- pred_df_resid %>% mutate(across(all_of(state_vars), ~pmax(0, .)))
    obs_cols_needed <- unname(obs_map[state_vars]); driver_cols <- c("temperature", "phosphate", "ammonium", "AphidDensity", "NP_ratio", "pH"); join_cols <- c("Time", "Pond"); if(has_treatment_col_local) join_cols <- c(join_cols, "Treatment")
    obs_group_for_merge <- group %>% select(all_of(join_cols), all_of(obs_cols_needed), any_of(driver_cols))
    pred_df_resid <- pred_df_resid %>% mutate(Time = as.numeric(Time), Pond = as.character(Pond)); obs_group_for_merge <- obs_group_for_merge %>% mutate(Time = as.numeric(Time), Pond = as.character(Pond))
    if(has_treatment_col_local) { pred_df_resid$Treatment <- as.character(pred_df_resid$Treatment); obs_group_for_merge$Treatment <- as.character(obs_group_for_merge$Treatment) }
    merged <- dplyr::left_join(pred_df_resid, obs_group_for_merge, by = join_cols)
    resid_df <- merged
    for(sv in state_vars){ resid_col_name <- paste0("Resid_", sv); pred_col_name <- paste0("Pred_", sv); obs_col_name <- obs_map[sv]
    if(obs_col_name %in% colnames(resid_df)) { pred_val <- resid_df[[sv]]; obs_val <- resid_df[[obs_col_name]]; if (is.numeric(pred_val) && is.numeric(obs_val)) { resid_df[[resid_col_name]] <- obs_val - pred_val } else { resid_df[[resid_col_name]] <- NA }; resid_df[[pred_col_name]] <- pred_val
    } else { resid_df[[resid_col_name]] <- NA; resid_df[[pred_col_name]] <- resid_df[[sv]] } }
    select_cols <- c("Time", "Pond", "Treatment", paste0("Resid_", state_vars), paste0("Pred_", state_vars), driver_cols)
    resid_df <- resid_df %>% select(any_of(select_cols)) %>% filter(if_any(starts_with("Resid_"), ~!is.na(.)))
    all_residuals[[length(all_residuals) + 1]] <- resid_df
  }
  return(list(predictions = bind_rows(all_preds), residuals = bind_rows(all_residuals)))
}

# -------------------------------
# Diagnostic Plotting Function
# -------------------------------
plot_diagnostics <- function(residuals_df, predictions_df, data_df, model_name) {
  # ... (same as previous version) ...
  if (is.null(residuals_df) || nrow(residuals_df) == 0) { message("No residuals data to plot for ", model_name); return(NULL) }
  if (is.null(predictions_df) || nrow(predictions_df) == 0) { message("No prediction data to plot for ", model_name); return(NULL) }
  plots <- list(); obs_map <- c(D = "DuckweedCoverage", A = "chlA", M = "dmagna")
  has_treatment_col_local <- "Treatment" %in% colnames(data_df)
  obs_cols_present <- intersect(unname(obs_map), colnames(data_df))
  if(length(obs_cols_present)==0) {message("No observation columns found in data_df for ", model_name); return(NULL)}
  obs_long <- data_df %>% select(Time, Pond, any_of("Treatment"), all_of(obs_cols_present)) %>%
    pivot_longer(all_of(obs_cols_present), names_to="ObsVar", values_to="Observed") %>%
    mutate(Species = names(obs_map)[match(ObsVar, obs_map)]) %>% filter(!is.na(Species))
  state_vars_present <- intersect(c("D", "A", "M"), colnames(predictions_df))
  if(length(state_vars_present) == 0) { message("No state variables found in prediction data for ", model_name); return(NULL)}
  pred_long <- predictions_df %>% select(Time, Pond, any_of("Treatment"), all_of(state_vars_present)) %>%
    pivot_longer(all_of(state_vars_present), names_to="Species", values_to="Predicted")
  join_cols_plot <- c("Time","Pond","Species"); if(has_treatment_col_local) join_cols_plot <- c(join_cols_plot, "Treatment")
  fit_data <- full_join(obs_long, pred_long, by=join_cols_plot)
  fit_plot_aes <- if(has_treatment_col_local) aes(x=Time, y=Observed, color=Treatment) else aes(x=Time, y=Observed)
  plots$fit <- ggplot() + geom_line(data=pred_long, aes(x=Time, y=Predicted, linetype=Species), linewidth=0.7) +
    geom_point(data=obs_long, fit_plot_aes, size=1.5, alpha=0.6) + facet_wrap(~Pond, scales="free_y") +
    labs(title=paste("Fit Plot:", model_name), y="Value", x="Time (weeks)") + theme_bw(base_size=10) +
    scale_linetype_manual(values=c("D"="solid", "A"="dashed", "M"="dotted")) + theme(legend.position="top")
  if(has_treatment_col_local) plots$fit <- plots$fit + scale_color_brewer(palette="Set2")
  
  resid_cols_present <- intersect(paste0("Resid_", c("D","A","M")), colnames(residuals_df))
  resid_long <- NULL # Initialize
  if(length(resid_cols_present) > 0) {
    resid_long_tmp <- residuals_df %>% select(Time, Pond, any_of("Treatment"), all_of(resid_cols_present)) %>%
      pivot_longer(all_of(resid_cols_present), names_to="Species", values_to="Residual", names_prefix="Resid_") %>% filter(!is.na(Residual))
    if(nrow(resid_long_tmp) > 0) {
      resid_long <- resid_long_tmp # Assign if valid rows exist
      plots$res_time <- ggplot(resid_long, aes(x=Time, y=Residual)) + geom_point(alpha=0.6) + geom_hline(yintercept=0, linetype="dashed", color="red") + geom_smooth(method="loess", se=FALSE, color="blue", linewidth=0.5) + facet_grid(Species ~ Pond, scales="free_y") + labs(title="Residuals vs Time", x="Time (weeks)") + theme_bw(base_size=10)
    } else { plots$res_time <- NULL }
  } else { plots$res_time <- NULL }
  
  pred_cols_present <- intersect(paste0("Pred_", c("D","A","M")), colnames(residuals_df))
  if(length(resid_cols_present) == 0 || length(pred_cols_present) == 0) { plots$res_pred <- NULL } else {
    resid_pred_long <- residuals_df %>% select(Time, Pond, any_of("Treatment"), all_of(resid_cols_present), all_of(pred_cols_present)) %>%
      pivot_longer(all_of(resid_cols_present), names_to="Species", values_to="Residual", names_prefix="Resid_") %>%
      pivot_longer(all_of(pred_cols_present), names_to="Species_P", values_to="Predicted", names_prefix="Pred_") %>%
      filter(Species == Species_P) %>% filter(!is.na(Residual) & !is.na(Predicted))
    if(nrow(resid_pred_long) > 0) { plots$res_pred <- ggplot(resid_pred_long, aes(x=Predicted, y=Residual)) + geom_point(alpha=0.6) + geom_hline(yintercept=0, linetype="dashed", color="red") + geom_smooth(method="loess", se=FALSE, color="blue", linewidth=0.5) + facet_wrap(~Species, scales="free") + labs(title="Residuals vs Predicted") + theme_bw(base_size=10)
    } else { plots$res_pred <- NULL } }
  
  # Residuals vs Temperature
  if(!is.null(resid_long) && "temperature" %in% colnames(residuals_df)) {
    join_cols_res <- c("Time", "Pond"); if(has_treatment_col_local) join_cols_res <- c(join_cols_res, "Treatment")
    resid_temp_df <- resid_long %>% left_join(residuals_df %>% select(all_of(join_cols_res), temperature) %>% distinct(), by = join_cols_res)
    plots$res_temp <- ggplot(resid_temp_df, aes(x=temperature, y=Residual)) + geom_point(alpha=0.6) + geom_hline(yintercept=0, linetype="dashed", color="red") + geom_smooth(method="loess", se=FALSE, color="blue", linewidth=0.5) + facet_wrap(~Species, scales="free_y") + labs(title="Residuals vs Temperature", x="Temperature (C)") + theme_bw(base_size=10)
  } else { plots$res_temp <- NULL }
  
  # Residuals vs pH
  if(!is.null(resid_long) && "pH" %in% colnames(residuals_df)) {
    join_cols_res <- c("Time", "Pond"); if(has_treatment_col_local) join_cols_res <- c(join_cols_res, "Treatment")
    if ("pH" %in% colnames(residuals_df)) {
      resid_pH_df <- resid_long %>% left_join(residuals_df %>% select(all_of(join_cols_res), pH) %>% distinct(), by = join_cols_res)
      if(nrow(resid_pH_df %>% filter(!is.na(pH))) > 0) { # Check if pH data exists after join
        plots$res_pH <- ggplot(resid_pH_df, aes(x=pH, y=Residual)) +
          geom_point(alpha=0.6) +
          geom_hline(yintercept=0, linetype="dashed", color="red") +
          geom_smooth(method="loess", se=FALSE, color="blue", linewidth=0.5) +
          facet_wrap(~Species, scales="free_y") +
          labs(title="Residuals vs pH", x="pH") + theme_bw(base_size=10)
      } else { plots$res_pH <- NULL; message("Note: No non-NA pH values found for residual plotting.") }
    } else { plots$res_pH <- NULL; message("Warning: pH column missing in residuals_df for plotting.") }
  } else { plots$res_pH <- NULL }
  
  # Observed vs Predicted
  if(nrow(fit_data %>% filter(!is.na(Observed) & !is.na(Predicted))) > 0) {
    plots$obs_pred <- ggplot(fit_data %>% filter(!is.na(Observed) & !is.na(Predicted)), aes(x=Observed, y=Predicted)) + geom_point(alpha=0.6) + geom_abline(intercept=0, slope=1, color="red", linetype="dashed") + facet_wrap(~Species, scales="free") + labs(title="Observed vs Predicted") + theme_bw(base_size=10)
  } else { plots$obs_pred <- NULL }
  
  plots <- plots[!sapply(plots, is.null)]; return(plots)
}

# -------------------------------
# Species Time Series Plotting Function
# -------------------------------
plot_species_timeseries <- function(predictions_df, data_df, model_name) {
  if (is.null(predictions_df) || nrow(predictions_df) == 0) { message("No prediction data to plot species time series for ", model_name); return(NULL) }
  obs_map <- c(D = "DuckweedCoverage", A = "chlA", M = "dmagna")
  has_treatment_col_local <- "Treatment" %in% colnames(data_df)
  obs_cols_present <- intersect(unname(obs_map), colnames(data_df))
  if(length(obs_cols_present)==0) {message("No observation columns found for species time series plot"); return(NULL)}
  obs_long <- data_df %>% select(Time, Pond, any_of("Treatment"), all_of(obs_cols_present)) %>%
    pivot_longer(all_of(obs_cols_present), names_to="ObsVar", values_to="Observed") %>%
    mutate(Species = names(obs_map)[match(ObsVar, obs_map)]) %>% filter(!is.na(Species))
  state_vars_present <- intersect(c("D", "A", "M"), colnames(predictions_df))
  if(length(state_vars_present) == 0) { message("No state variables found in prediction data for species time series plot"); return(NULL)}
  pred_long <- predictions_df %>% select(Time, Pond, any_of("Treatment"), all_of(state_vars_present)) %>%
    pivot_longer(all_of(state_vars_present), names_to="Species", values_to="Predicted")
  p <- ggplot() + geom_line(data=pred_long, aes(x=Time, y=Predicted, group=1), linewidth=0.7, color="blue") +
    geom_point(data=obs_long, aes(x=Time, y=Observed, group=1), size=1.5, alpha=0.6, color="red") +
    facet_grid(Species ~ Pond, scales="free_y") + labs(title=paste("Species Time Series by Pond:", model_name), y="Value", x="Time (weeks)") +
    theme_bw(base_size=9) + theme(legend.position="none", strip.text = element_text(face="bold", size=7), axis.text.x = element_text(angle = 45, hjust = 1, size=7))
  return(p)
}

# -------------------------------
# Observed Interaction Plotting Function
# -------------------------------
# (No changes needed here)
plot_observed_interactions <- function(data_df, model_name) {
  # ... (same as previous version) ...
  message("--- Generating Observed Interaction Plots ---")
  plots <- list(); has_treatment_col_local <- "Treatment" %in% colnames(data_df); plot_aes <- if(has_treatment_col_local) aes(color=Treatment) else aes()
  create_scatter <- function(data, x_var, y_var, x_lab=x_var, y_lab=y_var) { if (!all(c(x_var, y_var) %in% colnames(data))) return(ggplot() + labs(title=paste(y_lab, "vs", x_lab, "(Missing Data)")) + theme_void()); valid_data <- data %>% filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]])); if (nrow(valid_data) < 5) return(ggplot() + labs(title=paste(y_lab, "vs", x_lab, "(Insufficient Data)")) + theme_void()); ggplot(valid_data, aes(x=.data[[x_var]], y=.data[[y_var]])) + geom_point(plot_aes, alpha=0.6, size=1.5) + geom_smooth(method="loess", se=FALSE, color="blue", linewidth=0.8) + labs(x=x_lab, y=y_lab) + theme_bw(base_size=9) + theme(legend.position = "none") }
  plots$M_vs_A <- create_scatter(data_df, "chlA", "dmagna", x_lab="Algae (chlA)", y_lab="Daphnia"); plots$D_vs_A <- create_scatter(data_df, "chlA", "DuckweedCoverage", x_lab="Algae (chlA)", y_lab="Duckweed (%)"); plots$Aphid_vs_D <- create_scatter(data_df, "DuckweedCoverage", "AphidDensity", x_lab="Duckweed (%)", y_lab="Aphid Density")
  if("NP_ratio" %in% colnames(data_df)){ plots$NPratio_vs_D <- create_scatter(data_df, "NP_ratio", "DuckweedCoverage", x_lab="N:P Ratio (molar)", y_lab="Duckweed (%)"); plots$NPratio_vs_A <- create_scatter(data_df, "NP_ratio", "chlA", x_lab="N:P Ratio (molar)", y_lab="Algae (chlA)") } else { plots$Phos_vs_D <- create_scatter(data_df, "phosphate", "DuckweedCoverage", x_lab="Phosphate", y_lab="Duckweed (%)"); plots$Ammo_vs_D <- create_scatter(data_df, "ammonium", "DuckweedCoverage", x_lab="Ammonium", y_lab="Duckweed (%)"); plots$Phos_vs_A <- create_scatter(data_df, "phosphate", "chlA", x_lab="Phosphate", y_lab="Algae (chlA)"); plots$Ammo_vs_A <- create_scatter(data_df, "ammonium", "chlA", x_lab="Ammonium", y_lab="Algae (chlA)") }
  plots$Temp_vs_D <- create_scatter(data_df, "temperature", "DuckweedCoverage", x_lab="Temperature (C)", y_lab="Duckweed (%)"); plots$Temp_vs_A <- create_scatter(data_df, "temperature", "chlA", x_lab="Temperature (C)", y_lab="Algae (chlA)"); plots$Temp_vs_M <- create_scatter(data_df, "temperature", "dmagna", x_lab="Temperature (C)", y_lab="Daphnia")
  plots$pH_vs_D <- create_scatter(data_df, "pH", "DuckweedCoverage", x_lab="pH", y_lab="Duckweed (%)"); plots$pH_vs_A <- create_scatter(data_df, "pH", "chlA", x_lab="pH", y_lab="Algae (chlA)"); plots$pH_vs_M <- create_scatter(data_df, "pH", "dmagna", x_lab="pH", y_lab="Daphnia")
  combined_plot <- wrap_plots(plots, ncol=4) + plot_annotation(title = paste("Observed Variable Relationships:", model_name));
  return(combined_plot)
}

# -------------------------------
# Main Workflow (v19b)
# -------------------------------

message("\n\n===== Processing All Data Together")

# Define parameters for the target model
param_init_global_list <- calculate_initial_guesses(df_processed)
params_v19b <- param_init_global_list$param_names

# Store results for the single global fit
fit_summary_global <- list()
residual_plots_global <- list()
species_ts_plot_global <- NULL
observed_interactions_plot_global <- NULL

# --- Define Global Bounds (v19b) ---
# (Bounds defined in previous version are still suitable)
param_names_global <- param_init_global_list$param_names # 35 params
n_params <- length(param_names_global)
log_lower_bounds_global <- rep(log(1e-7), n_params); log_upper_bounds_global <- rep(log(1e+6), n_params)
param_indices <- setNames(1:n_params, param_names_global)

# Temperature bounds
if("T_min_D" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_min_D"]] <- log(4); log_upper_bounds_global[param_indices["T_min_D"]] <- log(10)
if("T_opt_D" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_opt_D"]] <- log(23); log_upper_bounds_global[param_indices["T_opt_D"]] <- log(32)
if("T_max_D" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_max_D"]] <- log(33); log_upper_bounds_global[param_indices["T_max_D"]] <- log(40)
if("T_min_A" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_min_A"]] <- log(3); log_upper_bounds_global[param_indices["T_min_A"]] <- log(10)
if("T_opt_A" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_opt_A"]] <- log(18); log_upper_bounds_global[param_indices["T_opt_A"]] <- log(31)
if("T_max_A" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_max_A"]] <- log(32); log_upper_bounds_global[param_indices["T_max_A"]] <- log(40)
if("T_min_M" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_min_M"]] <- log(4); log_upper_bounds_global[param_indices["T_min_M"]] <- log(12)
if("T_opt_M" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_opt_M"]] <- log(18); log_upper_bounds_global[param_indices["T_opt_M"]] <- log(27)
if("T_max_M" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_max_M"]] <- log(28); log_upper_bounds_global[param_indices["T_max_M"]] <- log(35)
if("T_crit_D" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_crit_D"]] <- log(5); log_upper_bounds_global[param_indices["T_crit_D"]] <- log(15)
if("T_crit_M" %in% names(param_indices)) log_lower_bounds_global[param_indices["T_crit_M"]] <- log(6); log_upper_bounds_global[param_indices["T_crit_M"]] <- log(18)

# N:P Optima bounds
if("optNP_D" %in% names(param_indices)) log_lower_bounds_global[param_indices["optNP_D"]] <- log(10); log_upper_bounds_global[param_indices["optNP_D"]] <- log(30)
if("optNP_A" %in% names(param_indices)) log_lower_bounds_global[param_indices["optNP_A"]] <- log(10); log_upper_bounds_global[param_indices["optNP_A"]] <- log(30)
# N:P Sensitivity bounds
if("sensitivity_D" %in% names(param_indices)) log_lower_bounds_global[param_indices["sensitivity_D"]] <- log(2); log_upper_bounds_global[param_indices["sensitivity_D"]] <- log(50)
if("sensitivity_A" %in% names(param_indices)) log_lower_bounds_global[param_indices["sensitivity_A"]] <- log(2); log_upper_bounds_global[param_indices["sensitivity_A"]] <- log(50)

# pH bounds
pH_M_params <- c("pH_min_M", "pH_opt_M", "pH_max_M"); if(all(pH_M_params %in% names(param_indices))) { log_lower_bounds_global[param_indices[pH_M_params]] <- log(5.0); log_upper_bounds_global[param_indices[pH_M_params]] <- log(10.0) }

# Shading effect bounds
if("k_T_shade" %in% names(param_indices)) { log_lower_bounds_global[param_indices["k_T_shade"]] <- log(1e-3); log_upper_bounds_global[param_indices["k_T_shade"]] <- log(0.9) }

# Other bounds
if("e_M" %in% names(param_indices)) { log_lower_bounds_global[param_indices["e_M"]] <- log(1e-6); log_upper_bounds_global[param_indices["e_M"]] <- log(1.0) }
m_params <- grep("^m_[DAM]$", param_names_global, value = TRUE); if(length(m_params)>0) { log_lower_bounds_global[param_indices[m_params]] <- log(1e-8); log_upper_bounds_global[param_indices[m_params]] <- log(5.0) }
m_lowT_params <- grep("^m_._lowT$", param_names_global, value = TRUE); if(length(m_lowT_params)>0) { log_lower_bounds_global[param_indices[m_lowT_params]] <- log(1e-8); log_upper_bounds_global[param_indices[m_lowT_params]] <- log(1.0) }
m_t_params <- grep("^m_D_t$", param_names_global, value = TRUE); if(length(m_t_params)>0) { log_lower_bounds_global[param_indices[m_t_params]] <- log(1e-7); log_upper_bounds_global[param_indices[m_t_params]] <- log(0.05) }
if("m_M_density" %in% names(param_indices)) { log_lower_bounds_global[param_indices["m_M_density"]] <- log(1e-7); log_upper_bounds_global[param_indices["m_M_density"]] <- log(1.0) }
if("r_D0" %in% names(param_indices)) log_upper_bounds_global[param_indices["r_D0"]] <- log(20.0)
if("r_A0" %in% names(param_indices)) log_upper_bounds_global[param_indices["r_A0"]] <- log(20.0)
if("c_M0" %in% names(param_indices)) log_upper_bounds_global[param_indices["c_M0"]] <- log(5.0)
if("beta_DD" %in% names(param_indices)) { log_lower_bounds_global[param_indices["beta_DD"]] <- log(1/150); log_upper_bounds_global[param_indices["beta_DD"]] <- log(1/50) }
if("beta_AA" %in% names(param_indices)) { log_lower_bounds_global[param_indices["beta_AA"]] <- log(1e-4); log_upper_bounds_global[param_indices["beta_AA"]] <- log(0.1) }
if("beta_DA" %in% names(param_indices)) { log_lower_bounds_global[param_indices["beta_DA"]] <- log(1e-5); log_upper_bounds_global[param_indices["beta_DA"]] <- log(0.1) }
if("beta_AD" %in% names(param_indices)) { log_lower_bounds_global[param_indices["beta_AD"]] <- log(1e-5); log_upper_bounds_global[param_indices["beta_AD"]] <- log(0.1) }
if("gamma_X" %in% names(param_indices)) { log_lower_bounds_global[param_indices["gamma_X"]] <- log(1e-4); log_upper_bounds_global[param_indices["gamma_X"]] <- log(0.5) }

log_bounds_global <- list(lower = log_lower_bounds_global, upper = log_upper_bounds_global)

# --- Define the v19b Model to Run ---
model_info <- list(name="v19b_CorrectNP_LitParams_Shade_SimpleNP_MpH_MDDM",
                   obj=objective_sse_v19b, # Use v19b objective
                   ode=model_func_v19b,    # Use v19b ODE func
                   params=params_v19b,     # Use v19b param list
                   states=c("D", "A", "M"))

# --- Run Optimization and Diagnostics ---
model_name <- model_info$name
fit_result <- run_optimization(model_name, model_info$obj, model_info$ode, model_info$params,
                               df_processed, log_bounds_global, param_names_global)
fit_result$model_name <- model_name
fit_summary_global[[model_name]] <- fit_result

if(fit_result$success) {
  sim_resid_data <- simulate_and_get_residuals(fit_result, df_processed, state_vars = model_info$states)
  if (!is.null(sim_resid_data)) {
    diag_plots <- plot_diagnostics(sim_resid_data$residuals, sim_resid_data$predictions, df_processed, paste(model_name, "(All Data)"))
    residual_plots_global[[model_name]] <- diag_plots
    if(!is.null(diag_plots$fit)) print(diag_plots$fit)
    
    species_ts_plot_global <- plot_species_timeseries(sim_resid_data$predictions, df_processed, paste(model_name, "(All Data)"))
    if(!is.null(species_ts_plot_global)) print(species_ts_plot_global)
    
    observed_interactions_plot_global <- plot_observed_interactions(df_processed, paste(model_name, "(All Data)"))
    if(!is.null(observed_interactions_plot_global)) print(observed_interactions_plot_global)
    
  } else {
    message("Simulation/Residual calculation failed for ", model_name)
    residual_plots_global[[model_name]] <- NULL; species_ts_plot_global <- NULL; observed_interactions_plot_global <- NULL
  }
} else {
  residual_plots_global[[model_name]] <- NULL; species_ts_plot_global <- NULL; observed_interactions_plot_global <- NULL
}

# -------------------------------
# Display All Diagnostic Plots (v19b)
# -------------------------------
message("\n\n===== Displaying Diagnostic Plots for Global DAM Fit (v19b) =====")
if(length(residual_plots_global) > 0) {
  model_name_iter <- names(residual_plots_global)[1] # Use different variable name
  diag_plots <- residual_plots_global[[model_name_iter]]
  if (!is.null(diag_plots) && length(diag_plots) > 0) {
    message("--- Diagnostics for Model: ", model_name_iter, " (All Data) ---")
    p_fit <- if(!is.null(diag_plots$fit)) diag_plots$fit else ggplot() + labs(title="Fit Plot N/A") + theme_void()
    p_obs_pred <- if(!is.null(diag_plots$obs_pred)) diag_plots$obs_pred else ggplot() + labs(title="Obs/Pred N/A") + theme_void()
    p_res_time <- if(!is.null(diag_plots$res_time)) diag_plots$res_time else ggplot() + labs(title="Resid/Time N/A") + theme_void()
    p_res_pred <- if(!is.null(diag_plots$res_pred)) diag_plots$res_pred else ggplot() + labs(title="Resid/Pred N/A") + theme_void()
    p_res_temp <- if(!is.null(diag_plots$res_temp)) diag_plots$res_temp else ggplot() + labs(title="Resid/Temp N/A") + theme_void()
    p_res_pH <- if(!is.null(diag_plots$res_pH)) diag_plots$res_pH else ggplot() + labs(title="Resid/pH N/A") + theme_void()
    
    combined_diag_plot <- (p_fit + p_obs_pred) / (p_res_time) / (p_res_pred + p_res_temp + p_res_pH) +
      plot_annotation(title = paste("Diagnostic Plots (", model_name_iter, " - All Data - Unweighted)")) &
      theme(plot.title = element_text(size=10), axis.title = element_text(size=9))
    
    tryCatch({ print(combined_diag_plot) }, error = function(e) { message("Error displaying combined plot for ", model_name_iter, ": ", e$message)})
    
    # Optional: Save the plot
    # file_name <- paste0("diagnostics_", gsub("[^A-Za-z0-9_]", "_", model_name_iter), "_GlobalFit.png")
    # ggsave(filename = file_name, plot=combined_diag_plot, width=16, height=18, units="in", dpi=150)
  } else { message("--- No diagnostic plots available for Model: ", model_name_iter, " ---") }
} else { message("No diagnostic plots generated for the global fit.") }

# Display species time series plot
if(!is.null(species_ts_plot_global)){
  message("\n\n===== Displaying Species Time Series Plot (Faceted by Pond AND Species) =====")
  tryCatch({ print(species_ts_plot_global) }, error = function(e) { message("Error displaying species time series plot: ", e$message)})
  # Optional: Save the plot
  # ggsave(filename = paste0("species_timeseries_", gsub("[^A-Za-z0-9_]", "_", model_name_iter), "_GlobalFit.png"), plot=species_ts_plot_global, width=16, height=12, units="in", dpi=150)
}

# Display observed interactions plot
if(!is.null(observed_interactions_plot_global)){
  message("\n\n===== Displaying Observed Interaction Plots (including pH) =====")
  tryCatch({ print(observed_interactions_plot_global) }, error = function(e) { message("Error displaying observed interaction plot: ", e$message)})
  # Optional: Save the plot
  # ggsave(filename = paste0("observed_interactions_", gsub("[^A-Za-z0-9_]", "_", model_name_iter), "_GlobalFit.png"), plot=observed_interactions_plot_global, width=12, height=14, units="in", dpi=150)
}

# -------------------------------
# Final Parameter Summary Output
# -------------------------------
message("\n--- Global DAM Fit (v19b) Phase Complete ---")
message("Review the console output for optimization success/failure messages.")
message("Examine the generated plots (if any) for residual patterns, particularly for Daphnia (M).")
message("Inspect the 'fit_summary_global' list object for parameter estimates and convergence details.")

message("\nOptimized Parameter Summary (Global DAM Fit - v19b - Unweighted SSE):")
for(model_name_iter in names(fit_summary_global)) { # Use different loop variable
  if(fit_summary_global[[model_name_iter]]$success) {
    message("\n--- Model: ", model_name_iter, " ---")
    opt_params <- fit_summary_global[[model_name_iter]]$params_opt
    print(format(opt_params, digits=4, scientific=FALSE))
  } else {
    message("\n--- Model: ", model_name_iter, " (Fit Failed) ---")
  }
}

# --- End of Script ---