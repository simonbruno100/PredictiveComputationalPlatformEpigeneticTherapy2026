

#################

graphics.off()  # Close all open graphics devices
rm(list = ls())     # clear all
cat("\014")     # clc


##################  


# Load necessary libraries
library(deSolve)
library(ggplot2)
library(dplyr)
library(rstan)
library(foreach)
library(doParallel)
library(tidyverse)
library(readxl)
library(viridis)
library(DiceKriging)
library(DiceOptim)
library(GA)  # Genetic Algorithm package
library(future.apply)
plan(multisession)   # portable parallelization


set.seed(1234) # Ensure reproducibility

# Output location

output_location <- "/Users/..." # to be modified
output_location_figures <- "/Users/..." # to be modified
if(!file.exists(output_location)) dir.create(output_location)
if(!file.exists(output_location_figures)) dir.create(output_location_figures)




# ========================================
# k_values
# ========================================

k_values <- c(
  3.799507e-01, 8.717916e+00, 1.173375e+00, 1.203015e+01, 4.520690e-02,
  1.272923e+00, 1.108645e-01, 9.394812e+00, 1.566855e-04, 4.560950e-01,
  1.952322e+01, 1.333967e+00, 1.256836e-01, 3.627773e-01, 4.110200e-03,
  1.227741e-02, 1.603587e-01, 1.003410e+01, 1.619134e+00, 1.640542e-05,
  9.486298e-02, 4.403414e-01, 2.477224e-01, 1.341020e+01, 1.360458e+00,
  1.336714e+00, 1.655805e-01, 2.233409e+00, 1.002968e+00, 1.147803e-01,
  0.01, # a1
  0.005, # b1
  0.01, # a2
  0.005  # b2
)

# Normal: 0,0,0,0
# high res:0.01,0.005,0.01,0.005
# low res:0.01,0.05,0.01,0.05

# Convert to 1-row parameter table
params_samples <- as.data.frame(t(k_values))
colnames(params_samples) <- paste0("k", 1:34)


################################################################################

# Define the pharmacokinetic model function

define_PK_EZHi_fun <- function(D0_e, dose_times, total_time, times_between_doses) {
  
  
  pk_EZHi_model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      # Differential equations for the one-compartment model with first-order absorption
      dA_GI <- -ka * A_GI  # Amount in gut decreases due to absorption
      dC_plasma <- (ka * A_GI) / Vd - ke * C_plasma  # Change in concentration in plasma
      
      list(c(dA_GI, dC_plasma))
    })
  }
  
  # Parameters
  parameters <- c(
    ka = 1.49,   # Absorption rate constant (1/hour) (this is calculated, otherwise 2.4 from the other paper)
    ke = 0.223, # Elimination rate constant (1/hour)
    Vd = 1230   # Volume of distribution (L)
  )
  
  
  # Initial conditions
  initial_conditions <- c(
    A_GI = D0_e[1]*0.33,   # Initial amount of drug in GI tract (mg)* Bioavailability
    C_plasma = 0  # Initial concentration in plasma (mg/L)
  )
  
  # Initialize output storage
  output_df_EZHi <- data.frame()
  
  # Current state of the system
  current_state <- initial_conditions
  
  # Current time to start the loop
  current_time <- 0
  
  
  #index for dose
  iijjii <- 2
  
  # Loop over each dosing time
  for (dose_time in dose_times) {
    # Define time sequence for current loop
    time_segment <- seq(current_time, dose_time, by = times_between_doses)
    
    # Simulate the model for the current segment
    output_segment <- ode(
      func = pk_EZHi_model,
      times = time_segment,
      y = current_state,
      parms = parameters,
      method = "bdf",
      maxsteps = 5000  # Increased maximum number of steps
    )
    
    # Convert output to a data frame and store it
    output_df_EZHi <- rbind(output_df_EZHi, as.data.frame(output_segment))
    
    # Update current state to the last row of the output segment
    current_state_first <- tail(output_segment, 1)[-1]  # Remove time column
    
    # Add dose to the gut compartment
    current_state_first[1] <- current_state_first[1] + D0_e[iijjii]*0.33 #bioavailability
    
    current_state <- c(
      A_GI = current_state_first[1],   # Initial amount of drug in GI tract (mg)
      C_plasma = current_state_first[2]  # Initial concentration in plasma (mg/L)
    )
    
    # Update current time
    current_time <- dose_time
    
    
    #Update index for dose
    iijjii <- iijjii + 1
  }
  
  # Simulate the final segment after the last dose until the end of the total time
  time_segment <- seq(current_time, total_time, by = times_between_doses)
  output_segment <- ode(
    func = pk_EZHi_model,
    times = time_segment,
    y = current_state,
    parms = parameters,
    method = "bdf",
    maxsteps = 5000  # Increased maximum number of steps
  )
  
  
  output_df_EZHi <- rbind(output_df_EZHi, as.data.frame(output_segment))
  
  #TO PASS FROM mg/L to micromolar
  Me <- 653.66
  output_df_EZHi <- output_df_EZHi %>%
    mutate(C_plasma = (C_plasma * 10^6) / Me)
  
  
  
  # Return the data frame
  return(output_df_EZHi)
  
}



################################################################################
################################################################################




##### AKTi PK FUNCTION
define_PK_AKTi_fun <- function(D0_a, dose_times, total_time, times_between_doses) {
  
  
  # ============================================================
  # PARAMETERS (your original values)
  # ============================================================
  
  parameters <- list(
    CLI = 162,    # Clearance from Central compartment (L/h)
    V2I = 1230,   # Volume of distribution of Central compartment (L)
    V3I = 2590,   # Volume of distribution Peripheral 1 (L)
    V4I = 4340,   # Volume of distribution Peripheral 2 (L)
    Q3I = 76.6,   # Intercompartmental clearance (Central <-> P1)
    Q4I = 3.26,   # Intercompartmental clearance (Central <-> P2)
    kaI = 1.84,   # First-order absorption rate constant (1/h)
    
    t0 = 0.348,          # Zero-order absorption duration (h)
    dose_times = dose_times,
    D0_a = D0_a
  )
  
  # ============================================================
  # INITIAL CONDITIONS (Depot receives zero-order input)
  # ============================================================
  
  state <- c(
    A_depot_2 = 0,
    A_central = 0,
    A_peripheral1 = 0,
    A_peripheral2 = 0
  )
  
  # ============================================================
  # PK MODEL (CORRECT MULTI-DOSE ZERO-ORDER)
  # ============================================================
  
  pk_model <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
      
      # -----------------------------
      # Zero-order input rate R_zero
      # -----------------------------
      R_zero <- 0
      
      for (i in seq_along(dose_times)) {
        t_i <- dose_times[i]
        D_i <- D0_a[i]
        
        if (time >= t_i && time < (t_i + t0)) {
          R_zero <- R_zero + (D_i / t0)
        }
      }
      
      # -----------------------------
      # ODEs
      # -----------------------------
      dA_depot_2   <- R_zero - kaI * A_depot_2
      
      dA_central <- kaI * A_depot_2 -
        (CLI / V2I) * A_central -
        (Q3I / V2I) * A_central + (Q3I / V3I) * A_peripheral1 -
        (Q4I / V2I) * A_central + (Q4I / V4I) * A_peripheral2
      
      dA_peripheral1 <- (Q3I / V2I) * A_central - (Q3I / V3I) * A_peripheral1
      dA_peripheral2 <- (Q4I / V2I) * A_central - (Q4I / V4I) * A_peripheral2
      
      list(c(dA_depot_2, dA_central, dA_peripheral1, dA_peripheral2))
    })
  }
  
  # ============================================================
  # SOLVE THE SYSTEM
  # ============================================================
  
  times <- seq(0, total_time, by = times_between_doses)
  
  output_df_AKTi <- as.data.frame(ode(
    y = state,
    times = times,
    func = pk_model,
    parms = parameters,
    method = "bdf",
    maxsteps = 5000
  ))
  
  # ============================================================
  # CONCENTRATIONS
  # ============================================================
  
  output_df_AKTi$C_total <- (output_df_AKTi$A_central + output_df_AKTi$A_peripheral1 + output_df_AKTi$A_peripheral2) /
    (parameters$V2I + parameters$V3I + parameters$V4I)
  
  output_df_AKTi$C_central <- output_df_AKTi$A_central / parameters$V2I
  
  # Convert to µM
  MW <- 458
  output_df_AKTi <- output_df_AKTi %>% mutate(C_total = (C_total * 1e6) / MW)
  
  
  
  # ============================================================
  # RETURN
  # ============================================================
  
  return(output_df_AKTi)
}


################################################################################
################################################################################

################################################################################
# AKTi schedule:
# - no dosing during first 5 pretreatment days
# - then 3 weeks
#   days_on = c(days in week1, days in week2, days in week3)
#   dose_a  = c(dose_week1, dose_week2, dose_week3)
################################################################################

generate_dosing_schedule_D0a <- function(dose_a, days_on,
                                         pretreat_days = 5, treat_days = 21) {
  
  total_days <- pretreat_days + treat_days
  schedule <- rep(0, total_days)
  
  # ---- Week 1 ----
  if (days_on[1] > 0) {
    schedule[(pretreat_days + 1):(pretreat_days + days_on[1])] <- dose_a[1]
  }
  
  # ---- Week 2 ----
  if (days_on[2] > 0) {
    week2_start <- pretreat_days + 8
    schedule[week2_start:(week2_start + days_on[2] - 1)] <- dose_a[2]
  }
  
  # ---- Week 3 ----
  if (days_on[3] > 0) {
    week3_start <- pretreat_days + 15
    schedule[week3_start:(week3_start + days_on[3] - 1)] <- dose_a[3]
  }
  
  return(schedule)
}




################################################################################
# EZH2i schedule:
# dose_e = c(dose_pretreat, dose_week1, dose_week2, dose_week3)
################################################################################

generate_dosing_schedule_D0e <- function(dose_e,
                                         pretreat_days = 5, treat_days = 21) {
  
  total_days <- pretreat_days + treat_days
  schedule <- numeric(total_days)
  
  # Pretreatment: days 1–5
  schedule[1:pretreat_days] <- dose_e[1]
  
  # Week 1: days 6–12
  schedule[(pretreat_days + 1):(pretreat_days + 7)] <- dose_e[2]
  
  # Week 2: days 13–19
  schedule[(pretreat_days + 7 + 1):(pretreat_days + 14)] <- dose_e[3]
  
  # Week 3: days 20–26
  schedule[(pretreat_days + 14 + 1):(pretreat_days + 21)] <- dose_e[4]
  
  return(schedule)
}



################################################################################
# BID administration 
################################################################################

expand_daily_doses <- function(d) {
  rep(d, each = 2)   # turns 26 days → 52 BID doses
}





#############################################################################





# SAVE FINALs

# Define file paths
tested_treatments_file <- "/Users/sbruno/Partners HealthCare Dropbox/Simone Bruno/Projects/Cichowski, Karen/Amy data (project 1)/R codes/data-analysis/GA_drug2BEFORE_autocatalysis_no_A_rep_E_2doses_per_day_CostF_dyingcells_AND_toxicity_2025_12_tested_treatments_2_bis.rds"
best_treatment_file <- "/Users/sbruno/Partners HealthCare Dropbox/Simone Bruno/Projects/Cichowski, Karen/Amy data (project 1)/R codes/data-analysis/GA_drug2BEFORE_autocatalysis_no_A_rep_E_2doses_per_day_CostF_dyingcells_AND_toxicity_2025_12_best_treatment_2_bis.rds"

# # # Save the objects
# saveRDS(tested_treatments, file = tested_treatments_file)
# saveRDS(best_treatment, file = best_treatment_file)

### To reload them later
tested_treatments <- readRDS(tested_treatments_file)
best_treatment <- readRDS(best_treatment_file)
# 
# 
set.seed(1234) # Ensure reproducibility




# ###################### SUMMARY PLOTs



tested_treatments <- tested_treatments %>%
  mutate(
    EZHi_total = round(Dose_EZHi_pt + Dose_EZHi_w1 + Dose_EZHi_w2 + Dose_EZHi_w3, 0),
    AKTi_total = round(Dose_AKTi_w1 + Dose_AKTi_w2 + Dose_AKTi_w3, 0),
    AKTi_days_total = AKTi_Days_Week1 + AKTi_Days_Week2 + AKTi_Days_Week3
  )

############################################################
# 2) Remove NEGATIVE cost-function solutions
############################################################


tested_treatments <- tested_treatments %>%
  filter(Cost_Function <= 0)


# tested_treatments <- tested_treatments[1:4025, ]

############################################################
# 3) Prepare BEST treatment with numeric axes
############################################################

best_treatment <- best_treatment %>%
  mutate(
    EZHi_total = round(Dose_EZHi_pt + Dose_EZHi_w1 + Dose_EZHi_w2 + Dose_EZHi_w3, 0),
    AKTi_total = round(Dose_AKTi_w1 + Dose_AKTi_w2 + Dose_AKTi_w3, 0),
    AKTi_days_total = AKTi_Days_Week1 + AKTi_Days_Week2 + AKTi_Days_Week3
  )

############################################################
# 4) FIX: ensure best_treatment is inside plotting range
############################################################

# Extend axis limits if needed
x_min <- min(tested_treatments$EZHi_total, best_treatment$EZHi_total)
x_max <- max(tested_treatments$EZHi_total, best_treatment$EZHi_total)

y_min <- min(tested_treatments$AKTi_total, best_treatment$AKTi_total)
y_max <- max(tested_treatments$AKTi_total, best_treatment$AKTi_total)


################################################################################
################################################################################
################################################################################
# RE-EVALUATE SELECTED SCHEDULES FOR WATERFALL PLOT
################################################################################

cat("Re-running selected treatment schedules on all patients...\n")

##################################################################
# SELECT N SCHEDULES
##################################################################


set.seed(1234) # Ensure reproducibility

N <- 5
n_patients_waterfall <- 1



selected_schedule_indices <- c(4034,  3646, 2827, 91)

selected_schedules <- tested_treatments[selected_schedule_indices, ]




##################################################################
# SIMULATION FUNCTION
##################################################################


simulate_schedule_for_patient <- function(D0_e_single, D0_a_single, days_on_vec, patient_params) {
  
  D0_e <- expand_daily_doses(generate_dosing_schedule_D0e(D0_e_single))
  D0_a <- expand_daily_doses(generate_dosing_schedule_D0a(D0_a_single, days_on_vec))
  
  dose_times <- sort(c(
    1, 1+12,
    1*24,1*24+12,
    2*24,2*24+12,
    3*24,3*24+12,
    4*24,4*24+12,
    5*24,5*24+12,
    6*24,6*24+12,
    7*24,7*24+12,
    8*24,8*24+12,
    9*24,9*24+12,
    10*24,10*24+12,
    11*24,11*24+12,
    12*24,12*24+12,
    13*24,13*24+12,
    14*24,14*24+12,
    15*24,15*24+12,
    16*24,16*24+12,
    17*24,17*24+12,
    18*24,18*24+12,
    19*24,19*24+12,
    20*24,20*24+12,
    21*24,21*24+12,
    22*24,22*24+12,
    23*24,23*24+12,
    24*24,24*24+12,
    25*24,25*24+12
  ))
  
  total_time <- 26 * 24
  
  # pkA <- define_PK_AKTi_fun(D0_a, dose_times[-1], total_time, 0.1)
  pkA <- define_PK_AKTi_fun(D0_a, dose_times, total_time, 0.1)
  pkE <- define_PK_EZHi_fun(D0_e, dose_times[-1], total_time, 0.1)
  
  pkA <- pkA[!duplicated(pkA$time), ]
  pkE <- pkE[!duplicated(pkE$time), ]
  
  AKTi_interp <- approxfun(pkA$time, pkA$C_total, rule = 2)
  EZHi_interp <- approxfun(pkE$time, pkE$C_plasma, rule = 2)
  
  times <- seq(0, 624, by = 0.1)
  y_initial <- c(y1=10.6765,y2=0,y3=0,y4=1,y5=1,y6=0,y7=0,y8=0,y9=0,y10=0)
  
  PD_model <- function(time, y, parms) {
    with(as.list(c(y, parms)), {
      k35 <- AKTi_interp(time)
      k36 <- EZHi_interp(time)
      dy1_dt <- (k1 * (y6 / (k2 + y6)) * exp(-k16 * time) * y8) - k15 * y1 * (1 - (y1 / 100))+ k32 * exp(-k16 * time) * y9  + k34 * exp(-k16 * time) * y10
      dy2_dt <- k3 * (k4 / (k4 + y4)) * y3 - k5 * y2
      dy3_dt <- (k6 +k25*y3) * (50 - y3) - (k29+k7 * (y5 / (k8 + y5)) *(1+k26*(50-y3)))* y3
      dy4_dt <- k9 - k10 * (k35 / (k11 + k35)) * y4 - k5 * y4
      dy5_dt <- k12 - k13 * (k36 / (k14 + k36)) * y5 - k5 * y5
      dy6_dt <- k19 * (k20 / (k20 + y4)) * y7 - k21 * y6
      dy7_dt <- (k22 + k27 * y7) * (50 - y7) - (k30 + k23 * (y5 / (k24 + y5)) * (1 + k28 * (50 - y7))) * y7
      dy8_dt <- k17 * (y2 / (k18 + y2)) * (100 - y1 - y8 - y9 - y10) - (k33 + k1 * (y6 / (k2 + y6)) * exp(-k16 * time) - k15 * (y1 / 100)) * y8
      dy9_dt  <- k31 * (100 - y1 - y8 - y9 - y10)- (k32 * exp(-k16 * time) - k15 * (y1 / 100)) * y9
      dy10_dt <- k33 * y8 - (k34 * exp(-k16 * time) - k15 * (y1 / 100)) * y10
      list(c(dy1_dt, dy2_dt, dy3_dt, dy4_dt, dy5_dt, dy6_dt, dy7_dt, dy8_dt, dy9_dt, dy10_dt))
    })
  }
  
  pd_output <- ode(y = y_initial, times = times, func = PD_model, parms = patient_params,method = "lsoda")
  
  return(tail(pd_output[, "y1"], 1))
}


##################################################################
# RUN ALL SELECTED SCHEDULES (PARALLEL)
##################################################################


# Always use the single patient
patient_indices_waterfall <- 1  

waterfall_data <- data.frame()
NN <- length(selected_schedule_indices)

for (s in seq_len(NN)) {
  
  row <- selected_schedules[s, ]
  
  D0_e_single <- c(row$Dose_EZHi_pt, row$Dose_EZHi_w1, row$Dose_EZHi_w2, row$Dose_EZHi_w3)
  D0_a_single <- c(row$Dose_AKTi_w1, row$Dose_AKTi_w2, row$Dose_AKTi_w3)
  days_on_vec <- c(row$AKTi_Days_Week1, row$AKTi_Days_Week2, row$AKTi_Days_Week3)
  
  
  
  results_list <- future_lapply(patient_indices_waterfall, function(i) {
    
    patient_par <- as.list(params_samples[i, ])
    
    y1_val <- simulate_schedule_for_patient(
      D0_e_single, D0_a_single, days_on_vec, patient_par
    )
    
    cost_val <- -(y1_val / 100)
    
    data.frame(
      schedule_rank = s,
      schedule_id = selected_schedule_indices[s],
      patient_id = i,
      y1_T624 = y1_val,
      cost = cost_val
    )
  })
  
  waterfall_data <- rbind(waterfall_data, do.call(rbind, results_list))
}


##################################################################
# BUILD waterfall_cost
##################################################################

waterfall_cost <- waterfall_data %>%
  arrange(cost) %>%
  mutate(sim_index = row_number())


##################################################################
# BUILD schedule_info FOR LEGEND
##################################################################

schedule_info <- tested_treatments %>%
  mutate(schedule_id = row_number()) %>%
  filter(schedule_id %in% unique(waterfall_cost$schedule_id)) %>%
  mutate(
    EZHi_total = Dose_EZHi_pt + Dose_EZHi_w1 + Dose_EZHi_w2 + Dose_EZHi_w3,
    AKTi_total = Dose_AKTi_w1 + Dose_AKTi_w2 + Dose_AKTi_w3
  ) %>%
  group_by(schedule_id) %>%
  summarise(
    avg_cost = mean(Cost_Function),
    EZHi_total = mean(EZHi_total),
    AKTi_total = mean(AKTi_total),
    .groups = "drop"
  ) %>%
  arrange(avg_cost) %>%                     # order schedules by cost (best → worst)
  mutate(
    legend_label = paste0(
      "ID ", schedule_id,
      " | cost=", round(avg_cost, 3),
      " | EZHi=", round(EZHi_total),
      " | AKTi=", round(AKTi_total)
    ),
    legend_label = factor(legend_label,     # enforce legend order by avg_cost
                          levels = legend_label)
  )

waterfall_cost <- waterfall_cost %>%
  left_join(schedule_info, by = "schedule_id")


##################################################################
#  FINAL WATERFALL PLOT
##################################################################

ggplot(waterfall_cost,
       aes(x = sim_index, y = -cost, fill = legend_label)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Waterfall Plot — Cost Function",
    subtitle = paste("Top", N, "Schedules Re-evaluated on",
                     n_patients_waterfall, "Patients"),
    x = "Simulated patients (ordered)",
    y = "Improvement (−cost = more dying cells)",
    fill = "Schedule (AvgCost | EZHi | AKTi)"
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.05, end = 0.9)+
  coord_cartesian(ylim = c(0, 0.8))












##################################################################
# FIXED SCHEDULE ORDER (DEFINE ONCE)
##################################################################

# Define reference order
schedule_order <- c(
  3646, 4034, 2827, 91
)


##################################################################
# APPLY FIXED ORDER TO WATERFALL DATA
##################################################################

waterfall_cost_fixed <- waterfall_cost %>%
  mutate(
    schedule_id_fixed = factor(schedule_id, levels = schedule_order)
  ) %>%
  arrange(schedule_id_fixed) %>%
  mutate(x_index = row_number())

##################################################################
# FIXED-ORDER WATERFALL PLOT
##################################################################

 ggplot(
  waterfall_cost_fixed,
  aes(x = x_index, y = -cost, fill = legend_label)
) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Waterfall Plot — Cost Function",
    subtitle = paste("Schedules shown in fixed reference order"),
    x = "Schedules (fixed order)",
    y = "Improvement (−cost = more dying cells)",
    fill = "Schedule (AvgCost | EZHi | AKTi)"
  ) +
  scale_fill_viridis_d(option = "plasma", begin = 0.05, end = 0.9) +
  coord_cartesian(ylim = c(0, 0.8))
