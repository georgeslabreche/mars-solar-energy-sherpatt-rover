library(mars)
Sys.setenv(NET_FLUX_FUNCTION_TYPE = "polynomial")

# Spacecraft properties.
spacecrafts = list(
  "VL1" = list(
    "latitude" = 22.3,
    "longitude" = -49.97,
    "beta_optimal" = 6.5
  ),
  "VL2" = list(
    "latitude" = 47.7,
    "longitude" = 134.29,
    "beta_optimal" = 22
  )
)

# Expected data
base_url = "https://raw.githubusercontent.com/georgeslabreche/mars/master/tests/testthat/data/"
expected_data = list(
  "VL1" = list(
    "tau" = read.csv(paste(base_url, "discretized_tau_at_vl1_fig_3_1991_update.csv", sep="")),
    "insolation" =
      list(
        "horizontal" = read.csv(paste(base_url, "daily_insolation_on_horizontal_surface_at_vl1_table_v_update_1990.csv", sep="")),
        "beta_equals_phi" = read.csv(paste(base_url, "daily_insolation_on_an_inclined_surface_for_beta_equals_phi_at_vl1_table_ii_1993.csv", sep="")),
        "beta_optimal" = read.csv(paste(base_url, "daily_insolation_on_optimal_inclined_angle_beta_at_vl1_table_iii_1993.csv", sep=""))
      )
  ),
  "VL2" = list(
    "tau" = read.csv(paste(base_url, "discretized_tau_at_vl2_fig_3_1991_update.csv", sep="")),
    "insolation" =
      list(
        "horizontal" = read.csv(paste(base_url, "daily_insolation_on_horizontal_surface_at_vl2_table_v_update_1990.csv", sep="")),
        "beta_equals_phi" = read.csv(paste(base_url, "daily_insolation_on_an_inclined_surface_for_beta_equals_phi_at_vl2_table_ii_1993.csv", sep="")),
        "beta_optimal" = read.csv(paste(base_url, "daily_insolation_on_optimal_inclined_angle_beta_at_vl2_table_iii_1993.csv", sep=""))
      )
  )
)

#' Calculate and plot differences between expected and calculated daily insolations surface.
#'
#' @param title 
#' @param insolation_type 
#' @param spacecraft 
#' @param Ls_seq 
#' @param beta_equals_phi 
#'
#' @return
plot_insolation_differences = function(title=title, insolation_type, spacecraft, Ls_seq, beta_equals_phi=FALSE){
  
  # Location.
  phi = spacecrafts[[spacecraft]][["latitude"]]
  longitude = spacecrafts[[spacecraft]][["longitude"]]
  
  # Orientation.
  gamma_c = 0
  
  # Inclination.
  beta = NULL

  # Optical depth as measured by the spacecraft (discretized).
  measured_taus = expected_data[[spacecraft]][["tau"]]
  
  # Get the expected insolations as measured by the spacecraft.
  # Also set the inclination angle when appropriate.
  expected_insolations = NULL
  if(insolation_type %in% c("Hbh", "Hdh", "Hh")){
    expected_insolations = expected_data[[spacecraft]][["insolation"]][["horizontal"]]
    
  }else if(insolation_type %in% c("Hb", "Hd", "Hal", "H")){
    if(isTRUE(beta_equals_phi)){
      expected_insolations = expected_data[[spacecraft]][["insolation"]][["beta_equals_phi"]]
      beta = phi
      
    }else{
      expected_insolations = expected_data[[spacecraft]][["insolation"]][["beta_optimal"]]
      beta = spacecrafts[[spacecraft]][["beta_optimal"]]
    }
  }else{
    stop("Invalid insolation type.")
  }
  
  calculated_insolations = c()
  taus = c()
  diffs = c()
  diffs_pc = c()
  
  for(Ls in Ls_seq){
    # Measured tau.
    tau_measured = measured_taus[measured_taus$Ls == Ls, "tau"]
    
    # Expected insolation.
    insolation_expected = expected_insolations[expected_insolations$Ls == Ls, insolation_type]
    
    # Calculated insolation.
    insolation_calculated = NULL
    if(insolation_type == "Hbh"){ # Beam daily insolation on a horizontal surface.
      insolation_calculated = H_bh(Ls=Ls, phi=phi, tau=tau_measured)
      
    }else if(insolation_type == "Hdh"){ # Diffuse daily insolation on a horizontal surface.
      insolation_calculated = H_dh(Ls=Ls, phi=phi, longitude=longitude, tau=tau_measured)
      
    }else if(insolation_type == "Hh"){ # Global daily insolation on a horizontal surface.
      insolation_calculated = H_h(Ls=Ls, phi=phi, longitude=longitude, tau=tau_measured)
      
    }else if(insolation_type == "Hb"){ # Beam daily insolation on an inclined surface.
      insolation_calculated = H_bi(Ls=Ls, phi=phi, tau=tau_measured, beta=beta, gamma_c=gamma_c)
        
    }else if(insolation_type == "Hd"){ # Diffuse daily insolation on an inclined surface.
      insolation_calculated = H_di(Ls=Ls, phi=phi, longitude=longitude, tau=tau_measured, beta=beta, gamma_c=gamma_c)
        
    }else if(insolation_type == "Hal"){ # Albedo daily insolation on an inclined surface.
      insolation_calculated = H_ali(Ls=Ls, phi=phi, longitude=longitude, tau=tau_measured, beta=beta, gamma_c=gamma_c)
        
    }else if(insolation_type == "H"){ # Global daily insolation on an inclined surface.
      insolation_calculated = H_i(Ls=Ls, phi=phi, longitude=longitude, tau=tau_measured, beta=beta, gamma_c=gamma_c)
        
    }else{
      stop("Invalid insolation type.")
    }
    
    # Collect data to be plotted later.
    calculated_insolations = c(calculated_insolations, insolation_calculated)
    taus = c(taus, tau_measured)
    
    diff = insolation_expected-insolation_calculated
    diffs = c(diffs, diff)
    
    diff_pc = ((insolation_expected-insolation_calculated) / insolation_calculated) * 100
    diffs_pc = c(diffs_pc, diff_pc)
    
  }
  
  # Plot the results.
  dev.new()
  par(mfrow=c(2,2))
  
  # Expected vs calculated insolations.
  plot(Ls_seq, expected_insolations[, insolation_type],
       xlab="Ls [deg]", ylab=paste(insolation_type, "[Wh/m2]"),
       col="red", type="l", main=title)
  lines(Ls_seq, calculated_insolations, col="blue")
  
  # Optical depths.
  plot(Ls_seq, taus,
       xlab="Ls [deg]", ylab="tau",
       col="orange", type="l", main=title)
  
  # Difference between expected and calculated insolations [Wh/m2].
  plot(taus, diffs,
       xlab="tau", ylab="difference [Wh/m2]",
       col="red", main=title)
  
  # Difference between expected and calculated insolations [%].
  plot(taus, diffs_pc,
       xlab="tau", ylab="difference [%]",
       col="red", main=title)
}

#' Calculate and plot differences between expected and calculated daily insolations on a horizontal surface.
#'
#' @param location 
#'
#' @return
plot_insolation_on_horizontal_surface = function(location){
  
  Ls_seq = seq(0, 355, 5)
  
  # Global insolation.
  plot_insolation_differences(
    title = paste("Hh at", location),
    insolation_type = "Hh",
    spacecraft = location,
    Ls_seq = Ls_seq)

  # Beam insolation.
  plot_insolation_differences(
    title = paste("Hbh at", location),
    insolation_type = "Hbh",
    spacecraft = location,
    Ls_seq = Ls_seq)

  # Diffuse insolation.
  plot_insolation_differences(
    title = paste("Hdh at", location),
    insolation_type = "Hdh",
    spacecraft = location,
    Ls_seq = Ls_seq)
}

#' Calculate and plot differences between expected and calculated daily insolations on a inclined surface.
#'
#' @param location 
#' @param beta_equals_phi 
#'
#' @return
plot_insolation_on_inclined_surface = function(location, beta_equals_phi=FALSE){
 
  Ls_seq = seq(0, 355, 5)
  if(isTRUE(beta_equals_phi)){
    Ls_seq = c(Ls_seq, 360)
  }
  
  # Beam insolation.
  plot_insolation_differences(
    title = paste("Hbi at", location),
    insolation_type = "Hb",
    spacecraft = location,
    Ls_seq = Ls_seq,
    beta_equals_phi = beta_equals_phi)
  
  # Diffuse insolation.
  plot_insolation_differences(
    title = paste("Hdi at", location),
    insolation_type = "Hd",
    spacecraft = location,
    Ls_seq = Ls_seq,
    beta_equals_phi = beta_equals_phi)
  
  # Albedo insolation.
  plot_insolation_differences(
    title = paste("Hali at", location),
    insolation_type = "Hal",
    spacecraft = location,
    Ls_seq = Ls_seq,
    beta_equals_phi = beta_equals_phi)
  
  # Global insolation.
  plot_insolation_differences(
    title = paste("Hi at", location),
    insolation_type = "H",
    spacecraft = location,
    Ls_seq = Ls_seq,
    beta_equals_phi = beta_equals_phi)
}

############################################
# Insolation on horizontal surface at VL1. #
############################################
#plot_insolation_on_horizontal_surface(locaiton="VL1")

############################################
# Insolation on horizontal surface at VL2. #
############################################
#plot_insolation_on_horizontal_surface(locaiton="VL2")

##########################################################
# Insolation on inclined surface with beta = phi at VL1. #
##########################################################
#plot_insolation_on_inclined_surface(location="VL1", beta_equals_phi=TRUE)

##########################################################
# Insolation on inclined surface with beta = phi at VL2. #
##########################################################
#plot_insolation_on_inclined_surface(location="VL2", beta_equals_phi=TRUE)

##########################################################
# Insolation on inclined surface with beta = phi at VL1. #
##########################################################
#plot_insolation_on_inclined_surface(location="VL1", beta_equals_phi=FALSE)

##########################################################
# Insolation on inclined surface with beta = phi at VL2. #
##########################################################
#plot_insolation_on_inclined_surface(location="VL2", beta_equals_phi=FALSE)
