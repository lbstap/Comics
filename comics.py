###########################################################################################################################
###                                COnceptual Model of Ice-sheet ChangeablenesS (COMICS)                                ###
###########################################################################################################################
# Execute: $ python comics.py
#
# Program to calculate transient ice volume of an ice sheet, depending on a prescribed
# control parameter - equilibrium ice volume (C-Veq) relation and evolution of the control parameter
# Basic procedure:
# If the ice volume is smaller than the equilibrium volume, it grows by the growth rate.
# If the ice volume is larger than the equilibrium volume, it shrinks by the decay rate.
#
# Cite:
#
#  Stap, L.B., Knorr, G., and Lohmann, G.: Anti-phased Miocene ice volume and CO2 changes by transient Antarctic ice sheet
#    variability, Paleoceanography and Paleoclimatology, in review.
#
# -------------------------------------------------------------------------------------------------------------------------
# The program consists of four phases:
# 1) Set-up
# 2) Initialization
# 3) Calculation of time evolution (model core)
# 4) Plot the main results
# -------------------------------------------------------------------------------------------------------------------------
# Written by: L.B. Stap at Alfred-Wegener-Institut, Helmholtz-Zentrum fuer Polar- und Meeresforschung, Bremerhaven, Germany
# Email: lennert.stap<at>gmail.com
# July 2019, last change: October 2020


import numpy as np
import matplotlib.pyplot as plt
import random
import math
import sys


#################################################
### 1) SET-UP                                 ###
#################################################
# Select the settings for your simulation in (i)
# -----------------------------------------------
# i)   Settings
# ii)  C-Veq relations
# iii) Control parameter evolution
# iv)  Growth and decay rates


# I) SETTINGS:

eq_profile                   = 'simple'                  # Control parameter - equilibrium ice volume relation, options: simple, hysteresis, PISM, custom
set_growth_rate              = 0.002                     # Sets the growth rate, not used if ctl_profile = PISM
set_decay_rate               = 0.004                     # Sets the decay rate, not used if ctl_profile = PISM
ctl_profile                  = 'triangle'                # Evolution of the control parameter, options: triangle, custom
number_of_timesteps          = 400                       # Length of the control parameter evolution. 40 kyr = 400, 100 kyr = 1000, 400 kyr = 4000
time_max                     = 5*number_of_timesteps     # Model integration length
start_option                 = 'zero'                    # Initial ice volume, default option is custom start, other options: 'warm', 'cold'
timestep_length              = 1                         # Length of a timestep (for plotting purposes only, does not affect the ice volume results), for PISM based runs: 100 yrs
just_show_last_cycle         = 'no'                     # For plot 3 (ice volume against control parameter), opt for 'yes' if you just want to plot the final forcing cycle

# Comparing two runs:
choose_second_run            = 'no'                      # If you want to compare two evolutions of ice volume, opt for 'yes'
second_run_option            = 'same'                    # Control parameter evolution for the second run, default option is the same, other options: 'amplitude_reduced' for reduced amplitude, 'period_reduced' for reduced period (increased frequency)
start_option_2nd             = 'zero'                    # Initial ice volume for the second run, default option is custom start, other options: 'warm', 'cold'
period_reduction_factor      = 2                         # Reduction factor of the period, in case of 'period_reduced'
amplitude_reduction_factor   = 0.5                       # Reduction factor of the amplitude, in case of 'amplitude_reduced'
amplitude_center             = 0.5                       # Mean of the control parameter variability, in case of 'amplitude_reduced'
multiply_decay               = 1.0                       # Multiply the specified decay rate by this factor for the second run
multiply_growth              = 1.0                       # Multiply the specified growth rate by this factor for the second run


# II) C-Veq RELATIONS:

# Simple symmetric linear C-Veq relation without hysteresis
if eq_profile == 'simple':
   def equilibrium_icevolume_top(control): # top of hysteresis curve
     return 1.-control

   def equilibrium_icevolume_bottom(control): # bottom of the hysteresis curve
     return 1.-control

   largest_volume=1. # For plotting purposes

# Simple piecewise symmetric linear C-Veq relation with hysteresis
elif eq_profile == 'hysteresis':
   def equilibrium_icevolume_top(control): # top of hysteresis curve
     if control > 0.5 :
       return (1.-control)*1.4
     else :
	   return 1.-(control*0.6)

   def equilibrium_icevolume_bottom(control): # bottom of the hysteresis curve
     if control > 0.5 :
       return (1.-control)*0.6
     else :
	   return 1.-(control*1.4)

   largest_volume=1. # For plotting purposes

# PISM C-Veq relation, derived from the PISM results of Stap et al. (2019), GRL
elif eq_profile == 'PISM':
  def equilibrium_icevolume_top(control): # top of hysteresis curve
    if control < 0.2671 :
      eq_vol_top = 18.1 - ((18.1 - 16.7)*(control)/0.2671)
    elif control >= 0.2671 and control < 0.5342 :
      eq_vol_top = 16.7 - ((16.7 - 6.5)*(control - 0.2671)/0.2671)
    elif control >= 0.5342 and control < 0.7671 :
      eq_vol_top = 6.5 - ((6.5 - 3.3)*(control - 0.5342)/0.2329)
    else:
      eq_vol_top = 3.3 - ((3.3 - 1.8)*(control - 0.7671)/0.2329)
    return eq_vol_top

  def equilibrium_icevolume_bottom(control): # top of hysteresis curve
    if control < 0.2671 :
      eq_vol_bot = 18.1 - ((18.1 - 15.1)*(control)/0.2671)
    elif control >= 0.2671 and control < 0.5342 :
      eq_vol_bot = 15.1 - ((15.1 - 5.1)*(control - 0.2671)/0.2671)
    elif control >= 0.5342 and control < 0.7671 :
      eq_vol_bot = 5.1 - ((5.1 - 2.2)*(control - 0.5342)/0.2329)
    else:
      eq_vol_bot = 2.2 - ((2.2 - 1.3)*(control - 0.7671)/0.2329)
    return eq_vol_bot

  largest_volume=20. # For plotting purposes

# Customize your C-Veq relation
elif eq_profile == 'custom':
   def equilibrium_icevolume_top(control): # top of hysteresis curve
     return 1.-control

   def equilibrium_icevolume_bottom(control): # bottom of the hysteresis curve
     return 1.-control

else:
  sys.exit("ERROR: Undefined C-Veq relation")

# III) CONTROL PARAMETER EVOLUTION:

# Simple 1 -> 0 -> 1 (normal) or 0 -> 1 -> 0 (reverse) triangle control parameter evolution, over the course of $number_of_timesteps$ time steps
if ctl_profile == 'triangle':
  def control(t):
    if (t % number_of_timesteps) < (number_of_timesteps/2):
      control=1.-((t % number_of_timesteps)*1./(number_of_timesteps/2))
    else:
      control=(((t % number_of_timesteps)- (number_of_timesteps/2))*1./(number_of_timesteps/2))
    if start_option == 'cold':
      return (control*-1.)+1. # reverse (cold start)
    else:
      return control # normal (warm, or zero start)

# Customize your control parameter evolution, e.g. by using sine waves
elif ctl_profile == 'custom':
  def control(t):
    precession       = ((math.sin(((t-0)*2.*math.pi)/230 )+1.)/2.)*0.00
    obliquity        = ((math.sin(((t-0)*2.*math.pi)/410 )+1.)/2.)*0.50
    eccentricity     = ((math.sin(((t-0)*2.*math.pi)/1050)+1.)/2.)*0.00
    eccentricity_mod = ((math.sin(((t-0)*2.*math.pi)/4000)+1.)/2.)*0.00
    control = precession + obliquity + eccentricity + eccentricity_mod
    return control

else:
  sys.exit("ERROR: Undefined control parameter evolution")

# IV) GROWTH AND DECAY RATES:

def growth_rate(control,icevolume,equilibrium_icesheet_bottom):
  if eq_profile == 'PISM':
    return 0.004+0.012*math.sqrt((max(0,equilibrium_icesheet_bottom-icevolume-9))) # tuned PISM profile
  else:
    return set_growth_rate

def decay_rate(control,icevolume,equilibrium_icesheet_top):
  if eq_profile == 'PISM':
    return 0.01+(icevolume-equilibrium_icesheet_top)*0.01 # tuned PISM profile
  else:
    return set_decay_rate


#############################################################################
### 2) INITIALIZATION                                                     ###
#############################################################################
# Initializing the time parameter, the control parameter,
# the hysteresis curves (top and bottom branch), and the transient ice volume

time = []
time = [0 for i in range(time_max)]
time[0] = 1

control_profile = []
control_profile = [0 for i in range(time_max)]
if start_option == 'cold':
  control_profile[0] = 0.0
else:
  control_profile[0] = 1.0

equilibrium_profile_top = []
equilibrium_profile_top = [0 for i in range(time_max)]
equilibrium_profile_top[0] = 0.0

equilibrium_profile_bottom = []
equilibrium_profile_bottom = [0 for i in range(time_max)]
equilibrium_profile_bottom[0] = 0.0

icevolume = []
icevolume = [0 for i in range(time_max)]
if start_option == 'warm':
  icevolume[0] = 1.3 # warm start (PISM derived)
elif start_option == 'cold':
  icevolume[0] = 23.0 # cold start (PISM derived)
else:
  icevolume[0] = 0.0 # default custom start, e.g. 0.0

if choose_second_run == 'yes':

  control_2nd = []
  control_2nd = [0 for i in range(time_max)]
  if start_option_2nd == 'cold':
    control_2nd[0] = 0.0
  else:
    control_2nd[0] = 1.0

  equilibrium_profile_top_2nd = []
  equilibrium_profile_top_2nd = [0 for i in range(time_max)]
  equilibrium_profile_top_2nd[0] = 0.0

  equilibrium_profile_bottom_2nd = []
  equilibrium_profile_bottom_2nd = [0 for i in range(time_max)]
  equilibrium_profile_bottom_2nd[0] = 0.0

  icevolume_2nd = []
  icevolume_2nd = [0 for i in range(time_max)]
  if start_option_2nd == 'warm':
    icevolume_2nd[0] = 1.3 # warm start (PISM derived)
  elif start_option_2nd == 'cold':
    icevolume_2nd[0] = 23.0 # cold start (PISM derived)
  else:
    icevolume_2nd[0] = 0.0 # default custom start, e.g. 0.0

#####################################################
### 3) CALCULATION OF TIME EVOLUTION (MODEL CORE) ###
#####################################################

for t in range(1,time_max):

  time[t]=t * timestep_length
  control_profile[t] = control(t)
  equilibrium_profile_top[t] = equilibrium_icevolume_top(control_profile[t])
  equilibrium_profile_bottom[t] = equilibrium_icevolume_bottom(control_profile[t])

  if icevolume[t-1] > equilibrium_icevolume_top(control_profile[t]):
    icevolume[t] = icevolume[t-1] - decay_rate(control_profile[t],icevolume[t-1],equilibrium_profile_top[t])
  elif icevolume[t-1] < equilibrium_icevolume_bottom(control_profile[t]):
    icevolume[t] = icevolume[t-1] + growth_rate(control_profile[t],icevolume[t-1],equilibrium_profile_bottom[t])
  else:
    icevolume[t] = icevolume[t-1]

  if choose_second_run == 'yes':

    if second_run_option == 'amplitude_reduced':
      control_2nd[t] = (amplitude_reduction_factor*control(t))+(amplitude_center-(amplitude_reduction_factor/2.))
    elif second_run_option == 'period_reduced':
      control_2nd[t] = control(period_reduction_factor*t)
    else:
      control_2nd[t] = control(t)

    if not start_option == start_option_2nd:
      control_2nd[t] = (control_2nd[t]*-1.)+1.

    equilibrium_profile_top_2nd[t] = equilibrium_icevolume_top(control_2nd[t])
    equilibrium_profile_bottom_2nd[t] = equilibrium_icevolume_bottom(control_2nd[t])

    if icevolume_2nd[t-1] > equilibrium_icevolume_top(control_2nd[t]):
      icevolume_2nd[t] = icevolume_2nd[t-1] - multiply_decay*decay_rate(control_2nd[t],icevolume_2nd[t-1],equilibrium_profile_top_2nd[t])
    elif icevolume_2nd[t-1] < equilibrium_icevolume_bottom(control_2nd[t]):
      icevolume_2nd[t] = icevolume_2nd[t-1] + multiply_growth*growth_rate(control_2nd[t],icevolume_2nd[t-1],equilibrium_profile_bottom_2nd[t])
    else:
      icevolume_2nd[t] = icevolume_2nd[t-1]


#######################################
### 4) PLOT THE MAIN RESULTS        ###
#######################################
# i)   Control parameter vs. time
# ii)  Ice volume vs. time
# iii) Ice volume vs. control parameter

# i) Plot of the control parameter against time
fig1 = plt.figure(1,figsize=(12,7))
plt.plot(time[1:],control_profile[1:], color="black", linewidth=2.0)
if choose_second_run == 'yes':
  plt.plot(time[1:],control_2nd[1:], color="black", linewidth=2.0, dashes=[5,5])
plt.xlabel('Time', fontsize=30)
plt.ylabel('Control parameter', fontsize=30)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.ylim(0, 1)
fig1.subplots_adjust(bottom=0.20)

# ii) Plot of ice volume against time
fig2 = plt.figure(2,figsize=(12,7))
plt.plot(time[1:],icevolume[1:], color="black", linewidth=2.0)
plt.plot(time[1:],equilibrium_profile_top[1:], color="red", linewidth=2.0)
plt.plot(time[1:],equilibrium_profile_bottom[1:], color="blue", linewidth=2.0)
if choose_second_run == 'yes':
  plt.plot(time[1:],icevolume_2nd[1:], color="black", linewidth=2.0, dashes=[5,5])
  plt.plot(time[1:],equilibrium_profile_top_2nd[1:], color="red", linewidth=2.0, dashes=[5,5])
  plt.plot(time[1:],equilibrium_profile_bottom_2nd[1:], color="blue", linewidth=2.0, dashes=[5,5])
plt.xlabel('Time', fontsize=30)
plt.ylabel('Ice volume', fontsize=30)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
fig2.subplots_adjust(bottom=0.20)

# iii) Plot of the ice volume against the control parameter
fig3 = plt.figure(3,figsize=(12,7))
cyclestart = 1
cyclestop  = (number_of_timesteps)/2
plt.plot(control_profile[cyclestart:cyclestop],equilibrium_profile_bottom[cyclestart:cyclestop], color="blue", linewidth=2.0)
plt.plot(control_profile[cyclestart:cyclestop],equilibrium_profile_top[cyclestart:cyclestop], color="red", linewidth=2.0)
plt.fill_between(control_profile[cyclestart:cyclestop], 0., equilibrium_profile_bottom[cyclestart:cyclestop], color='blue', alpha='0.2')
plt.fill_between(control_profile[cyclestart:cyclestop], equilibrium_profile_top[cyclestart:cyclestop], largest_volume, color='red', alpha='0.2')
if just_show_last_cycle == 'yes' and time_max > number_of_timesteps:
  plt.plot(control_profile[(time_max-number_of_timesteps)-1:],icevolume[(time_max-number_of_timesteps)-1:], color="black", linewidth=2.0)
else:
  plt.plot(control_profile[1:],icevolume[1:], color="black", linewidth=2.0)
if choose_second_run == 'yes':
  if just_show_last_cycle == 'yes' and time_max > number_of_timesteps:
    plt.plot(control_2nd[(time_max-number_of_timesteps)-1:],icevolume_2nd[(time_max-number_of_timesteps)-1:], color="black", linewidth=2.0, dashes=[5,5])
  else:
    plt.plot(control_2nd[1:],icevolume_2nd[1:], color="black", linewidth=2.0, dashes=[5,5])
plt.xlabel('Control parameter', fontsize=30)
plt.ylabel('Ice volume', fontsize=30)
plt.xticks(fontsize = 25)
plt.yticks(fontsize = 25)
plt.ylim(0, largest_volume)
fig3.subplots_adjust(bottom=0.20)

# Save and show the plots
fig1.savefig('plot1_time_C.pdf')
fig2.savefig('plot2_time_V.pdf')
fig3.savefig('plot3_C_V.pdf')
plt.show()
