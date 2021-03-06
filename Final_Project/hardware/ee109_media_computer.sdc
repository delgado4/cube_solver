## Generated SDC file "ee109_media_computer.sdc"

## Copyright (C) 1991-2012 Altera Corporation
## Your use of Altera Corporation's design tools, logic functions 
## and other software and tools, and its AMPP partner logic 
## functions, and any output files from any of the foregoing 
## (including device programming or simulation files), and any 
## associated documentation or information are expressly subject 
## to the terms and conditions of the Altera Program License 
## Subscription Agreement, Altera MegaCore Function License 
## Agreement, or other applicable license agreement, including, 
## without limitation, that your use is for the sole purpose of 
## programming logic devices manufactured by Altera and sold by 
## Altera or its authorized distributors.  Please refer to the 
## applicable agreement for further details.


## VENDOR  "Altera"
## PROGRAM "Quartus II"
## VERSION "Version 12.1 Build 175 10/24/2012 SJ Full Version"

## DATE    "Fri May 24 13:03:52 2013"

##
## DEVICE  "EP4CE115F29C7"
##


#**************************************************************
# Time Information
#**************************************************************

set_time_format -unit ns -decimal_places 3



#**************************************************************
# Create Clock
#**************************************************************

create_clock -name {CLOCK_50} -period 20.000 -waveform { 0.000 10.000 } [get_ports {CLOCK_50}]
create_clock -name {TD_CLK27} -period 37.037 -waveform { 0.000 18.518 } [get_ports {TD_CLK27}]


#**************************************************************
# Create Generated Clock
#**************************************************************

#create_generated_clock -name {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|clk[0]} -source [get_pins {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|inclk[0]}] -duty_cycle 50.000 -multiply_by 1 -master_clock {CLOCK_50} [get_pins {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|clk[0]}] 
#create_generated_clock -name {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|clk[1]} -source [get_pins {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|inclk[0]}] -duty_cycle 50.000 -multiply_by 1 -phase -54.000 -master_clock {CLOCK_50} [get_pins {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|clk[1]}] 
#create_generated_clock -name {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|clk[2]} -source [get_pins {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|inclk[0]}] -duty_cycle 50.000 -multiply_by 1 -divide_by 2 -phase 180.000 -master_clock {CLOCK_50} [get_pins {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|clk[2]}] 
#create_generated_clock -name {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|clk[3]} -source [get_pins {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|inclk[0]}] -duty_cycle 50.000 -multiply_by 2 -divide_by 3 -master_clock {CLOCK_50} [get_pins {NiosII|external_clocks|DE_Clock_Generator_System|auto_generated|pll1|clk[3]}] 
#create_generated_clock -name {NiosII|external_clocks|DE_Clock_Generator_Audio|auto_generated|pll1|clk[1]} -source [get_pins {NiosII|external_clocks|DE_Clock_Generator_Audio|auto_generated|pll1|inclk[0]}] -duty_cycle 50.000 -multiply_by 14 -divide_by 31 -master_clock {TD_CLK27} [get_pins {NiosII|external_clocks|DE_Clock_Generator_Audio|auto_generated|pll1|clk[1]}] 

derive_pll_clocks

set_clock_groups -exclusive -group [get_clocks pll_inst|altpll_component|auto_generated|pll1|clk[0]] -group [get_clocks pll_inst|altpll_component|auto_generated|pll1|clk[1]] -group [get_clocks pll_inst|altpll_component|auto_generated|pll1|clk[2]] -group [get_clocks pll_inst|altpll_component|auto_generated|pll1|clk[3]]



#**************************************************************
# Set Clock Latency
#**************************************************************



#**************************************************************
# Set Clock Uncertainty
#**************************************************************
derive_clock_uncertainty



#**************************************************************
# Set Input Delay
#**************************************************************



#**************************************************************
# Set Output Delay
#**************************************************************



#**************************************************************
# Set Clock Groups
#**************************************************************



#**************************************************************
# Set False Path
#**************************************************************



#**************************************************************
# Set Multicycle Path
#**************************************************************



#**************************************************************
# Set Maximum Delay
#**************************************************************



#**************************************************************
# Set Minimum Delay
#**************************************************************



#**************************************************************
# Set Input Transition
#**************************************************************

