# Simplified SDRE

SDRE is basically a sub-optimal algorithm. To the best of our knowledge, we are the first to eliminate the sub-optimality of the algorithm and proposed "Simplified SDRE" for both regulation and tracking problems. Simplified SDRE was verified by implementing it in wind energy conversion systems and robotic hands. Please go through the following papers to know the details about existing SDRE and our proposed Simplified SDRE. 

## Simplified SDRE for Regulation
https://ieeexplore.ieee.org/document/8834201

## Simplified SDRE for Tracking
https://ieeexplore.ieee.org/document/9123230

# Running the Code

Simplified SDRE for Regulation (written in MATLAB):
1) First, run the file- Simp_SDRE_Regulation_P1_lqrnss.m
2) Then, run the file- Simplified_SDRE_Regulation_Part2.m 

Simplified SDRE for Tracking (written in MATLAB):
1) First, run the file- Simp_SDRE_Tracking_P1_lqtnss.m
2) Then, run other file- Simplified_SDRE_Trcking_Part2.m

Existing SDRE (written in Python):
* SDRE_Optimal_Control_Simulator.ipynb file works as a SDRE simulator. You just need to input the system you are working with, your desired output, and the controller parameters. The simulator will give you the control signal and the tracking error. 


## Programming Platform:
* MATLAB/Simulink
* Python




