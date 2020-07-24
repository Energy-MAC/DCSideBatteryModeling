# Grid-Coupled Dynamic Response of Battery-Driven Voltage Source Converters

This repository accompanies the paper *Grid-Coupled Dynamic Response of Battery-Driven Voltage Source Converters*
https://arxiv.org/abs/2007.11776

*Abstract* With the increasing interest in converter-fed islanded microgrids, particularly for resilience, it is becoming more critical to understand the dynamical behavior of these systems. This paper takes a holistic view of grid-forming converters and considers control approaches for both modeling and regulating the DC-link voltage when the DC-source is a battery energy storage system. We are specifically interested in understanding the performance of these controllers, subject to large load changes, for decreasing values of the DC-side capacitance. We consider a fourth, second, and zero-order model of the battery; and establish that the zero-order model captures the dynamics of interest for the timescales considered for disturbances examined. Additionally, we adapt a grid search for optimizing the controller parameters of the DC/DC controller and show how the inclusion of AC side measurements into the DC/DC controller can improve its dynamic performance. This improvement in performance offers the opportunity to reduce the DC-side capacitance given an admissible DC voltage transient deviation, thereby, potentially allowing for more reliable capacitor technology to be deployed.

- The file [parameter_search.ipynb](https://github.com/Energy-MAC/DCSideBatteryModeling/blob/master/parameter_search.ipynb) contains eigenvalue analysis used in the Small-signal tunning section of the results.

- The File [DAE_simulation.ipynb](https://github.com/Energy-MAC/DCSideBatteryModeling/blob/master/DAE_simulation.ipynb) contains the time domain simulation results relevant to the Large-signal tunning of the results.


Corresponding author: ciaranr_r@berkeley.edu
