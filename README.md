# Reconstruction of RIR from spectrum having missing entries
This is an example implementation of the ADMM algorithm for reconstructing RIR (room impulse response) proposed in the following paper.

Kohei Yatabe and Akiko Sugahara, "Simulation of room impulse response using frequency-domain FEM and convex-optimization-based post-processing," submitted.

## Files
 - demo.m
	 - script for running the ADMM algorithm
 - RIRreconstADMM.m
	 - function of the ADMM algorithm
 - RIR_FDTD.mat
	 - data used in demo.m
