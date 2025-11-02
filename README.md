ENUME 2025 — MATLAB Assignments (A/B/C)

Project overview
This repository contains three MATLAB assignments and their code:
A: Estimate left-ankle dorsiflexion angle from 3D motion-capture by building LT/LP frames; add Gaussian noise, denoise with a low-pass filter, and select the optimal cutoff via MAE.
B: Fourier-series approximation of an exponential pulse; compute coefficients via a rectangle rule and via integral(), compare coefficient/reconstruction errors versus step size and truncation order; observe the Gibbs phenomenon.
C: Extended Lotka–Volterra predator–prey model with a quadratic density term; implement ode45, Gear2, and RK4; perform parameter fitting (fminsearch) and compare accuracy–runtime tradeoffs.

Environment

MATLAB R2023a or newer.

Assignment A needs the provided dataset and a low-pass filter helper (e.g., typical_gait.mat and lpfilt.m or equivalent).

Assignment C needs animals_*.csv; other functions (ode45, odeset, deval, fminsearch) are built-in.

Quick start
A (Gait angle + filtering): place typical_gait.mat in data, run the main script to compute clean/noisy angles, sweep cutoff frequencies, compute MAE, and output the best cutoff and comparison plots.
B (Fourier approximation): run the main script to compute a_n, b_n (rectangle rule vs integral), reconstruct signals for multiple truncation orders M, and report coefficient/reconstruction error curves.
C (Extended LV, fitting, and method comparison): solve with ode45/Gear2/RK4; fit parameters to the given data; sweep tolerances or step sizes to generate accuracy–time curves and summarize which solver is best under which conditions.
