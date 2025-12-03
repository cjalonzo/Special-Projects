# Transmission Line Simulator (MATLAB)

This project implements a time- and frequency-domain transmission line simulator using MATLAB.  
The simulator models distributed transmission line parameters (R, L, G, C) using an N-section lumped LC/RG ladder and computes wave propagation, 
attenuation, reflections, and impedance characteristics.

---

## Features

###  **Distributed Parameter Modeling**
- Models transmission lines using per-unit-length **R, L, G, C** parameters.
- Supports skin-effect-dependent resistance (optional frequency-dependent model).
- Includes dielectric loss modeling via the **G** parameter.

###  **Time-Domain Simulation**
- Simulates pulse propagation along a discretized line.
- Captures reflections for open, short, and matched terminations.
- Produces waveform outputs analogous to **TDR/TDT measurements**.

###  **Frequency-Domain Analysis**
- Computes:
  - **Characteristic impedance (Z₀)**
  - **Propagation constant (γ = α + jβ)**
  - **Frequency-dependent attenuation due to R and G**
  - Optional skin effect correction:  
    \( R(f) = R_0 \sqrt{f/f_0} \)

###  **Visualization**
- Plots:
  - Voltage vs. time at various points along the line
  - Impedance and propagation constant vs. frequency
  - Simulated time-domain reflection behavior of digital transmission lines, producing waveforms conceptually similar to TDR responses used in impedance-discontinuity analysis.

---

## Skin Effect Modeling

Skin effect is incorporated by modifying the series resistance per unit length as a function of frequency:

```matlab
R_f = R0 * sqrt(f / f_ref);
