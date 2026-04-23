# Tracer Mass Conservation Report - SAS Convection Schemes

## Overview
This report summarizes the verification of the structural upgrade to the UFS CCPP-physics repository to fix tracer mass conservation leaks in the Simplified Arakawa-Schubert (SAS) convection schemes. The previous advective finite differencing and clipping method has been replaced with strictly conservative flux-form and hole-filling methods, supported by a global proportional mass fixer. Additionally, the transport routines have been optimized to handle multiple species simultaneously via Fortran interfaces (overloading).

## 1. Standalone Transport Module Verification
The `tracer_transport_mod` was tested using a 1D column with 10 layers. A consistent initial tracer distribution containing a negative value (a "hole") was used to compare how each method handles non-physical values and maintains total mass.

### Scenario: Fillable Hole (Positive Column Mass)
Initial state: $q = 10^{-3}$ everywhere except $q_5 = -10^{-4}$.

| Method | Initial Mass | Final Mass | Mass Change (Delta) | Min(q) | Status |
| :--- | ---: | ---: | ---: | ---: | :--- |
| **Flux-Form** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | +1.000E-04 | **Conservative** |
| **Hole-Filling Only** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | 0.000E+00 | **Conservative** |
| **Hole-Fill + Fixer** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | 0.000E+00 | **Conservative** |
| **Original Advective** | 9.0754742955 | 9.1774460188 | +1.01972E-01 | +1.000E-10 | **Leaky** (Clipped) |

*Note: The Hole-Filling routine alone is perfectly conservative because it redistributes mass to fill the deficit. The Advective method generates mass by clipping the negative value.*

## 2. Integrated SAS Scheme Verification (Loaded Gun Profile)
The refactored deep and shallow convection schemes were verified using the "Loaded Gun" severe convection profile to ensure correct integration and end-to-end conservation under realistic atmospheric conditions.

| Scheme | Transport Method | Initial Mass | Final Mass | Delta | Status |
| :--- | :--- | ---: | ---: | ---: | :--- |
| **Deep** | Flux-Form | 137.5130558688 | 137.5130558688 | 0.00000E+00 | **Conservative** |
| **Deep** | Hole-Filling | 137.5130558688 | 137.5130558688 | 0.00000E+00 | **Conservative** |
| **Deep** | Advective (Leaky) | 137.5130558688 | 137.5130558688 | 0.00000E+00 | **Conservative*** |
| **Shallow** | Flux-Form | 137.5130558688 | 137.5130558688 | 0.00000E+00 | **Conservative** |
| **Shallow** | Hole-Filling | 137.5130558688 | 137.5130558688 | 0.00000E+00 | **Conservative** |
| **Shallow** | Advective (Leaky) | 137.5130558688 | 137.5130558688 | 0.00000E+00 | **Conservative*** |

*\*Note: In this test state, the Advective method appears conservative because no tracer values fall below the clipping threshold. Standalone tests (Section 1) confirm its leakage when transport creates "holes".*

## 3. Multi-Tracer Overloading Verification
Verified the ability to handle multiple tracers (3D arrays) in a single call, optimizing memory bandwidth.

| Tracer Index | Initial Mass | Final Mass | Mass Change (Delta) | Status |
| ---: | ---: | ---: | ---: | :--- |
| **Tracer 1** | 10.1971621298 | 10.1971621298 | 0.00000E+00 | **Conservative** |
| **Tracer 2** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | **Conservative** |
| **Tracer 3** | 5.6084391714 | 5.6084391714 | 0.00000E+00 | **Conservative** |

## Conclusion
The structural upgrade successfully introduces conservative transport mechanisms and multi-tracer efficiency. Both the **Flux-Form** and the hybrid **Hole-Filling (Conservative Advection + Hole-Fixer)** methods provide machine-precision mass conservation, eliminating the unphysical mass generation inherent in the previous clipping approach.
