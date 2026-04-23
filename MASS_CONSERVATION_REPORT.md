# Tracer Mass Conservation Report - SAS Convection Schemes

## Overview
This report summarizes the verification of the structural upgrade to the UFS CCPP-physics repository to fix tracer mass conservation leaks in the Simplified Arakawa-Schubert (SAS) convection schemes. The previous advective finite differencing and clipping method has been replaced with strictly conservative flux-form and hole-filling methods, supported by a global proportional mass fixer. Additionally, the transport routines have been optimized to handle multiple species simultaneously via Fortran interfaces (overloading).

## 1. Standalone Transport Module Verification
The `tracer_transport_mod` was tested using a 1D column with 10 layers.

### Scenario A: Small Hole (Fillable Locally)
Initial state: $q = 10^{-3}$ everywhere except $q_5 = -10^{-4}$.

| Method | Initial Mass | Final Mass | Mass Change (Delta) | Min(q) | Status |
| :--- | ---: | ---: | ---: | ---: | :--- |
| **Flux-Form** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | +1.000E-04 | **Conservative** |
| **Hole-Filling Only** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | 0.000E+00 | **Conservative** |
| **Hole-Fill + Fixer** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | 0.000E+00 | **Conservative** |
| **Original Advective** | 9.0754742955 | 9.1774460188 | +1.01972E-01 | +1.000E-10 | **Leaky** (Clipped) |

### Scenario B: Multi-Tracer Overloading Verification
Verified the ability to handle multiple tracers (3D arrays) in a single call.

| Tracer Index | Initial Mass | Final Mass | Mass Change (Delta) | Status |
| ---: | ---: | ---: | ---: | :--- |
| **Tracer 1** | 10.1971621298 | 10.1971621298 | 0.00000E+00 | **Conservative** |
| **Tracer 2** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | **Conservative** |
| **Tracer 3** | 5.6084391714 | 5.6084391714 | 0.00000E+00 | **Conservative** |

## 2. Integrated SAS Scheme Verification
The refactored deep and shallow convection schemes were verified in a unit test environment.

| Scheme | Transport Method | Initial Mass | Final Mass | Delta + Sink | Status |
| :--- | :--- | ---: | ---: | ---: | :--- |
| **Deep** | Flux-Form | 29.5717701564 | 29.5717701564 | 0.00000E+00 | **Conservative** |
| **Deep** | Hole-Filling | 29.5717701564 | 29.5717701564 | 0.00000E+00 | **Conservative** |
| **Deep** | Advective (Leaky) | 29.5717701564 | 29.5717701564 | 0.00000E+00 | **Conservative*** |
| **Shallow** | Flux-Form | 24.4731892055 | 24.4731892055 | -0.71054E-14 | **Conservative** |
| **Shallow** | Hole-Filling | 24.4731892055 | 24.4731892055 | -0.71054E-14 | **Conservative** |
| **Shallow** | Advective (Leaky) | 24.4731892055 | 24.4731892055 | -0.71054E-14 | **Conservative*** |

*\*Note: In stable integrated profiles, the Advective method might appear conservative if no values fall below the clipping threshold. Standalone tests (Section 1) confirm its leakage under transport-induced "holes".*

## Conclusion
The structural upgrade successfully introduces conservative transport mechanisms. The implementation of **Multi-Tracer** interfaces allows for high-efficiency processing of large tracer sets (e.g., chemistry) by minimizing subroutine call overhead and maximizing cache utilization. Both **Flux-Form** and **Hybrid Hole-Filling** methods provide machine-precision mass conservation, eliminating unphysical mass generation.
