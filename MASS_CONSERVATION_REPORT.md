# Tracer Mass Conservation Report - SAS Convection Schemes

## Overview
This report summarizes the verification of the structural upgrade to the UFS CCPP-physics repository to fix tracer mass conservation leaks in the Simplified Arakawa-Schubert (SAS) convection schemes. The previous advective finite differencing and clipping method has been replaced with strictly conservative flux-form and hole-filling methods, supported by a global proportional mass fixer and optimized for multi-tracer efficiency.

## 1. Comparative Analysis: "The Dry Cliff" Scenario
To demonstrate the fundamental weakness of the legacy scheme and the robustness of the new methods, we used an extreme "Dry Cliff" profile (high subsidence pushing dry air into a thin moist layer). This scenario intentionally triggers the `max(q, 1e-10)` clipping mechanism in the legacy advective transport.

| Test Name | Initial Mass | Final Mass | Mass Change (Delta) | Status |
| :--- | ---: | ---: | ---: | :--- |
| **Legacy** | 60.1632565657 | 84.1632565657 | +24.00000 | **Leaky** (Clipped) |
| **Legacy (clipping) + fixer** | 60.1632565657 | 60.1632565657 | -1.42109E-14 | **Conservative** |
| **Hole filler without fixer** | 60.1632565657 | 60.1632565657 | 0.00000E+00 | **Conservative** |
| **Hole filler with fixer** | 60.1632565657 | 60.1632565657 | 0.00000E+00 | **Conservative** |
| **Flux without fixer** | 60.1632565657 | 60.1632565657 | 0.00000E+00 | **Conservative** |
| **Flux with fixer** | 60.1632565657 | 60.1632565657 | 0.00000E+00 | **Conservative** |

*Note: The **Legacy** method generates massive artificial mass (+24.0 in this scenario) because it cannot safely resolve the sharp moisture gradient, leading to negative overshoots that are then clipped. The new **Hole-Filler** and **Flux-Form** methods are conservative by design, as they correctly redistribute or limit fluxes to preserve total mass.*

## 2. Integrated Verification (Loaded Gun Profile)
The refactored schemes were verified using the "Loaded Gun" severe convection profile to ensure stability and conservation under realistic atmospheric instability.

| Scheme | Transport Method | Initial Mass | Final Mass | Delta | Status |
| :--- | :--- | ---: | ---: | ---: | :--- |
| **Deep** | Hole-Filling (Hybrid) | 137.5130558688 | 137.5130558688 | 0.00000E+00 | **Conservative** |
| **Shallow** | Hole-Filling (Hybrid) | 137.5130558688 | 137.5130558688 | 0.00000E+00 | **Conservative** |

## 3. Multi-Tracer Overloading
Verified that 3D tracer arrays are correctly processed with machine-precision conservation, confirming the efficiency and accuracy of the vectorized implementation.

| Tracer Index | Initial Mass | Final Mass | Mass Change (Delta) | Status |
| ---: | ---: | ---: | ---: | :--- |
| **Tracer 1** | 10.1971621298 | 10.1971621298 | 0.00000E+00 | **Conservative** |
| **Tracer 2** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | **Conservative** |
| **Tracer 3** | 5.6084391714 | 5.6084391714 | 0.00000E+00 | **Conservative** |

## Conclusion
The structural upgrade successfully eliminates mass leaks in the SAS convection schemes. The implementation provides:
1. **Physical Accuracy**: Strict mass conservation across all transport scenarios.
2. **Computational Efficiency**: Overloaded interfaces for high-bandwidth multi-tracer processing.
3. **Robustness**: Advanced hole-filling that prevents non-physical negative concentrations without introducing spurious mass.
