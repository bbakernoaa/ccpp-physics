# Tracer Mass Conservation Report - SAS Convection Schemes

## Overview
This report summarizes the verification of the structural upgrade to the UFS CCPP-physics repository to fix tracer mass conservation leaks in the Simplified Arakawa-Schubert (SAS) convection schemes. The previous advective finite differencing and clipping method has been replaced with strictly conservative flux-form and hole-filling methods, supported by a proportional mass fixer.

## Standalone Transport Module Verification
The `tracer_transport_mod` was tested using a 1D column with 10 layers. A consistent initial tracer distribution containing a negative value (a "hole") was used to compare how each method handles non-physical values and maintains total mass.

### Scenario: Fillable Hole
Initial state: $q = 10^{-3}$ everywhere except $q_5 = -10^{-4}$. Total column mass is positive.

| Method | Initial Mass | Final Mass | Mass Change (Delta) | Min(q) | Status |
| :--- | ---: | ---: | ---: | ---: | :--- |
| **Flux-Form** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | +1.000E-04 | **Conservative** |
| **Hole-Filling Only** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | 0.000E+00 | **Conservative** |
| **Hole-Fill + Fixer** | 9.0754742955 | 9.0754742955 | 0.00000E+00 | 0.000E+00 | **Conservative** |
| **Original Advective** | 9.0754742955 | 9.1774460188 | +1.01972E-01 | +1.000E-10 | **Leaky** (Clipped) |

*Note: In this scenario, the Hole-Filling routine alone is perfectly conservative because there is enough mass in the column to fill the deficit locally. The Advective method generates mass by clipping the negative value to a minimum floor.*

## Integrated SAS Scheme Verification
The refactored deep and shallow convection schemes were verified in a unit test environment to ensure correct integration of the new transport module.

### Shallow Convection (`shalcnv.F`)
Tested with `transport_opt = 2` (Hole-Filling Hybrid) and `use_mass_fixer = .true.`.

| Tracer | Initial Mass | Final Mass | Delta |
| :--- | :--- | :--- | :--- |
| **Specific Humidity (q1)** | 24.4731891115 | 24.4731891115 | 0.00000E+00 |

### Deep Convection (`sascnvn.F`)
Tested with `transport_opt = 2` (Hole-Filling Hybrid) and `use_mass_fixer = .true.`.

| Tracer | Initial Mass | Final Mass | Delta + Sink (Precip) |
| :--- | :--- | :--- | :--- |
| **Specific Humidity (q1)** | 29.5717701764 | 29.5717701764 | 0.00000E+00 |

## Conclusion
The structural upgrade successfully introduces conservative transport mechanisms. Both the **Flux-Form** and the hybrid **Hole-Filling (Conservative Advection + Hole-Fixer)** methods provide machine-precision mass conservation. The **Mass Fixer** further ensures that any residual errors or physical sinks (like precipitation) are accounted for proportionally across the column, eliminating the unphysical mass generation caused by the previous clipping method.
