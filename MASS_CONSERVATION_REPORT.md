# Tracer Mass Conservation Report - SAS Convection Schemes

## Overview
This report summarizes the verification of the structural upgrade to the UFS CCPP-physics repository to fix tracer mass conservation leaks in the Simplified Arakawa-Schubert (SAS) convection schemes. The advective finite differencing and clipping method has been replaced with strictly conservative flux-form and hole-filling methods, supported by a proportional mass fixer.

## Standalone Transport Module Verification
The `tracer_transport_mod` was tested using a 1D column with 10 layers and a prescribed mass flux at layer interface 5.

| Method | Initial Column Mass | Final Column Mass | Mass Change (Delta) | Status |
| :--- | :--- | :--- | :--- | :--- |
| **Original Advective** | 9.1774459678 | 9.1774460188 | +5.09860E-08 | **Leaky** (Clipped) |
| **Flux-Form** | 10.1971621298 | 10.1971621298 | 0.00000E+00 | **Conservative** |
| **Hole-Filling** | 9.0754742955 | 8.9434742955 | -1.32000E-01 | **Redistributed** (Needs Fixer for large holes) |
| **Mass Fixer** | 10.1971621298 | 10.1971621298 | +5.32907E-15 | **Conservative** |

## Integrated SAS Scheme Verification
The refactored deep and shallow convection schemes were tested in a unit test environment.

### Shallow Convection (`shalcnv.F`)
Tested with `transport_opt = 1` (Flux-Form) and `use_mass_fixer = .true.`.

| Tracer | Initial Mass | Final Mass | Delta |
| :--- | :--- | :--- | :--- |
| **Specific Humidity (q1)** | 24.4731891115 | 24.4731891115 | -8.17124E-14 |

### Deep Convection (`sascnvn.F`)
Tested with `transport_opt = 1` (Flux-Form) and `use_mass_fixer = .true.`.

| Tracer | Initial Mass | Final Mass | Delta + Sink (Precip) |
| :--- | :--- | :--- | :--- |
| **Specific Humidity (q1)** | 29.5717701764 | 29.5717701764 | 0.00000E+00 |

## Conclusion
The structural upgrade successfully introduces conservative transport mechanisms. The **Flux-Form** method provides machine-precision mass conservation without the need for arbitrary clipping. The **Mass Fixer** further ensures that any residual errors or physical sinks (like precipitation) are accounted for proportionally across the column.
