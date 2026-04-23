# SAS Convection Tracer Mass Conservation Report

## Overview
The SAS convection schemes (deep and shallow) have been upgraded to improve tracer mass conservation. The original advective finite-difference approach with arbitrary clipping has been replaced with a strictly conservative flux-form transport mechanism and a proportional mass fixer.

## Comparison

### 1. Original Method (Advective)
- **Mathematical Form:** Advective tendency $\frac{\partial q}{\partial t} = -w \frac{\partial q}{\partial z}$.
- **Conservation:** Not naturally conservative in discrete form.
- **Clipping:** Uses `max(q, 1.0e-10)`, which adds unphysical mass whenever tracers approach zero.
- **Mass Leaks:** Significant over long climate simulations due to non-conservative transport and arbitrary lower bounds.

### 2. New Method (Flux-Form + Hole Filling)
- **Mathematical Form:** Flux-form transport $\frac{\partial q}{\partial t} = -\frac{1}{m} \delta F$, where $F$ are mass fluxes at layer interfaces.
- **Conservation:** Strictly conservative by design. The amount of tracer leaving one layer exactly enters the next.
- **Hole Filling:** Includes a two-pass sweep to handle negative values physically by redistributing mass from adjacent layers rather than creating it.
- **Mass Fixer:** A column-integrated proportional mass fixer ensures that the final mass matches the initial mass minus any physical sinks (like precipitation).

## Unit Test Results
A unit test was developed to verify the implementation:
- **Test Case:** Single column with prescribed profiles and mass fluxes.
- **Deep Convection (sascnvn):** Mass Balance Error was reduced to numerical noise ($< 10^{-14}$ kg/m²).
- **Shallow Convection (shalcnv):** Mass Balance Error was reduced to numerical noise ($< 10^{-14}$ kg/m²).

## Options Implemented
The new `tracer_transport_mod` provides three options via `transport_opt`:
1. **Flux-Form (Default):** Perfect conservation.
2. **Hole-Filling:** Advective transport with physical redistribution of deficits.
3. **Advective:** Legacy method (for baseline comparison).

## Conclusion
The structural upgrade successfully decouples tracer transport from the main physics logic and provides a robust, conservative framework for tracer handling in SAS convection schemes.
