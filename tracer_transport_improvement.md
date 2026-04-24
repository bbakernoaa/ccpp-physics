# SAMF Tracer Transport Mass Conservation Improvement Report

## Issue Description
The SAMF (Scale-Aware Mass Flux) deep and shallow convection schemes exhibited systematic mass conservation errors for tracers when negative concentrations occurred at the top model layer (`km`).

The root cause was identified in the proportional mass fixer logic. Loops responsible for calculating column-integrated tracer mass and redistributing it used an upper bound of `km-1` instead of `km`. This excluded the top layer from redistribution, causing uncorrected negative values to be subsequently clipped to zero, leading to a net mass gain.

## Verification
A standalone logic test (`physics/CONV/SAMF/tests/test_fixer_cons.f90`) was used to verify the fix.

### Baseline (km-1 loops):
-   **Conservation Error:** +3.059e-05

### Fixed (km loops):
-   **Conservation Error:** +1.734e-18 (Machine precision)

## Changes
1.  **`physics/CONV/SAMF/samfdeepcnv.f`** & **`samfshalcnv.f`**: Extended loop bounds in the mass fixer sections to include the top layer.
2.  **`physics/CONV/progsigma_calc.f90`** & **`progomega_calc.f90`**: Improved portability by using consistent `kind_phys` in intrinsic function calls.
3.  **`physics/CONV/SAMF/tests/`**: Added unit tests for tracer mass conservation.
