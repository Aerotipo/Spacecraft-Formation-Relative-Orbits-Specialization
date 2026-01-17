### *commit 19*
quiz 4 "general equations of motion", exercise 2. The propagators for: (i) inertial ODE integrator, (ii) Kepler Equations, (iii) Relative EOM ODE integrator seem to be working all fine. The integrators (i) and (ii) lead to similar results within 6 digits precision.
Nevertheless, when mapping from inertial (ECI) frame to Hill-frame, there are discrepancies between [(i)-(ii)] and (iii) on the deputy inertial position and velocity. The discrepancy is proportional to the integration time. This subject requires further investigation.
Recommendation: investigate the mathematical formulation of the general relative EOM to see if there is something missing in the code, relative to the formulation.

### *commit 26*
<p align="left">quiz 16 "OE Difference Solutions". The problem was solved using 3 methods:
- Method 1: Non-linear mapping between orbit elements and Hill state vector X
- Method 2: Using the EO diffrence solution for general elliptic chief orbits
- Method 3: Linearized mapping about chief orbit elements using matrix [A(oe)_chief]
The submission was done using the results from Method 1. Nevertheless, discrepancies between methods 2 & 3 with respect to Method 1 (and to each other) were observed. The equations were checked thoroughly and confirmed to be similar to those presented in the slides. More investigation is needed to explain/correct this discrepancy. As it stands right now, I don't trust the results of method 2 and method 3. For method 3 the code for constructing the matrix [A(oe__chief] was validated in a previous quiz.</p>
