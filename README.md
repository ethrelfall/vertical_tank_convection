# vertical_tank_convection

Numerical work in support of vertical natural convection experiments.

This repo contains Firedrake scripts for finding the steady-state flow and temperature profile in a tank of fluid with a temperature difference between the left and the right walls.  All scripts are the same basic code but they are set up with different example flows.

All scripts execute in at most a few core-minutes on a laptop.

Note on units: the scripts use units scaled to the free-fall velocity, $u_0 = \sqrt{g \beta  \Delta T}$, (see Appendix B.2 of M6c.3 ExCALIBUR report), so that

$$ 
\begin{align}
\dot{u} + u \cdot \nabla u &= - \nabla p + \sqrt{\frac{Pr}{Ra}} \nabla^2 u + T \hat{y};\\
\dot{T} + u \cdot \nabla T &= \frac{1}{\sqrt{Ra Pr}} \nabla^2 T;\\
\nabla \cdot u &= 0.
\end{align}
$$

In these units, I believe the heat flux is given by $F = \frac{1}{\sqrt{Pr}} \nabla T$.  This agrees with the code outputs (for steady-state solutions, the flux into the LHS must equal the flux out of the RHS).

In all cases, the Nusselt number is computed as the ratio of heat flux across the slow with flow to the same quantity in the purely conducting case.

**Leeds_vertical_convection_MIT.py** - reproduces the flow from https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/mit_benchmark/mit_result.html assuming that the flow is steady (which is at least very nearly true).  The script gives a Nusselt number of $4.577$, which compares well the the values in the reference.  Note that in order to obtain the solution at the values of Rayleigh and Prandtl numbers required by the benchmark, it is necessary to use continuation i.e. do multiple runs in parameter space to approach the desired values.

**Leeds_vertical_convection_water_constPr.py** - treats the case of a water-filled slot with a constant Prandtl number ($Pr=7.0$) corresponding to water at approximately room temperature.  The Nusselt number is 1.731.  This example does not require continuation.

**Leeds_vertical_convection_water_varyPr.py** - treats the case of a water-filled slot with a Prandtl number varying linearly with temperature ($Pr=7.0$ on the cold side and $Pr=2.0$ on the hot side).  The Nusselt number is 1.736 - little different to the constant $Pr$ case.  This example does not require continuation.  Note that with $Pr$ varying in this way, the computation becomes less accurate, as is seen in discrepancy between flux in and flux out.  The discrepancy can be much reduced by increasing element order.  Note that the no-flow case has analytic solution $T=\frac{7}{5} - \frac{1}{5} \left ( \sqrt{2} + (\sqrt{7}-\sqrt{2}) x \right )^2$ and the amplitude of the fluxes is $\frac{2}{5}(\sqrt{7}-\sqrt{2}) \approx 0.492615$.

**Leeds_vertical_convection_vortices.py** - this is included as a fun example rather than a physical demonstration.  The Rayleigh number is $10^4$ and the Prandtl number is $10^{-1.4}$ and for these parameters a well-defined pattern of secondary flow vortices is seen; the author's (quite possibly flawed) thinking was to lower the Prandtl number in order to widen the thermal boundary layer and thus get the flows associated to the two long edges to interact.  It is true that similar secondary flows are seen in Elder's tank experiments, but the author has not found it easy to reproduce these flows using Firedrake.  This example uses continuation to obtain the solution.  For reference, a value of $1.676$ is obtained for the Nusselt number.







