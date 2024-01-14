# vertical_tank_convection

Numerical work in support of vertical natural convection experiments.

This repo contains Firedrake scripts for finding the steady-state flow and temperature profile in a tank of fluid with a temperature difference between the left and the right walls.  All scripts are the same basic code but they are set up with different example flows.

All scripts execute in at most a few core-minutes on a laptop.

Leeds_vertical_convection_MIT.py - reproduces the flow from https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/mit_benchmark/mit_result.html assuming that the flow is steady (which is at least very nearly true).  The script gives a Nusselt number of 


