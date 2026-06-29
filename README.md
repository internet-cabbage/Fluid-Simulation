# Fluid-Simulation

## Benchmark info

### Setup:
  - dt = 0.0035
  - tSteps = 10000
  - N = 1000

### Pre-optimisation (yoshida integration, with rudimentary parralelisation):

Total time to run: 23.513921 
Integration time: 23.136364 
Collision time: 0.270424 

Distance optimisation:

If the distance between two bodies is greater than 2.0, then the force calculation is skipped as the force at that distance would be negligible.

Total time to run: 7.907529 
Integration time: 7.565380 
Collision time: 0.255236 
