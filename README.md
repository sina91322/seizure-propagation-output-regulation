# Suppressing Seizure Propagation in a Coupled Neural Mass Model Using Output Regulation

This repository contains the MATLAB and Maple code associated with the manuscript:

**"Suppressing seizure propagation in a coupled neural mass model using output regulation"**

The code implements an observer-based nonlinear output-regulation strategy for suppressing seizure propagation in a bidirectionally coupled two-column Jansen–Rit neural mass model.

The control objective is to regulate the inter-column propagation state associated with pathological transmission from the epileptic cortical column to the initially normal column.

## Repository Structure

### MATLAB code

- `simulate_seizure_propagation.m`  
  Simulates the coupled two-column Jansen–Rit model and reproduces the seizure-propagation scenario presented in Fig. 3.

- `symbolic_regulator_map_A1.m`  
  Implements the steady-state regulator mappings as functions of the excitatory synaptic gain `A1`.

- `design_controller_gain.m`  
  Designs the state-feedback gain using the nominal LQR formulation.

- `K_nominal.mat`  
  Contains the nominal state-feedback gain used in the closed-loop simulations.

- `observer_A1_gains.mat`  
  Contains the observer gains used for state and parameter estimation.

- `simulate_nominal_closed_loop.m`  
  Runs the nominal closed-loop simulation for `A1 = 7` and reproduces the results presented in Fig. 4.

- `simulate_A1_robustness.m`  
  Performs the closed-loop robustness simulations for:
  `A1 = 6.5, 7.0, 7.5`.

- `symbolic_A1_robustness_results.mat`  
  Contains the saved results of the robustness simulations.

- `plot_A1_robustness.m`  
  Generates the robustness plots presented in Fig. 5.

- `verify_controller_stability.m`  
  Performs the common quadratic Lyapunov stability analysis for the closed-loop controller dynamics over the bounded post-transient operating region considered in the manuscript.

- `verify_observer_stability.m`  
  Performs the common quadratic Lyapunov stability analysis for the augmented observer error dynamics.

### Symbolic derivation

- `regulator_equations_A1_symbolic.mw`  
  Maple worksheet used for the symbolic/semi-analytical derivation of the regulator equations and steady-state mappings.

Maple is not required for running the MATLAB simulations because the resulting regulator mappings are implemented directly in the MATLAB code.

## Main Simulation Parameters

The principal simulation parameters used in the manuscript are:

- Epileptic column: `A1 = 7`, `B1 = 22`
- Initially normal column: `A2 = 3.25`, `B2 = 22`
- Inter-column coupling: `k12 = 4000`, `k21 = 1`
- Near-zero propagation-state reference: `epsilon = 1e-3`
- Stochastic input mean: `mu = 90`
- Stochastic input standard deviation: `sigma = 30`
- Simulation duration: `20 s`
- Controller activation time: `5 s`

For the robustness analysis, `A1` is varied over the interval `[6.5, 7.5]`, while the controller and observer gains are kept unchanged.

## Recommended Execution Order

To reproduce the main numerical results, run the MATLAB scripts in the following order:

1. `simulate_seizure_propagation.m`
2. `symbolic_regulator_map_A1.m`
3. `design_controller_gain.m`
4. `simulate_nominal_closed_loop.m`
5. `simulate_A1_robustness.m`
6. `plot_A1_robustness.m`
7. `verify_controller_stability.m`
8. `verify_observer_stability.m`

## Reproducibility

The provided scripts reproduce the principal numerical simulations and stability analyses reported in the manuscript.

The nonlinear simulations use Gaussian stochastic cortical inputs. Therefore, individual stochastic trajectories may vary depending on the random-number initialization.

The robustness simulations investigate `A1 = 6.5, 7.0, 7.5` without retuning the controller or observer gains.

The common quadratic stability analyses use:

- 128 vertices for the closed-loop controller polytope.
- 512 vertices for the augmented observer error-dynamics polytope.

## Requirements

- MATLAB
- Maple (optional; only required for inspecting or reproducing the symbolic derivation)

Additional MATLAB toolbox and optimization requirements are specified by the corresponding scripts.



