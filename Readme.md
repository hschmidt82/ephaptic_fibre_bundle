
# Spike Propagation Model for Ephaptic Coupling Effects in Axonal Fibre Bundles.

This is the MATLAB code used in Schmidt, Hahn, Deco and Kn&ouml;sche (PLOS Comp Biol, 2020).

*Run_trials.m* generates the delay data using the *fSpikerun.m* function. The script uses the MATLAB Parallel Toolbox; if you do not have access to it, replace the *parfor*-loop by a *for*-loop. The output of this script is saved in *Results.mat*.

*Plot_results.m* visualises the results saved in *Results.mat*.


