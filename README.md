boats0d-review
--------------

BOATS-1.0: The BiOeconomic mArine Trophic Size-spectrum model. A
bioenergetically-constrained coupled fisheries-economics model for global
studies of harvesting and climate change.

Version
-------

Zero-dimensional version (for a single patch of ocean / a single site) submitted
for review to Geoscientific Model Development.

Source
------

https://zenodo.org/badge/latestdoi/12481/davidcarozza/boats0d-review

http://dx.doi.org/10.5281/zenodo.27700


Authors
-------

David A. Carozza  (david.carozza@gmail.com; corresponding author)

Daniele Bianchi   (dbianchi@atmos.ucla.edu)

Eric D. Galbraith (eric.galbraith@mcgill.ca)

Summary
-------

BOATS is written in MATLAB version R2012a. Here we provide the
script, functions, and forcing data required to run BOATS. BOATS is implemented using a
single MATLAB structure, named boats, that stores the model parameters, initial
conditions, output, and diagnostics. The boats structure is passed among the
various functions in order to set the model parameters, initialize the model
from a restart state, integrate the model through time, save the output and a
restart file, and make plots of model variables.

Usage and examples
------------------

Run the BOATS model with the run script boats0d_main.m in MATLAB as:

boats0d_main

To change the ocean site at which the model is run, change the lat and lon parameters. For
example, to run the model at a site in the Newfoundland-Labrador Shelf Large
Marine Ecosystem, change the lat and lon parameters in boats0d_main.m as:

lat = 138; lon = 310;

To modify a parameter, use the boats_change_input.m function. For example, to
apply a trophic efficiency of 0.15, add the following code after the first call
to boats_change_input.m in boats0d_main.m:

boats = boats_change_input(boats,'te',0.15);

Files
-----

boats0d_main.m

BOATS run script. Includes typical workflow from initialization, model integration,
to output, diagnostics, and plots. Calls functions that operate on the structure boats.

boats0d_parameters.m

Set the standard input parameters for boats structure. Parameters are stored
in the boats.parameters field.

boats_change_input.m

Change parameter values.

boats0d_initialize.m

Load restart file and initialize boats structure.

boats0d_integrate.m

Integrate BOATS in time and save variables to boats structure.

boats0d_save_restart.m

Save final state of the fish biomass spectra (and effort) in the 
boats.restart field, and save to a restart file.

boats0d_plot_base.m

Plot fish (and harvest) time series and equilibrium values of variables.

boats0d_time_average.m

Remove the time dimension for a run by averaging over a
given time index range using generalmean.m - default is to use the last time step.

boats0d_add_diagnostics.m

Add diagnostics to boats.diagnostics. These are usually scalar or small
vector diagnostics that can be easily plotted in suites of runs.

sigmoid_And_length.m

Calculate values of a sigmoid function (Andersen) with input values in
terms of organism length.

sigmoid_And_mass.m

Calculate values of a sigmoid function (Andersen) with input values in
terms of organism mass.

generalmean.m

Calculate the temporal mean of a variable. Used by boats0d_time_average.m

parse_pv_pairs.m

Parses sets of property value pairs, allows defaults. Used by
boats0d_add_diagnostics.m, boats0d_parameters.m, and boats0d_time_average.m

best_MCV3.mat

MATLAB structure of biological parameters from the 5 best Monte Carlo
simulations of the Monte Carlo suite V3 (vectorized, i.e. only Large Marine
Ecosystems, version 3).

catchability_forcing.mat

MATLAB structure of catchability forcing scenarios.

price_forcing.mat

MATLAB structure of price forcing scenarios

data_monthly-review.mat

MATLAB structure of net primary production (NPP) and temperature (T) forcing data.

restart_sim_EBS.mat

MATLAB structure of restart file for spinup run with East Bering Sea NPP and T.

Licensing
---------

Copyright 2015 David Anthony Carozza
