boats0d-review

BOATS-1.0: The Bioeconomic Open Access Trophic Size-spectrum model. A
bioenergetically-constrained coupled fisheries-economics model for global
studies of harvesting and climate change.

Version
-------

Zero-dimensional version (for a single patch of ocean / a single site) submitted
for review to Geoscientific Model Development.

Authors
-------

David A. Carozza  (david.carozza@gmail.com; corresponding author)
Daniele Bianchi   (danbian@uw.edu)
Eric D. Galbraith (eric.galbraith@mcgill.ca)

Summary
-------

BOATS is written in MATLAB version R2012a. Here we provide a list of the
functions and forcing data required to run BOATS. BOATS is implemented using a
single MATLAB structure, named boats, that stores the model parameters, initial
conditions, output, and diagnostics. The boats structure is passed among the
various subroutines in order to set the model parameters, initialize the model
from a restart state, integrate the model through time, save the output and a
restart file, and make plots of model variables.

Usage
-----

Run the boats0d_main.m script in matlab by typing boats0d_main:

boats0d_main

To run the model at another ocean site, change the lat and lon parameters. For
example, to run the model at a site in the Netfoundland-Labrador Shelf Large
Marine Ecosystem, change the lat and lon parameters in boats0d_main as:

lat = 138; lon = 310;

To modify a parameter, use the boats_change_input.m function. For example, to
apply a trophic efficiency of 0.15, add the following code after the first call
to boats_change_input in boats0d_main:

boats = boats_change_input(boats,'te',0.15);

Files
-----

boats0d_main.m

BOATS run script. Includes typical workflow from initialization
to output and diagnostics. Calls functions that operate on the structure boats

boats0d_parameters.m

Set the standard input parameters for boats structure.
Parameters are stored in the field boats.parameters

boats_change_input.m

Change a parameter value.

boats0d_initialize.m

Load restart file and initialize boats structure.

boats0d_integrate.m

Integrate BOATS in time and saves variables to boats structure.

boats0d_save_restart.m

Save final state of the fish biomass spectra and effort
(if harvest on) in the boats.restart field, and save to a .mat restart file.

boats0d_plot_base.m

Plot fish and harvest (if harvest on) time series and equilibrium values of
variables.

boats0d_time_average.m

Remove the time dimension for a run by averaging over a
given time index range using generalmean.m - default is last time step.

boats0d_add_diagnostics.m

Adds diagnostics to boats.diagnostics. These are usually scalar or small
vector diagnostics that can be easily plotted in suites of runs

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
boats0d_add_diagnostics.m, boats0d_parameters.m, boats0d_time_average.m

best_MCV2_P_S.mat

Structure of biological parameters from best Monte Carlo
simulation V2, Pearson and Spearman

catchability_forcing.mat

Structure of catchability forcing scenarios

price_forcing.mat

Structure of price forcing scenarios

Licensing
---------

Copyright 2015 David Anthony Carozza