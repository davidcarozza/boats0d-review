% boats0d_main.m
%-----------------------------------------------------------------------------------------
% run script for 0d BOATS model
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Define specific parameters for this run
% Parameters not defined here will be set by boats0d_parameters.m
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% General parameters
 name 	      = 'SIM_EBS';	% name of simulation
 dtts 	      = 15;	% days (converted to seconds in boats0d_integrate.m)
 run_length   = 100;	% years
 iforcing_ESM = 1;      % 1 = annual forcing for a site; 10 = constant forcing with parameters npp and temp75

%-----------------------------------------------------------------------------------------
% Latitude and longitude of site
% lat = 87; lon = 279; % (Humboldt Current)
 lat = 154; lon = 195; % (East Bering Sea)
% lat = 70; lon = 12; % (Benguela)
% lat = 138; lon = 310; % (Newfoundland-Labrador Shelf)
% lat = 65; lon = 245; % (South Pacific Gyre)

%-----------------------------------------------------------------------------------------
% Load model parameters file from best of the Monte-Carlo runs
 load best_MCV3.mat
 ind_best = 1; % overall favoured parameter set
 
%-----------------------------------------------------------------------------------------
% Set names of load (lname_rest) and save (sname_rest) files

% Start from primary-production-dominated model (PP)
 lname_rest   = 'PP';
 sname_rest   = name; 

%-----------------------------------------------------------------------------------------
% Switches and forcing
 idoecon               = 0;         % Economics OFF
 iforcing_catchability = 1; 	    % See boats0d_integrate.m for scenarios
 iforcing_price        = 1;         % See boats0d_integrate.m for scenarios
 qcatch0               = 2*1e-5;    % Catchability

%-----------------------------------------------------------------------------------------
% Initialize parameters
 boats = boats0d_parameters;

%-----------------------------------------------------------------------------------------
% Update input parameters
 boats = boats_change_input(boats,'name',name,'dtts',dtts,'run_length',run_length,'iforcing_ESM',iforcing_ESM, ...
             'lat',lat,'lon',lon, ...
             'lname_rest',lname_rest,'sname_rest',sname_rest, ...
             'idoecon',idoecon,'iforcing_catchability',iforcing_catchability, ...
             'iforcing_price',iforcing_price,'qcatch0',qcatch0, ...                            
             'E_activation_A',best.E_activation_A(ind_best),'kappa_eppley',best.kappa_eppley(ind_best),'Prod_star',best.Prod_star(ind_best), ...
             'te',best.te(ind_best),'ppmr',best.ppmr(ind_best),'b_allo',best.b_allo(ind_best),'E_activation_m',best.E_activation_m(ind_best), ...
             'zeta1',best.zeta1(ind_best),'h_allo',best.h_allo(ind_best),'A00',best.A00(ind_best),'egg_surv',best.egg_surv(ind_best), ...
             'sel_slope',best.sel_slope(ind_best),'sel_pos_scale',best.sel_pos_scale(ind_best), ...
             'npp',1000,'temp75',20);

%-----------------------------------------------------------------------------------------                            
% Set initial conditions and loads restart if needed
 boats = boats0d_initialize(boats);

%-----------------------------------------------------------------------------------------
% Run model
 boats = boats0d_integrate(boats);
 
%-----------------------------------------------------------------------------------------
% Save restart
 boats = boats0d_save_restart(boats);

%-----------------------------------------------------------------------------------------
% Plot base time series and equilibrium distributions of fish and harvest (if needed)
 if (1)
   boats0d_plot_base(boats);
 end

%-----------------------------------------------------------------------------------------
% Remove the time dimension by averaging (if needed)
 if (0)
   boats = boats0d_time_average(boats);
 end

%-----------------------------------------------------------------------------------------
% Add diagnostics
 boats = boats0d_add_diagnostics(boats);

%-----------------------------------------------------------------------------------------
% END OF SCRIPT
