function boats = boats0d_parameters(varargin)

%-----------------------------------------------------------------------------------------
% boats0d_parameters.m
% Define the standard input parameters
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Conversion factors

 A.sperd        = 3600*24;					% seconds per day
 A.spery        = A.sperd*360;				% seconds per year
 A.gC_2_wetB    = 10;						% grams of wet fish biomass per gram of fish carbon
 A.mmolC_2_wetB = (12*A.gC_2_wetB)/1000;    % grams of wet fish biomass per mmol of fish carbon
 A.epsln        = 1e-50;					% small epsilon

%-------------------------------------------------------------------------------------
% Time integration parameters

 A.dtts         = 30;						% days per timestep
 A.run_length   = 5;						% duration of simulation (years)

%-----------------------------------------------------------------------------------------
% Forcing

 % Earth System Model forcing
 % 1==annual climatology, 2==monthly climatology, 3==time varying, 4==time varying no climate change (for CM2.6)
 A.iforcing_ESM          = 1;
 % 1==constant; Numerous other scenarios
 A.iforcing_catchability = 1;
 % Price forcing
 % Standard price change is 1950 price previous to 1950, SAUP price from 1950 to 2006, and change after 2006
 % 1==constant price; 2==standard and constant after 2006; 3==standard and increasing after 2006
 % 4==standard and decreasing after 2006
 A.iforcing_price        = 1;
 A.cap_npp               = 10000; 			% limit on npp (m mol C m-2 d-)
 A.start_year            = 0;
 A.start_year_price      = 0;

%-----------------------------------------------------------------------------------------
% Model choice indices
 A.idoecon               = 0;  				% 1==do economics

%-----------------------------------------------------------------------------------------
% Choose site (latitude and longitude coordinates)

% South Pacific Gyre
 A.lat          = 60;
 A.lon          = 250;
% EEZ Peru
% lat =  80;
% lon = 281;
% East Bering Sea
% lat = 145;
% lon = 199;

%----------------------------------------------------------------------------------------- 
% load data_monthly input data
 data_monthly_exist = exist('data_monthly','var');
 if (~data_monthly_exist) 
   load data_monthly-review.mat
 end

 A.data_monthly = data_monthly;
 [A.nmonths,A.nlat,A.nlon] = size(data_monthly.npp);

%-----------------------------------------------------------------------------------------
% Mass class structure (bin boundaries in logarithmic space)
 A.nfmass 	= 50;							% number of fish mass classes
 A.fmass_0 	= 10;							% initial mass class (g)
 A.fmass_e 	= 1e5;							% final mass class (g)
 A.ifmbound     = 1:1:(A.nfmass+1);       	% size class number bounds: 1,2,3, ...
 A.fmbound      = A.fmass_0 .* (A.fmass_e/A.fmass_0).^( (A.ifmbound-1)./(A.nfmass));
 for indm=1:A.nfmass
     A.fmass(indm) = (A.fmbound(indm)*A.fmbound(indm+1))^(0.5);
 end
 A.delfm        = diff(A.fmbound);			% width of mass classes (g)
 
%-----------------------------------------------------------------------------------------
% Group structure
 A.minf         = [0.01*(30/0.95)^3 0.01*(90/0.95)^3 1e5];	% asymptotic mass
 A.eta_alpha    = 0.25;										% mass at maturity as fraction of asymptotic mass 
 A.malpha       = A.eta_alpha*A.minf;						% maturity mass
 A.nfish        = length(A.minf);							% number of fish groups

%-----------------------------------------------------------------------------------------
% Group structure (2-d arrays)
% Mass of mass class (g)
 A.fmass_2d     = repmat(A.fmass,[A.nfish 1]);
% Asymptotic mass (g)
 A.minf_2d      = repmat(A.minf',[1 A.nfmass]);
% Maturity mass (g)
 A.malpha_2d    = repmat(A.malpha',[1 A.nfmass]);
% Width of mass classes (g)
 A.delfm_2d     = repmat(A.delfm,[A.nfish 1]);
% Width of mass classes from 2nd to final mass class (g)
 A.delfm_2end_2d = repmat(A.delfm(2:end),[A.nfish 1]);
 
%-----------------------------------------------------------------------------------------
% Temperature dependence
 A.E_activation_A = 0.45;        % Activation energy of metabolism (growth A) (eV) (Savage et al., 2004)
 A.E_activation_m = 0.45;        % Activation energy of metabolism (mortality) (eV) (Savage et al., 2004)
 A.k_Boltzmann    = 8.617e-05;   % Boltzmann Constant (eV K-1)
 A.temp_ref_A     = 10 + 273.15; % Reference temperature (K) (Andersen and Beyer, 2013, p. 18)

%-----------------------------------------------------------------------------------------
% Primary production
 A.kappa_eppley = 0.063;		    % Eppley constant (degC-1)
 A.Prod_star    = 0.37; 			    % Pivotal primary production (m mol C m-3 d-1)
 A.C_2_K        = 273.15;		    % deg C to Kelvin
 A.mc_phy_l     = 5.6234132519e-06; % mass of typical large phytoplankton (g)
 A.mc_phy_s     = 5.6234132519e-15; % mass of typical small phytoplankton (g)

%-----------------------------------------------------------------------------------------
% Ecology
 A.te          = 0.125;						% trophic efficiency
 A.ppmr        = 5000;						% predator to prey mass ratio
 A.tro_sca     = log10(A.te)/log10(A.ppmr);	% trophic scaling
 A.b_allo      = 0.66;						% allometric scaling
 A.zeta1       = 0.57;						% constant mortality scaling
 A.h_allo      = 0.5;						% mass scaling of mortality
 A.eff_a       = 0.8;						% efficiency of activity (Andersen and Beyer, 2013, p. 4)
 A.A00         = 4.46;						% allometric growth rate (Andersen and Beyer, 2013, p. 4)

%-----------------------------------------------------------------------------------------
% Reproduction
 A.m_egg       = 5.2e-4;                    % egg mass (g)
 A.frac_fem    = 0.5;	                    % fraction of individuals that allocate energy to reproduction (females)
 A.egg_surv    = 0.01;	                    % egg survival
 A.rep_slope   = 5;		                    % slope parameter of sigmoidal allocation to reproduction function
 A.rep_pos     = 1;		 				    % position parameter of sigmoidal allocation to reproduction function as fraction of malpha

%----------------------------------------------------------------------------------
% Masks
% NaN where mass is > asymptotic mass
 A.mask_notexist_2d = (repmat(A.fmass,[A.nfish 1]) > repmat(A.minf,[A.nfmass 1])');

%-----------------------------------------------------------------------------------------
% Economic parameters
 A.landedvalue_global = 8.4233e+10; 							 % SAUP 1990-2006 average ($)
 A.yield_global       = 7.9963e+13;								 % SAUP 1990-2006 average (g)
 A.price_global       = A.landedvalue_global/A.yield_global;     % Global price ($ g-1)
 A.cost_global        = A.landedvalue_global; 					 % Global total cost ($) Assume C = R (matches Lam et al., 2011)
 A.effort_global      = 14.6229e9;  							 % Global effort (W) (Watson et al. (2012) 1990-2006 average)

 A.cost_effort_0      = A.cost_global/(A.effort_global*A.spery); % Cost per unit effort ($ W-1 s-1)
 A.k_e                = 1e-6; 									 % Fleet dynamic parameter (W $-1 s-1) 1 W of effort per dollar of net revenue per second
 A.sel_pos_1          = 1;										 % Selectivity position shift 1
 A.sel_pos_2          = 0.5;									 % Selectivity position shift 2
 A.sel_pos_3          = 0.25;									 % Selectivity position shift 3
 A.sel_pos_scale      = 1;                                       % Selectivity position scale
 A.sel_slope          = 18;										 % Selectivity slope

 A.harvest_start      = 0;										 % Year of starting harvest [y]
 A.qcatch0            = 0; 										 % Base catchability
 A.price_0            = A.price_global; 						 % Base price (constant)

%-----------------------------------------------------------------------------------------
% Simulation parameters (for CMIP5 run)
 A.model         = 'data_monthly'; 				% Model
 A.sim_type      = 'SPINUP';					% Simulation type
 A.simulation    = 'Cm';						% Simulation (forcing type)
 A.model_version = 'V2';						% Model version	
 sim_variant     = 'MB';						% Simulation variant (parameter set)

 A.lname_rest    = 'none'; 						% name of restart-file to save ('none' for no restart save)
 A.sname_rest    = 'none'; 						% name of restart-file to load ('none' for no restart load)

 A.rc_q_ini      = 0;							% Discount rate of catchability to determine initial value
 A.q_discount_y  = 0;							% Number of years to discount catchability

 A.test          = 0;							% Test (0 == no test, 1 == test)

%-----------------------------------------------------------------------------------------
% Primary production and temperature forcing for constant forcing simulation
 A.npp_mean    = 2.0; 							% mmolC/m3/d
 A.npp         = 85.0; 						    % mmolC/m2/d
 A.temp75      = 20.0;							% degC
 A.temp400     = 14.5;							% degC

%-----------------------------------------------------------------------------------------
% Parse required variables, substituting defaults where necessary
 boats.parameters = parse_pv_pairs(A, varargin);

%-----------------------------------------------------------------------------------------
% END OF SCRIPT
