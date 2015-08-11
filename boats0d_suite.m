% boats0d_suite.m
%-----------------------------------------------------------------------------------------
% Run a suite of 0-D simulations
% Either 1 or 2 variables
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Add paths
  addpath('/archive/dcarozza/MATLAB_PATH')
  addpath('/archive/dcarozza/boats_gen')

%-----------------------------------------------------------------------------------------
% Set suite name
 clear Suite;
% Suite.name = 'NPP_vs_NH_HR';
% Suite.name = 'NPP_vs_T_HR';
 Suite.name = 'NPP_vs_Q_HR_2';
%Suite.name = 'EA_sen';
%Suite.name = 'Em_sen';
%Suite.name = 'kE_sen';
%Suite.name = 'PS_sen';
%Suite.name = 'te_sen';
%Suite.name = 'pr_sen';
%Suite.name = 'zeta1_sen';
%Suite.name = 'hallo_sen';
%Suite.name = 'A00_sen';
%Suite.name = 'eggsurv_sen';

%-----------------------------------------------------------------------------------------
% Set basic parameters - used by the whole run suite
% All other parameters will be as in boats0d_parameters.m
 dtts            = 15;           % days (converted to seconds in integration.m
 run_length      = 1000;          % years
 idoecon         = 1;
 qcatch0         = 10e-5;
 saveequilstates = 0;			 % store equilibrium states
 iforcing_ESM    = 10;
 
 model_version   = 'V3';

%-----------------------------------------------------------------------------------------
% Set restart files
% lname_rest   = 'RUN_H';   % name of restart-file to save ('none' for no restart save)
% sname_rest   = 'none';   % name of restart-file to load ('none' for no restart load)

% Start from spinup with no harvest (SPINUP_NOH)
 lname_rest   = 'PP';
 sname_rest   = 'RUN_H';

%-----------------------------------------------------------------------------------------
% Latitude and longitude of site 
 lat = 87; lon = 279; % (Peruvian Upwelling)
% lat = 65; lon = 245; % (South Pacific Gyre)
 %lat = 154; lon = 195; % (East Bering Sea)

%-----------------------------------------------------------------------------------------
% Load model parameters file from best of the Monte-Carlo runs

 model_variant = 'A';

%-------------------------------------------------------------------------------
% Load best parameters from MCV3
 load /archive/dcarozza/DATA/best_MCV3.mat
 ind_best = 1; % overall favoured parameter set

%-----------------------------------------------------------------------------------------
% Set parameter values to analyze
%-----------------------------------------------------------------------------------------


%-----------------------------------------------------------------------------------------
% NPP versus catchability
%-----------------------------------------------------------------------------------------

if (1)

 Suite.params 	= {'npp','qcatch0'};
 
 Suite.npp = 100:100:2000;
 Suite.qcatch0 = (0:1:50)*1e-5;

% Suite.npp = 100:400:2000;
% Suite.qcatch0 = (0:20:100)*1e-5;

end

%-----------------------------------------------------------------------------------------
% NPP versus temperature
%-----------------------------------------------------------------------------------------

if (0)

 Suite.params 	= {'npp','temp75'};
 
 Suite.npp = 50:50:2000;
 Suite.temp75 = -2:2:30;

% Suite.npp = 50:400:2000;
% Suite.temp75 = -2:8:30;

end

%-----------------------------------------------------------------------------------------
% NPP versus T with harvest
%-----------------------------------------------------------------------------------------

if (0)

 idoecon = 1;
 qcatch0 = 20e-5;

 Suite.params 	= {'npp','temp75'};
 
% Suite.npp = 50:50:2000;
% Suite.temp75 = -2:2:30;

 Suite.npp = 50:200:2000;
 Suite.temp75 = -2:4:30;

end 


if (0)

 Suite.params 	= {'npp'};
 
 Suite.npp = 1000:50:1200

end 

if (0)

 idoecon = 1;
 qcatch0 = 10e-5;
 Suite.params 	= {'npp'}; 
 Suite.npp = 1000:50:1200

end

if (0)

  Suite.params 	= {'rep_scale_pos'};
  Suite.rep_scale_pos = 0.1:0.2:2;

end

if (0)

  Suite.params 	= {'rep_scale_slope'};
  Suite.rep_scale_slope = -3:0.5:1;
  Suite.rep_scale_slope = 10.^(Suite.rep_scale_slope);

end

if (0)

  Suite.params 	= {'rep_scale_slope','rep_scale_pos'};
  Suite.rep_scale_slope = -3:0.5:1;
  Suite.rep_scale_slope = 10.^(Suite.rep_scale_slope);
  Suite.rep_scale_pos = 0.1:0.2:2;

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'sel_scale_pos_1'};
  Suite.sel_scale_pos_1 = 0.1:0.2:1.9;

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'sel_scale_slope_1'};
  Suite.sel_scale_slope_1 = -3:0.5:1;
  Suite.sel_scale_slope_1 = 10.^(Suite.sel_scale_slope_1);

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'sel_scale_pos_1','sel_scale_slope_1'};
  Suite.sel_scale_slope_1 = -3:0.5:1;
  Suite.sel_scale_slope_1 = 10.^(Suite.sel_scale_slope_1);
  Suite.sel_scale_pos_1 = 0.1:0.2:1.9;

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'sel_scale_pos_2','sel_scale_slope_2'};
  Suite.sel_scale_slope_2 = -3:0.5:1;
  Suite.sel_scale_slope_2 = 10.^(Suite.sel_scale_slope_2);
  Suite.sel_scale_pos_2 = 0.1:0.2:1.9;

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'sel_scale_pos_3','sel_scale_slope_3'};
  Suite.sel_scale_slope_3 = -3:0.5:1;
  Suite.sel_scale_slope_3 = 10.^(Suite.sel_scale_slope_3);
  Suite.sel_scale_pos_3 = 0.1:0.2:1.9;

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'rep_scale_pos','sel_scale_pos_2'};
  Suite.rep_scale_pos = 0.1:0.2:2;
  Suite.sel_scale_pos_2 = 0.1:0.2:1.9;

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'rep_scale_pos','sel_scale_slope_2'};
  Suite.rep_scale_pos = 0.1:0.2:2;
  Suite.sel_scale_slope_2 = -3:0.5:1;
  Suite.sel_scale_slope_2 = 10.^(Suite.sel_scale_slope_2);

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'rep_scale_slope','sel_scale_slope_2'};
  Suite.rep_scale_slope = -3:0.5:1;
  Suite.rep_scale_slope = 10.^(Suite.rep_scale_slope);
  Suite.sel_scale_slope_2 = -3:0.5:1;
  Suite.sel_scale_slope_2 = 10.^(Suite.sel_scale_slope_2);

end

if (0)

  idoecon = 1;
  qcatch0 = 10e-5;
  Suite.params 	= {'rep_scale_slope','sel_scale_pos_2'};
  Suite.rep_scale_slope = -3:0.5:1;
  Suite.rep_scale_slope = 10.^(Suite.rep_scale_slope);
  Suite.sel_scale_pos_2 = 0.1:0.2:1.9;

end

%-----------------------------------------------------------------------------------------
% Other parameters
%-----------------------------------------------------------------------------------------

% Suite.params 	= {'E_activation_A'};
% Suite.E_activation_A = [0.45-3*0.09:0.09:0.45+3*0.09];

% Suite.params 	= {'E_activation_m'};
% Suite.E_activation_m = [0.45-3*0.09:0.09:0.45+3*0.09];

% Suite.params 	= {'kappa_eppley'};
% Suite.kappa_eppley = [0.0631-3*0.009:0.009:0.0631+3*0.009];
 
% Suite.params 	= {'Prod_star'};
% Suite.Prod_star = [0.37-3*0.1:0.1:0.37+3*0.1];
 
% Suite.params 	= {'te'};
% Suite.te = [0.13-sqrt(3)*0.04:0.01:0.13+sqrt(3)*0.04];

% Suite.params 	= {'ppmr'};
% Suite.ppmr = [1000:1000:10000];

% Suite.params 	= {'b_allo'};
% Suite.b_allo = [0.7-sqrt(3)*0.05:0.02:0.7+sqrt(3)*0.05];

% Suite.params 	= {'zeta1'};
% Suite.zeta1 = [0.55-3*0.57:0.2:0.55+3*0.57];

% Suite.params 	= {'h_allo'};
% Suite.h_allo = [0.54-3*0.09:0.09:0.54+3*0.09];

%Suite.params  = {'A00'};
%Suite.A00      = [4.46-3*0.5:0.5:4.46+3*0.5];

% Suite.params 	= {'egg_surv'};
% Suite.egg_surv = [0.0252-sqrt(3)*0.0143:0.005:0.0252+sqrt(3)*0.0143];

%-----------------------------------------------------------------------------------------
% Set diagnostic parameters
%-----------------------------------------------------------------------------------------

 Suite.additional_processing = 1;	% 1 - does additional input processing before each run (see below)
 Suite.time_average = 1;		% 0 - no time averaging	
                    			% 1 - simple time_average
 Suite.add_diagnostics = 1;		% 0 -x =  no diagnostics
                                % 1 - adds diagnostics
 Suite.rmout = 1;               % removes the single Out cells from final suite
                                % otherwise output could be large
 AddName = 0;                   % 1 to add the parameter names to the Suite name

%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------
% Run suite of experiments
%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------
% Run baseline - for setup purposes

 Suite.base = boats0d_parameters;
              
 Suite.base = boats_change_input(Suite.base,'RunName',Suite.name,'lname_rest',lname_rest, ...
             'model_variant',model_variant,'model_version',model_version,'dtts',dtts, ...
             'iforcing_ESM',iforcing_ESM,'idoecon',idoecon,'run_length',run_length,'lat',lat,'lon',lon, ...
             'E_activation_A',best.E_activation_A(ind_best),'kappa_eppley',best.kappa_eppley(ind_best),'Prod_star',best.Prod_star(ind_best), ...
             'te',best.te(ind_best),'ppmr',best.ppmr(ind_best),'b_allo',best.b_allo(ind_best),'E_activation_m',best.E_activation_m(ind_best), ...
             'zeta1',best.zeta1(ind_best),'h_allo',best.h_allo(ind_best),'A00',best.A00(ind_best),'egg_surv',best.egg_surv(ind_best), ...
             'sel_slope',best.sel_slope(ind_best),'sel_pos_scale',best.sel_pos_scale(ind_best), ...
             'npp',2000,'temp75',10);

 if (idoecon)
 
   Suite.base = boats_change_input(Suite.base,'qcatch0',qcatch0)
 
 end
              
 Suite.base       =  boats0d_initialize(Suite.base);
 Suite.parameters = Suite.base.parameters;
 Suite.nparam     = length(Suite.params);
 Suite.dims       = zeros(1,Suite.nparam);
 Suite.AllParam   = cell(1,Suite.nparam);

 for ip = 1:Suite.nparam
    Suite.dims(ip) = length(eval(['Suite.' Suite.params{ip}]));
    Suite.AllParam{ip} = eval(['Suite.' Suite.params{ip}]);
 end

 Suite.nruns = prod(Suite.dims);

 if length(Suite.dims)>1
    Suite.Out = cell(Suite.dims);
 else
    Suite.Out = cell(1,Suite.dims);
 end

 %----------------------------------------------------------------------------------------
 % Starts main loop
 %----------------------------------------------------------------------------------------
 
 Tsuite = Suite.base;
 runindex = cell(Suite.nparam,1);
 for irun = 1:Suite.nruns
    disp('--------------------------------------------------------');
    disp(['Run number # ' num2str(irun) '/' num2str(Suite.nruns)]);
    Tsuite =  boats_change_input(Tsuite,'RunName',[Suite.name '-' num2str(irun)]);
    [runindex{:}] = ind2sub(Suite.dims,irun);
    for ipar = 1:Suite.nparam
       disp([ Suite.params{ipar} ' - Start ........  ' num2str(Suite.AllParam{ipar}(runindex{ipar}))]);
       % Updates the suite parameters
       Tsuite =  boats_change_input(Tsuite,Suite.params{ipar},Suite.AllParam{ipar}(runindex{ipar}));
    end
    Tsuite = boats0d_initialize(Tsuite);
    %-------------------------------------------------------------------------------------
    if (Suite.additional_processing)
       % Do additional calculations if needed
       % For example substitute parameters that depend on the input
       % Set depth-averaged npp to npp/75 m,  where 75 m is a typical euphotic zone depth 
       Tsuite.parameters.npp_mean = Tsuite.parameters.npp/75;
    end
    %-------------------------------------------------------------------------------------
    Suite.Out{irun} = boats0d_integrate(Tsuite);

    %------------------------------------------------------------------------------------
    % Save equilibrium states (es)

    if (saveequilstates)

      [fish1_equil_states{irun} fish2_equil_states{irun} fish3_equil_states{irun} ...
      harvest1_equil_states{irun} harvest2_equil_states{irun} harvest3_equil_states{irun}] = boats0d_unique_equistates(Suite.Out{irun});

    end

    % Time average of the results
    if Suite.time_average == 1
       Suite.Out{irun} = boats0d_time_average(Suite.Out{irun});
    end
    if Suite.add_diagnostics == 1
       Suite.Out{irun} = boats0d_add_diagnostics(Suite.Out{irun});
    end

 end
 
 Suite = rmfield(Suite,'base');
 %----------------------------------------------------------------------------------------
 if Suite.add_diagnostics==1
    % Process the diagnostics in the suite to matrix form    
    % First create a dummy structure that holds all diagnostics
    diagnames = fieldnames(Suite.Out{1}.diagnostics);
    AllDiag = cell(Suite.nruns,1);
    for indr=1:Suite.nruns
       AllDiag{indr} = Suite.Out{indr}.diagnostics;
    end
    if length(Suite.dims)>1
       AllDiag = reshape(AllDiag,Suite.dims);
    end
    AllDiag = catstruct(1,AllDiag{:});
    %AllDiag = catstruct(AllDiag{:});
    %AllDiag = [AllDiag{:}];
    if length(Suite.dims)>1
      AllDiag = reshape(AllDiag,Suite.dims);
    end
    % Add the diagnostics in matrix form
    % Deal with any additional dimension they might have
    for indn=1:length(diagnames)
       if isnumeric(AllDiag(1).(diagnames{indn}))
          tsize = size(AllDiag(1).(diagnames{indn}));
          tsize(tsize==1) = [];
          ndadd = length(tsize);
          if length(Suite.dims)>1
             Suite.diagnostics.(diagnames{indn}) = reshape([AllDiag.(diagnames{indn})],[tsize Suite.dims]);
             Suite.diagnostics.(diagnames{indn}) = shiftdim(Suite.diagnostics.(diagnames{indn}),ndadd);
          else
             Suite.diagnostics.(diagnames{indn}) = squeeze(reshape([AllDiag.(diagnames{indn})],[tsize 1 Suite.dims]));
             Suite.diagnostics.(diagnames{indn}) = shiftdim(Suite.diagnostics.(diagnames{indn}),ndadd);
          end
       end
    end
 end

 %----------------------------------------------------------------------------------------
 if Suite.rmout == 1
    % Removes all the single runs from the suite
    % to keep the diagnostics only
    Suite = rmfield(Suite,'Out');
 end

 % Remove data_monthly from the output
 Suite.parameters.data_monthly = [];
 
 % Rename the suite
 snewname = ['suite_' Suite.name];
 if AddName==1
    % Create a newname that includes all the parameters
    for indn=1:Suite.nparam
       snewname = [snewname '_' Suite.params{indn}];
    end
 end

%----------------------------------------------------------------------------------------
% Save suite
%----------------------------------------------------------------------------------------

 eval([snewname ' = Suite;']);
 % Save the suite
 Suite.Out = [];
 eval(['save ' snewname '.mat ' snewname ' -v7.3;']);

%----------------------------------------------------------------------------------------
% END OF SCRIPT
