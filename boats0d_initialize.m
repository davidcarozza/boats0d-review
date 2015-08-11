 function boats = boats0d_initialize(boats)

%-----------------------------------------------------------------------------------------
% boats0d_initialize.m
% Set forcing datasets and convert units (to m molC m-2 s-1)
% Set initial dfish and effort (if needed)
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Set parameters
%-----------------------------------------------------------------------------------------

 model   = boats.parameters.model;
 sperd   = boats.parameters.sperd;
 idoecon = boats.parameters.idoecon;
 nfish   = boats.parameters.nfish;
 nfmass  = boats.parameters.nfmass;
 cap_npp = boats.parameters.cap_npp;
 cap_npp = cap_npp * (1/sperd) * (1/12); % convert from mg C m-2 d-1 to mmol C m-2 s-1
 
 lat     = boats.parameters.lat;
 lon     = boats.parameters.lon;

%-----------------------------------------------------------------------------------------
% Set forcing structure (temperature and primary production)
% Convert units (to m molC m-2 s-1) 
%-----------------------------------------------------------------------------------------
% npp              m molC m-2 s-1
% npp_ed           m molC m-3 d-1 (Average over euphotic zone depth, needed for frac_lg calculation)
% temp_phyto       deg C
% temp_fish        deg K
% mclim_intpp      monthly climatology of integrated primary production
% intpp            time-varying integrated primary production
% mclim_thetao75_C monthly climatology of 75m meter average temperature (deg C)
% mclim_thetao75_K monthly climatology of 75m meter average temperature (K)
%-----------------------------------------------------------------------------------------

% Set data_monthly
 data_monthly = boats.parameters.data_monthly;
 
 %-----------------------------------------------------------------------------------------
% Set forcing with observations (data_monthly) to ESM_forcing
%-----------------------------------------------------------------------------------------

 if (strcmp(model,'data_monthly'))

   % store SEAWIFS NPP in ESM_forcing and convert from mmol C m-2 d-1 to mmol C m-2 s-1
   ESM_forcing.mclim_intpp            = squeeze(data_monthly.npp(:,lat,lon)) * (1/sperd);
   % limit npp to cap_npp
   ESM_forcing.mclim_intpp            = min(cap_npp,ESM_forcing.mclim_intpp);
   % NPP averaged over euphotic zone depth (ed) for frac_lg calculation
   ESM_forcing.mclim_intpp_ed         = ESM_forcing.mclim_intpp * (sperd / 75); % mmolC m-2 s-1 to mmolC m-3 d-1
   % store WOA temperature in ESM_forcing
   ESM_forcing.mclim_thetao75_C       = squeeze(data_monthly.temp75(:,lat,lon));
   ESM_forcing.mclim_thetao75_K       = squeeze(data_monthly.temp75(:,lat,lon)) + boats.parameters.C_2_K;

%-----------------------------------------------------------------------------------------
% Set forcing with climate model output to ESM_forcing
%-----------------------------------------------------------------------------------------
       
 else
 
   % directory where forcing files are saved
   dir_forcing = '/archive/dcarozza/DATA/CMIP5/';
   % load forcing data structure
   load([dir_forcing model '.mat']);
   % rename model forcing structure
   ESM_forcing = eval(model);
   % clear original model forcing structure
   eval(['clear ' model])
   % convert units from molC m-2 s-1 to mmolC m-2 s-1
   ESM_forcing.mclim_intpp            = (1000) * squeeze(ESM_forcing.mclim_intpp(lat,lon,:));
   ESM_forcing.intpp                  = (1000) * squeeze(ESM_forcing.intpp(lat,lon,:));
   % mclim_intpp
   % limit to cap_npp
   ESM_forcing.mclim_intpp            = min(cap_npp,ESM_forcing.mclim_intpp);
   % intpp
   % limit to cap_npp and use mask_land to keep land cells as NaNs
   ESM_forcing.intpp                  = min(cap_npp,ESM_forcing.intpp);
   % mclim_intpp_ed and intpp_ed
   % convert units from mmolC m-2 s-1 to mmolC m-3 d-1
   ESM_forcing.mclim_intpp_ed         = ESM_forcing.mclim_intpp * (sperd / 75);
   ESM_forcing.intpp_ed               = ESM_forcing.intpp * (sperd / 75);
   % mclim_thetao75_K and mclim_thetao75_C
   ESM_forcing.mclim_thetao75_K       = squeeze(ESM_forcing.mclim_thetao75(lat,lon,:));
   ESM_forcing.mclim_thetao75_C       = squeeze(ESM_forcing.mclim_thetao75(lat,lon,:)) - boats.parameters.C_2_K;
   % thetao75_K and thetao75_C
   ESM_forcing.thetao75_K             = squeeze(ESM_forcing.thetao75(lat,lon,:));
   ESM_forcing.thetao75_C             = squeeze(ESM_forcing.thetao75(lat,lon,:)) - boats.parameters.C_2_K;
 
 end % if (strcmp(model,'data_monthly'))
 
% Add ESM_forcing structure to boats
 boats.ESM_forcing = ESM_forcing;
 
%----------------------------------------------------------------------------------------
% Add other forcing structures to boats
%----------------------------------------------------------------------------------------
 
% Price (price_forcing)
 load price_forcing.mat
 boats.price_forcing = price_forcing;
% Catchability (catchability_forcing)
 load catchability_forcing.mat
 boats.catchability_forcing = catchability_forcing;

%-----------------------------------------------------------------------------------------
% Set initial dfish and effort (if needed)
%-----------------------------------------------------------------------------------------

switch boats.parameters.lname_rest

  %--------------------------------------------------------------------------------------
  % no initial condition given
  case 'none'
   
    disp(['loading restart.... name: ' boats.parameters.lname_rest]); 
    % set dfish to 1e-6 in each mass class
    %boats.initial.dfish = 1e-6 + zeros(nfish,nfmass);
    boats.initial.dfish = zeros(nfish,nfmass);
    % Economic harvesting
    % Set effort to zero in each group
    if idoecon==1
      boats.initial.effort = zeros(1,nfish);
    end

  %--------------------------------------------------------------------------------------
  % Initial dfish is set by analytical primary-production regime
  case 'PP'

    disp(['loading restart.... name: ' boats.parameters.lname_rest]);
    % set parameters
    tro_sca        = log10(boats.parameters.te)/log10(boats.parameters.ppmr);
    b_allo         = boats.parameters.b_allo;
    h_allo         = boats.parameters.h_allo;
    kappa_eppley   = boats.parameters.kappa_eppley;
    E_activation_A = boats.parameters.E_activation_A;
    k_Boltzmann    = boats.parameters.k_Boltzmann;
    zeta1          = boats.parameters.zeta1;
    Prod_star      = boats.parameters.Prod_star;
    mc_phy_l       = boats.parameters.mc_phy_l;
    mc_phy_s       = boats.parameters.mc_phy_s;

    %------------------------------------------------------------------------------------
    % Set npp and temperature maps for dfish initial state   
    % Use annual averages
    %------------------------------------------------------------------------------------
    npp          = squeeze(nanmean(ESM_forcing.mclim_intpp));        % mmolC m-2 s-1
    npp_ed       = squeeze(nanmean(ESM_forcing.mclim_intpp_ed));     % mmolC m-3 d-1
    temp_phyto   = squeeze(nanmean(ESM_forcing.mclim_thetao75_C));   % degC
    temp_fish    = squeeze(nanmean(ESM_forcing.mclim_thetao75_K));   % degK

    %------------------------------------------------------------------------------------
    % Calculate quantities required for dfish 

    s_over_p   = ( -1.0 + ( 1.0 + 4.0 .* npp_ed ./ (exp(kappa_eppley.*temp_phyto) .* ...
      Prod_star) ).^0.5) .* 0.5;
    frac_lg_du = s_over_p ./ (1.0 + s_over_p); % large fraction of PP as in Dunne et al. (2005)
    mphyto     = (mc_phy_l.^frac_lg_du) .* (mc_phy_s.^(1.0 - frac_lg_du));
  
    temp_dep_A = exp( (-E_activation_A/k_Boltzmann) .* (1./temp_fish - 1./boats.parameters.temp_ref_A));
    A          = (boats.parameters.A00/boats.parameters.spery)*temp_dep_A; % growth rate of Andersen and Beyer (2013, p. 18)
    mortality0 = (exp(zeta1)/3)*A;
    
    %-------------------------------------------------------------------------------------
    % calculate initial dfish
    dfish(1,:,:) = (1/nfish) * (1 - tro_sca) .* npp ./ ...
      ( mortality0 .* mphyto^(tro_sca) .* boats.parameters.minf_2d.^(h_allo + b_allo - 1)) ...
      .* boats.parameters.fmass_2d.^(tro_sca + h_allo - 1);
    % make non existent cells NaNs
    dfish(1,boats.parameters.mask_notexist_2d) = NaN;
    boats.initial.dfish = dfish(1,:,:);

    % Economic harvesting
    % Set effort to zero in each group
    if idoecon==1
      boats.initial.effort = zeros(1,nfish);
    end

  %--------------------------------------------------------------------------------------
  % Initial dfish is set by another condition
  %--------------------------------------------------------------------------------------
   
  otherwise
 
    %--------------------------------------------------------------------------------------
    % Initial dfish if specified restart file lname_rest IS NOT IN in the working directory
    if ~exist(['restart_' boats.parameters.lname_rest '.mat'],'file')
     
      % set dfish to 1e-6 in each mass class
      boats.initial.dfish = 1e-6 + zeros(nfish,nfmass);
  
      % Economic harvesting
      % Set effort to zero in each group
      if idoecon==1
        boats.initial.effort = zeros(1,nfish);
      end

     %--------------------------------------------------------------------------------------
     % Initial dfish if specified restart file lname_rest IS IN the working directory
     else
     
       disp(['loading restart.... name: ' boats.parameters.lname_rest]);
       % Load restart file
       tmp = load(['restart_' boats.parameters.lname_rest '.mat']);
       restart = tmp.restart;
       % Use dfish from restart
       boats.initial.dfish  = restart.dfish;
       
       % Economic harvesting
       if idoecon==1
         
         % Use effort there is a field named effort in the restart file
         if isfield(restart,'effort')
           boats.initial.effort  = restart.effort;
         else
         % Set effort to zero in each group
           boats.initial.effort = zeros(1,nfish);
         end % if isfield(restart,'effort')
             
       end % if idoecon==1
   end %if ~exist(['restart_' boats.parameters.lname_rest '.mat'],'file')
end % switch boats.parameters.lname_rest

%----------------------------------------------------------------------------------------
% END OF SCRIPT
