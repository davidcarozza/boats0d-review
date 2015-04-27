function boats = boats0d_integrate(boats)

%-----------------------------------------------------------------------------------------
% boats2d_integrate_vec.m
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Authors
%-----------------------------------------------------------------------------------------

% David A. Carozza* (david.carozza@gmail.com; corresponding author)
% Daniele Bianchi   (danbian@uw.edu)
% Eric D. Galbraith (eric.galbraith@mcgill.ca)

%-----------------------------------------------------------------------------------------
% Introduction
%-----------------------------------------------------------------------------------------

% Bioeconomic Open-Access Trophic Size-Spectrum (BOATS) model based on the
% McKendrick-von Foerster model for a mass-spectrum of fish biomass
% with an open access economic framework to calculate effort and harvest.
% Forced with monthly primary production and temperature data (from model or observations)

% Primary production in units of mmolC m-2 s-1
% Core unit of biomass is mmolC for fish and harvest plots in model
% Convert these to grams of wet biomass (g wetB) using mmolC_2_wetB
% Time in seconds
% Time step is 30 days (1 month)

% 3 fish groups defined by different asymptotic masses
% 1 effort group and 1 selectivity function per fish group
% 0-dimensional version of the code

%-----------------------------------------------------------------------------------------
% Other comments
%-----------------------------------------------------------------------------------------

% We set a lower limit on effort by adding a small constant (epsln = 1e-15)
% to each application of effort when calculating harvest, cost, and effort change.
% This guarantees that as effort approaches zero, the effort change equation approaches 
% the analytical simplification. This prevents dividing by zero.
% We cannot use the analytical simplification of the change in effort directly because we
% limit harvest and we need harvest directly to calculate revenue.

%-----------------------------------------------------------------------------------------
% MAIN CODE
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Load monthly observational data and grid structure
%-----------------------------------------------------------------------------------------

 data_monthly = boats.parameters.data_monthly;

%-----------------------------------------------------------------------------------------
% Define general parameters and variables
%-----------------------------------------------------------------------------------------

%-------------------------------------------------------------------------------------
% Site latitude and longitude

 lat = boats.parameters.lat;
 lon = boats.parameters.lon;

%-------------------------------------------------------------------------------------
% Conversion factors
 sperd        = boats.parameters.sperd;
 spery 		  = boats.parameters.spery;
% 12 is molC -> gC; /1000 is mg to g; *10 is gC to wetB; primary production data is in m molC
 gC_2_wetB    = boats.parameters.gC_2_wetB; % gC to grams wet biomass
 mmolC_2_wetB = (12*gC_2_wetB)/1000; % [wetB / mmolC]
 epsln 		  = boats.parameters.epsln;

%-------------------------------------------------------------------------------------
% Time integration parameters
 dtts       = boats.parameters.dtts * sperd;		% days * seconds/day = seconds
 run_length = boats.parameters.run_length;			% years
 time       = 1:dtts:run_length*spery;				% seconds
 time       = time(:);
 ntime      = length(time);
 year       = time / spery;

%-------------------------------------------------------------------------------------
% Forcing
 iforcing_ESM          = boats.parameters.iforcing_ESM;
 iforcing_price        = boats.parameters.iforcing_price;
 iforcing_catchability = boats.parameters.iforcing_catchability;

 start_year            = boats.parameters.start_year;
 start_year_price      = boats.parameters.start_year_price;

%-----------------------------------------------------------------------------------------
% Model choice indices
 idoecon        = boats.parameters.idoecon;       % economics index (1 = do economics)
 
%----------------------------------------------------------------------------------
% Mass class structure (bin boundaries in logarithmic space)
 nfmass   = boats.parameters.nfmass; 	% number of fish mass classes
 fmass_0  = boats.parameters.fmass_0;	% initial mass class (g)
 fmass_e  = boats.parameters.fmass_e;	% final mass class (g)
 ifmbound = boats.parameters.ifmbound;	% size class number bounds: 1,2,3, ...
 fmbound  = boats.parameters.fmbound;   % mass class bounds (center is geometric mean of all mass classes)
 delfm    = boats.parameters.delfm;     % width of mass classes
 fmass    = boats.parameters.fmass;
 
% mass and mass class with for boundary conditions only
% dummy previous mass class using same scaling
 fmass_bc = fmass(1) * fmass(1)/fmass(2);
 delfm_bc = fmbound(1) - fmbound(1) * fmbound(1)/fmbound(2);

%-------------------------------------------------------------------------------------
% Asymptotic masses for fish groups
 minf      = boats.parameters.minf;   % asymptotic mass (g)
 eta_alpha = boats.parameters.eta_alpha; % maturation to asymptotic mass ratio
 malpha    = boats.parameters.malpha; % mass of first reproduction (g)
 nfish     = boats.parameters.nfish; % number of groups (number of effort groups)
 
%-------------------------------------------------------------------------------------
% Save out frequently used arrays and constants
% repmats of minf, malpha, fmass, delfm, delfm_2end
 fmass_2d          = repmat(fmass,[nfish 1]);
 minf_2d           = repmat(minf',[1 nfmass]);
 malpha_2d         = repmat(malpha',[1 nfmass]); 
 delfm_2d          = repmat(delfm,[nfish 1]);
 delfm_2end_2d     = repmat(delfm(2:end),[nfish 1]);

%-------------------------------------------------------------------------------------
% Define mass classes of adult in each group (used to calculate ratio_Myers)
 transition = [12 30 43];

%-------------------------------------------------------------------------------------
% Temperature dependence
 E_activation_A	= boats.parameters.E_activation_A;    % eV  (Savage et al., 2004)
 E_activation_m = boats.parameters.E_activation_m;  % eV  (Save et al., 2004)
 k_Boltzmann 	= boats.parameters.k_Boltzmann;     % eV K-1
 temp_ref_A     = boats.parameters.temp_ref_A;      % K (Andersen and Beyer, 2013, p. 18)
 
%-------------------------------------------------------------------------------------
% Primary production
 kappa_eppley   = boats.parameters.kappa_eppley; % degC-1
 Prod_star      = boats.parameters.Prod_star;    % mmolC m-3 d-1
 C_2_K          = boats.parameters.C_2_K;        % celcius to kelvin
 mc_phy_l       = boats.parameters.mc_phy_l;     % mass of typical large phytoplankton
 mc_phy_s       = boats.parameters.mc_phy_s;     % mass of typical small phytoplankton
 
%-------------------------------------------------------------------------------------
% Ecology
 te 	     = boats.parameters.te;          % trophic efficiency
 ppmr 	     = boats.parameters.ppmr;        % predator to prey mass ratio
 tro_sca     = log10(te)/log10(ppmr);        % trophic scaling
 b_allo      = boats.parameters.b_allo;		 % allometric scaling
 mortality00 = boats.parameters.mortality00; % constant mortality scaling
 h_allo      = boats.parameters.h_allo;      % mass scaling of mortality
 eff_a       = boats.parameters.eff_a;       % efficiency of activity (Andersen and Beyer, 2013, p. 4)
 A00         = boats.parameters.A00;         % allometric growth rate (Andersen and Beyer, 2013, p. 4)
 A0          = boats.parameters.A00/boats.parameters.spery; % growth rate per second (Andersen and Beyer, 2013, p. 4)
 
%-------------------------------------------------------------------------------------
% Reproduction
 m_egg           = boats.parameters.m_egg;            % egg mass (g)
 frac_fem        = boats.parameters.frac_fem;         % fraction of individuals that allocate energy to reproduction (females)
 egg_surv        = boats.parameters.egg_surv;         % egg survival
 rep_slope       = boats.parameters.rep_slope;        % slope parameter of allocation to reproduction function
 rep_pos         = boats.parameters.rep_pos;          % position parameter of allocation to reproduction function as fraction of malpha

% calculate mass scaling of energy allocation to reproduction  
 rep_scale                 = nan(nfish,nfmass);
 rep_scale_Myers           = nan(nfish,nfmass); 
 for indf = 1:nfish
   rep_scale(indf,:)       = sigmoid_And_mass(fmass,rep_pos*malpha(indf),rep_slope);
   rep_scale_Myers(indf,:) = sigmoid_And_mass(fmass,malpha(indf),10);
 end
 
 rep_alloc_frac = rep_scale .* (1 - eff_a) ./ ( (fmass_2d./minf_2d).^(b_allo-1) - eff_a);
 
%-----------------------------------------------------------------------------------------
% Partition of primary production between groups
 part_PP_b = 1/nfish; % partition of primary production at boundary condition (recruitment)
 part_PP_g = 1/nfish; % partition of primary production of growth

%-----------------------------------------------------------------------------------------
% Masks
 mask_notexist = (repmat(fmass,[nfish 1]) > repmat(minf,[nfmass 1])'); % NaN where mass is > asymptotic mass
 
%-----------------------------------------------------------------------------------------
% Forcing data
 ESM_forcing          = boats.ESM_forcing;
 price_forcing        = boats.price_forcing;
 catchability_forcing = boats.catchability_forcing; 
 
%-----------------------------------------------------------------------------------------
% Initialize biological arrays
%-----------------------------------------------------------------------------------------

 dfish 			  = nan(ntime,nfish,nfmass);
 en_input_P 	  = nan(ntime,nfish,nfmass);
 gamma_growth_vb  = nan(ntime,nfish,nfmass);
 gamma_loss_vb    = nan(ntime,nfish,nfmass);
 en_input_vb      = nan(ntime,nfish,nfmass);
 en_input         = nan(ntime,nfish,nfmass);
 gamma	 		  = nan(ntime,nfish,nfmass);
 
 flux_in	 	  = nan(ntime,nfish,nfmass);
 flux_out	 	  = nan(ntime,nfish,nfmass);
 flux_fish_growth = nan(ntime,nfish,nfmass);
 flux_in_rep      = nan(ntime,nfish);
 flux_in_P        = nan(ntime,nfish);
 ena_regime       = nan(ntime,nfish,nfmass);
 rec_regime       = nan(ntime,nfish);
 mortality        = nan(ntime,nfish,nfmass);

 number_spawners  = nan(ntime,nfish);
 flux_in_spawners = nan(ntime,nfish);
 ratio_Myers      = nan(ntime,nfish);
   
%-----------------------------------------------------------------------------------------
% Biomass initial condition
%-----------------------------------------------------------------------------------------

 dfish(1,:,:) = boats.initial.dfish;
 
%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------
% Economics
%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------

if idoecon==1

%-----------------------------------------------------------------------------------------
% Define general parameters and variables
%-----------------------------------------------------------------------------------------
 
 cost_effort       = boats.parameters.cost_effort_0 * ones(1,nfish); % cost per unit effort
 k_e               = boats.parameters.k_e;                           % fleet dynamics parameter
 sel_pos_1         = boats.parameters.sel_pos_1;				     % selectivity position shift 1
 sel_pos_2         = boats.parameters.sel_pos_2;				     % selectivity position shift 2
 sel_pos_3         = boats.parameters.sel_pos_3;				     % selectivity position shift 3
 sel_slope_1       = boats.parameters.sel_slope_1;			         % selectivity slope 1
 sel_slope_2       = boats.parameters.sel_slope_2;			         % selectivity slope 2
 sel_slope_3 	   = boats.parameters.sel_slope_3;			         % selectivity slope 3

%-----------------------------------------------------------------------------------------
% Harvesting selectivity
  
 selectivity = nan(nfish,nfmass);
   
 selectivity(1,:) = sigmoid_And_length(fmass,sel_pos_1*malpha(1),sel_slope_1);
 selectivity(2,:) = sigmoid_And_length(fmass,sel_pos_2*malpha(2),sel_slope_2);
 selectivity(3,:) = sigmoid_And_length(fmass,sel_pos_3*malpha(3),sel_slope_3);

%-----------------------------------------------------------------------------------------
% Initialize economics arrays
 dharvest 	   = nan(ntime,nfish,nfmass);
 dfish_temp    = nan(ntime,nfish,nfmass);           
 effort 	   = nan(ntime,nfish);
 effort_change = nan(ntime,nfish);
 cost	 	   = nan(ntime,nfish);
 revenue	   = nan(ntime,nfish);
   
%-----------------------------------------------------------------------------------------
% Effort initial condition
%-----------------------------------------------------------------------------------------
 
 effort(1,:) =  boats.initial.effort;
   
end % ideoecon

%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------
% MAIN LOOP
%-----------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------

for indt = 1:ntime-1

  % local month
  local_month = ceil((mod(time(indt),spery)/spery)*12);

  %% WILL NEED TO RETHINK THIS. AFTER INITIALIZE_2D

  %---------------------------------------------------------------------------------------
  % Physical and ecological forcings (NPP in mmolC m-2 s-1)
  %---------------------------------------------------------------------------------------
  switch (iforcing_ESM)
    case 1 % force with annual climatology
      npp        = squeeze(nanmean(ESM_forcing.mclim_intpp));          % mmolC m-2 s-1
      npp_ed     = squeeze(nanmean(ESM_forcing.mclim_intpp_ed));     % mmolC m-3 d-1
      temp_phyto = squeeze(nanmean(ESM_forcing.mclim_thetao75_C));     % degC
      temp_fish  = squeeze(nanmean(ESM_forcing.mclim_thetao75_K));   % degK
    case 2 % force with monthly climatology
      npp         = squeeze(ESM_forcing.mclim_intpp(local_month));
      npp_ed      = squeeze(ESM_forcing.mclim_intpp_ed(local_month)); 
      temp_phyto  = squeeze(ESM_forcing.mclim_thetao75_C(local_month));
      temp_fish   = squeeze(ESM_forcing.mclim_thetao75_K(local_month));
    case 3 % force with time-varying changing climate model output
      npp         = squeeze(ESM_forcing.intpp(indt+start_year*12));
      npp_ed      = squeeze(ESM_forcing.intpp_ed(indt+start_year*12)); 
      temp_phyto  = squeeze(ESM_forcing.thetao75_C(indt+start_year*12));
      temp_fish   = squeeze(ESM_forcing.thetao75_K(indt+start_year*12));
    case 4 % force with time-varying constant climate model output (for CM2.6)
      npp         = squeeze(ESM_forcing.mclim_intpp(indt+start_year*12));
      npp_ed      = squeeze(ESM_forcing.mclim_intpp_ed(indt+start_year*12)); 
      temp_phyto  = squeeze(ESM_forcing.mclim_thetao75_C(indt+start_year*12));
      temp_fish   = squeeze(ESM_forcing.mclim_thetao75_K(indt+start_year*12));
    case 10 % constant forcing (in mg C m-2 d-1, converted to m molC m-2 s-1)
      npp         = boats.parameters.npp * (1/sperd) / 12; % mg C m-2 d-1 to mmolC m-2 s-1
      npp_ed      = (boats.parameters.npp / 12) /75;       % mg C m-2 d-1 to m molC m-2 d-1 since lambda0 is d-1
      temp_phyto  = boats.parameters.temp75;               % degC
      temp_fish   = boats.parameters.temp75 + C_2_K;       % degK
    end 
  
  %---------------------------------------------------------------------------------------
  % Large fraction of phytoplankton and representative phytoplankton mass
  %--------------------------------------------------------------------------------------- 

  s_over_p = ( -1.0 + ( 1.0 + 4.0 * npp_ed / (exp(kappa_eppley*temp_phyto) * Prod_star))^0.5) * 0.5;
  frac_lg_du = s_over_p / (1.0 + s_over_p); % large fraction of PP as in Dunne et al. (2005)
  mphyto = (mc_phy_l^frac_lg_du) * (mc_phy_s^(1.0 - frac_lg_du));

  %---------------------------------------------------------------------------------------
  % Growth rate
  %---------------------------------------------------------------------------------------
  
  %---------------------------------------------------------------------------------------
  % Based on primary production
  % growth rate = production distribution * mass / biomass distribution
  % multiply by growth partition function (part_PP_g)
  
  en_input_P(indt,:,:) = (npp/mphyto * (repmat(fmass,[nfish 1])./mphyto).^(tro_sca-1)) .* repmat(fmass,[nfish 1]) ./ (reshape(dfish(indt,:,:)+epsln,[nfish nfmass])) .* ...
   part_PP_g;
    
  %---------------------------------------------------------------------------------------
  % Based on allometric scaling (von Bertalanffy)
  % calculate temperature dependencies, then growth rate, the activity loss constant (ka)

  temp_dep_A   = exp( (-E_activation_A/k_Boltzmann) .* (1/temp_fish - 1/temp_ref_A));
  temp_dep_m   = exp( (-E_activation_m/k_Boltzmann) .* (1/temp_fish - 1/temp_ref_A));
  A 		   = A0*temp_dep_A;
  ka 		   = A*eff_a*minf_2d.^(b_allo-1); % s-1
     
  gamma_growth_vb(indt,:,:)  = A * fmass_2d.^b_allo;
  gamma_loss_vb(indt,:,:)    = ka .* fmass_2d;
  
  en_input_vb(indt,:,:)      = squeeze(gamma_growth_vb(indt,:,:) - gamma_loss_vb(indt,:,:));
    
  %---------------------------------------------------------------------------------------
  % input energy (energy available to growth and reproduction)
  % minimum of en_input_P and en_input_vb

  en_input(indt,:,:)           = min(en_input_P(indt,:,:),en_input_vb(indt,:,:));
  en_input(indt,mask_notexist) = NaN;

  %---------------------------------------------------------------------------------------
  % Somatic growth rate (growth rate of body tissues)

  gamma(indt,:,:)              = (1 - rep_alloc_frac).*squeeze(en_input(indt,:,:));
  %gamma(indt,:,:)              = squeeze(en_input(indt,:,:));
    
  %----------------------------------------------------------------------------------
  % flux out of a mass class
  
  flux_out(indt,:,:) = reshape(gamma(indt,:,:),[nfish nfmass]) .* reshape(dfish(indt,:,:),[nfish nfmass]) ./ ...
    delfm_2d;

  %---------------------------------------------------------------------------------------
  % Biomass growth
  % increase in biomass in growth to next mass class
  % arises in conversion from abundance-based MFVF to biomass-based
  
  flux_fish_growth(indt,:,:) = reshape(gamma(indt,:,:),[nfish nfmass]) .* reshape(dfish(indt,:,:),[nfish nfmass]) ./ ...
    fmass_2d;
    
  %----------------------------------------------------------------------------------
  % boundary condition (flux in to first mass class)
  %----------------------------------------------------------------------------------
    
  %---------------------------------------------------------------------------------------    
  % Boundary condition based on primary production
  % multiply by boundary condition partition function (part_PP_b)
  
  flux_in_P(indt,:) = (npp/mphyto * (fmass_bc/mphyto).^(tro_sca-1) * fmass_bc ./ delfm(1)) * part_PP_b;
    
  %---------------------------------------------------------------------------------------
  % Boundary condition based on recruitment (production and survival of eggs)

  flux_in_rep(indt,:) = frac_fem*(fmass_bc/m_egg)*(egg_surv) .* ...
    nansum( rep_alloc_frac .* squeeze(en_input(indt,:,:).*dfish(indt,:,:)) .* ...
    delfm_2d ./ fmass_2d,2) / delfm(1);
    
  %---------------------------------------------------------------------------------------
  % Boundary condition (Beverton-Holt form)

  flux_in(indt,:,1) = flux_in_P(indt,:) .* flux_in_rep(indt,:) ./ (flux_in_P(indt,:) + flux_in_rep(indt,:) + epsln);
        
  %---------------------------------------------------------------------------------------
  % Flux in of other mass classes
  
  flux_in(indt,:,2:end) = reshape(gamma(indt,:,1:end-1),[nfish (nfmass-1)]) .* ...
    reshape(dfish(indt,:,1:end-1),[nfish (nfmass-1)]) ./ ...
    delfm_2end_2d;

  %---------------------------------------------------------------------------------------
  % Calculate Myers ratio (Myers, 2001, ICES, Stock and recruitment: generalizations about
  % maximum reproductive rate, density dependence, and variability using meta-analytic approaches)
  % Found that new spawners per spawner is between 1 and 7
  
  % calculate number flux of new spawners
  for indf = 1:nfish
    flux_in_spawners(indt,indf) = gamma(indt,indf,transition(indf)-1) * (dfish(indt,indf,transition(indf)-1)) / fmass(transition(indf)-1);

  end

  % calculate number of spawners
  % use rep_scale_Myers since it is (essentially) 1 for spawners and 0 for non-spawners
  number_spawners(indt,:) = nansum( rep_scale_Myers .* squeeze(dfish(indt,:,:)) .* (delfm_2d ./ fmass_2d),2);
  % calculate ratio
  % Multiply by spery to have new spawners per year per spawner
  ratio_Myers(indt,:) = (flux_in_spawners(indt,:) + epsln) ./ (number_spawners(indt,:) + epsln) * spery;

  %---------------------------------------------------------------------------------------
  % Input energy (available energy to growth and reproduction regime
  % 1 if en_input_P less - then en_input_P determines input energy

  ena_regime(indt,:,:) = (reshape(en_input_P(indt,:,:),[nfish nfmass]) < ...
    reshape(en_input_vb(indt,:,:),[nfish nfmass]));
  
  %---------------------------------------------------------------------------------------
  % Recruitment regime
  % 1 if flux_in_P less - then flux_in_P determines recruitment
  rec_regime(indt,:) = (flux_in_P(indt,:) < flux_in_rep(indt,:));
    
  %----------------------------------------------------------------------------------
  % Mortality
  %----------------------------------------------------------------------------------
    
  %---------------------------------------------------------------------------------------  
  % Mass-specific mortality rate 
  % Calculate associated growth rate with mortality temperature dependence temp_dep_m
  % calculate mortality rate mortality0
  % mortality00 is the exp(zeta_1) term from the model description
    
  A          = A0*temp_dep_m;
  mortality0 = mortality00*A/3;

  %---------------------------------------------------------------------------------------  
  % Mortality rate 
  % Charnov et al. (2013)
  % minf_4d_p_hplusbm1 is minf_4d.^(h_allo + b_allo - 1)
  
  mortality(indt,:,:) = mortality0 .* minf_2d.^(h_allo + b_allo - 1) .* fmass_2d.^(-h_allo) .* reshape(dfish(indt,:,:),[nfish nfmass]);
    
  %---------------------------------------------------------------------------------------
  % Integrate fish and effort
  %---------------------------------------------------------------------------------------
    
  if (idoecon == 0) % No economic harvesting

    %-------------------------------------------------------------------------------------
    % Integrate dfish
    %-------------------------------------------------------------------------------------
    
    dfish(indt+1,:,:)  = dfish(indt,:,:) + ( flux_in(indt,:,:) - flux_out(indt,:,:) + flux_fish_growth(indt,:,:) - mortality(indt,:,:) ) * dtts;
    mask_dfish_neg = (squeeze(dfish(indt+1,:,:)) < 0);
    dfish(indt+1,mask_dfish_neg)  = 0;
  
  else
    
    %-------------------------------------------------------------------------------------
    % Update fish to calculate harvest
    %-------------------------------------------------------------------------------------
    
    dfish_temp(indt,:,:) = dfish(indt,:,:) + ( flux_in(indt,:,:) - flux_out(indt,:,:) + flux_fish_growth(indt,:,:) - mortality(indt,:,:) ) * dtts;
    mask_dfish_temp_neg = (squeeze(dfish_temp(indt,:,:)) < 0);
    dfish_temp(indt,mask_dfish_temp_neg)  = 0;
        
    %-------------------------------------------------------------------------------------
    % Catchability 
    %-------------------------------------------------------------------------------------
     
    harvest_start = boats.parameters.harvest_start; % year at which harvest begins
    
    if (time(indt) < harvest_start*spery) % catchability zero before start year
    
        qcatch0 = 0;
    
    else % set catchability forcing scenario
    
      switch (iforcing_catchability)
        case 1 % constant catchability
          qcatch0 = boats.parameters.qcatch0 * 1;
        case 2 % increasing by 3% per year
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p03_190y(indt);
        case 3 % like case 2 but increasing price for 2007-2100
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p03_40y_p04_150y(indt);
        case 4 % like case 3 but going from 1910 to 2100
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p03_40y_p05_150y(indt);
        case 5 % like case 3 but going from 1910 to 2100
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p03_40y_p05_250y(indt);
        case 10 
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p04_150y(indt);  
        case 11
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p0_10y_p04_190y(indt);
        case 12
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p0_10y_p05_190y(indt);
        case 13
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p0_50y_p05_190y(indt);
        case 14
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p03_250y(indt);
        case 15
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p04_250y(indt);
        case 16
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_p05_250y(indt);
        case 21
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_93y_k(indt);
        case 22
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_93y_i5(indt);
        case 23
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_93y_i2p5(indt);
        case 24
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_93y_d2p5(indt);        
        case 31
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_s_93y_k(indt);
        case 32
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_s_93y_i2p5(indt);
        case 33
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_s_93y_d2p5(indt);
        case 41
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_43y_k(indt);
        case 42
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_43y_i5(indt);
        case 43
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_43y_i2p5(indt);
        case 44
          qcatch0 = boats.parameters.qcatch0 * catchability_forcing.catch_43y_d2p5(indt);
        otherwise % error
          error(['Error: iforcing_catchability choice ' num2str(iforcing_catchability) ' not implemented']);
      end
    
      % adjust size of vector for dharvest calculation
      qcatch = qcatch0 * ones(1,nfish);

    end % time(indt) < harvest_start*spery

    %-------------------------------------------------------------------------------------
    % dharvest [nfish,nfmass]
    %-------------------------------------------------------------------------------------
    % qcatch * effort * selectivity * dfish_temp    
    % Use repmat and permute so that all 4 arrays are of dimension [nfish,nfmass]
    % Set upper limit (min) so that no more than the fish present can be caught (dfish_temp/dtts)
    %-------------------------------------------------------------------------------------

    dharvest(indt,:,:) = min(squeeze(dfish_temp(indt,:,:))/dtts, repmat(qcatch,[nfmass 1])' .*  ...
      repmat(effort(indt,:)+epsln,[nfmass 1])' .* selectivity .* reshape(dfish_temp(indt,:,:),[nfish nfmass]));
    
    %----------------------------------------------------------------------------------
    % Price forcing
    % input price ($ g-1), multiply by mmolC_2_wetB to get ($ mmolC-1)
    %----------------------------------------------------------------------------------
    
    switch (iforcing_price)
      case 1 % constant price
        price0      = boats.parameters.price_0*mmolC_2_wetB;
      case 2 % follow SAUP for 1950-2006; 1950 for 1850-1949; 2006 for 2007-2100
        price0      = price_forcing.price_SAUP_standard_future_k(indt+start_year_price*12)*mmolC_2_wetB;
      case 3 % like case 2 but increasing price for 2007-2100
        price0      = price_forcing.price_SAUP_standard_future_i(indt+start_year_price*12)*mmolC_2_wetB;
      case 4 % like case 2 but decreasing price for 2007-2100
        price0      = price_forcing.price_SAUP_standard_future_d(indt+start_year_price*12)*mmolC_2_wetB;
      case 5 % idealized prices for 40 years (5 groups of 8 years each)
        price0      = price_forcing.price_ideal_40y_5_groups(indt)*mmolC_2_wetB;
      otherwise % error
        error(['Error: iforcing_price choice ' num2str(iforcing_price) ' not implemented']);
    end
    % adjust size of vector for dharvest calculation
    price       = price0 * ones(nfish,nfmass);
    
    %-------------------------------------------------------------------------------------
    % revenue [nfish]
    %-------------------------------------------------------------------------------------
    % sum over mass (price * dharvest * delfm)
    % Use repmat and permute so that all 4 arrays are of dimension [nfish]
    %-------------------------------------------------------------------------------------    
    
    revenue(indt,:) = nansum(price .* reshape(dharvest(indt,:,:),[nfish nfmass]) .* delfm_2d,2);

    %-------------------------------------------------------------------------------------
    % cost [nfish]
    %-------------------------------------------------------------------------------------
    % cost_effort * effort
    %-------------------------------------------------------------------------------------

    cost(indt,:)    = cost_effort .* (effort(indt,:) + epsln);
    
    %-------------------------------------------------------------------------------------
    % effort_change [nfish]
    %-------------------------------------------------------------------------------------
         
    effort_change(indt,:) = k_e * (revenue(indt,:) - cost(indt,:)) ./ (effort(indt,:) + epsln);
    
    %-------------------------------------------------------------------------------------
    % integrate dfish [nfish,nfmass]
    %-------------------------------------------------------------------------------------

    dfish(indt+1,:,:) = dfish_temp(indt,:,:) - dharvest(indt,:,:) * dtts;
%    mask_dfish_neg  = (squeeze(dfish(indt+1,:,:)) < 0);
%    dfish(indt+1,mask_dfish_neg) = 0;

    mask_group_biomass_small = ( nansum(squeeze(dfish(indt+1,:,:)) .* delfm_2d,2) < 0.1);    
    mask_dfish_small = repmat(mask_group_biomass_small,[1 nfmass]);
    dfish(indt+1,mask_dfish_small) = 0;

    %-------------------------------------------------------------------------------------
    % integrate effort [nfish]
    %-------------------------------------------------------------------------------------
    
    effort(indt+1,:) = effort(indt,:) + effort_change(indt,:) * dtts;
    mask_effort_neg = (squeeze(effort(indt+1,:)) < 0);
    effort(indt+1,mask_effort_neg)  = 0;
      
    end % idoecon
    
 end % for indt = 1:ntime-1

%-----------------------------------------------------------------------------------------
% Post-Processing
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% array to calculate integrals from distributions (for dfish and dharvest)
 delfm_rep(1,1,:) = delfm;
 ddFF = repmat(delfm_rep,[ntime,nfish,1]);

%-----------------------------------------------------------------------------------------
% calculate fish array and convert to wetB
 fish = dfish .* ddFF * mmolC_2_wetB;
 
% Economic harvesting
 if (idoecon==1)
   
   % calculate temporary fish and convert to wetB
   fish_temp = dfish_temp .* ddFF * mmolC_2_wetB;

   % calculate fish harvest and convert to wetB
   harvest = dharvest .* ddFF * mmolC_2_wetB;
 end

%-----------------------------------------------------------------------------------------
% Save arrays to boats structure
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Structural variables

 boats.year = year;
 boats.time = time;
 boats.ntime = ntime;
 boats.nfish = nfish;
 boats.nfmass = nfmass;
 boats.fmass = fmass;
 boats.malpha = malpha;
 boats.delfm = fmass;
 boats.delfm_2d = delfm_2d;
 boats.fmass_2d      = fmass_2d; 
 boats.mask_notexist = mask_notexist;
 boats.rep_scale = rep_scale;
 boats.rep_alloc_frac = rep_alloc_frac;

%-----------------------------------------------------------------------------------------
% Model variables

 boats.en_input_P       = en_input_P;
 boats.temp_dep_A       = temp_dep_A;
 boats.temp_dep_m       = temp_dep_m;
 boats.gamma_growth_vb  = gamma_growth_vb;
 boats.gamma_loss_vb    = gamma_loss_vb;
 boats.en_input_vb      = en_input_vb; 
 boats.en_input         = en_input;
 boats.gamma            = gamma;
 boats.flux_out         = flux_out;
 boats.flux_fish_growth = flux_fish_growth;
 boats.flux_in_P        = flux_in_P;
 boats.flux_in_rep      = flux_in_rep;
 boats.flux_in          = flux_in;
 boats.flux_in_spawners = flux_in_spawners;
 boats.number_spawners  = number_spawners;
 boats.ratio_Myers      = ratio_Myers;
 boats.ena_regime       = ena_regime;
 boats.rec_regime       = rec_regime;
 boats.mortality0       = mortality0;
 boats.mortality        = mortality;
 boats.dfish            = dfish;
 boats.fish             = fish;
 
%-----------------------------------------------------------------------------------------
% Economic harvesting variables 

 if idoecon==1
   boats.dfish_temp       = dfish_temp;
   boats.fish_temp        = fish_temp; 
   boats.selectivity   = selectivity;
   boats.dharvest      = dharvest;
   boats.harvest       = harvest;
   boats.revenue       = revenue;
   boats.cost          = cost;
   boats.effort_change = effort_change;
   boats.effort        = effort;
 end

%-----------------------------------------------------------------------------------------
% END OF SCRIPT
