function boats = boats0d_add_diagnostics(boats,varargin)

%-----------------------------------------------------------------------------------------
% boats0d_add_diagnostics.m
% Adds simple diagnostics to a run
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Define standard input
 A.placeholder = nan;

%-----------------------------------------------------------------------------------------
% Parse required variables, substituting defaults where necessary
 A = parse_pv_pairs(A, varargin);
%-----------------------------------------------------------------------------------------
 
 disp(['Adding diagnostics']);

%-----------------------------------------------------------------------------------------
% Works with the time-averaged boats output
 if ~isfield(boats,'timeaverage')
    tboats = boats0d_time_average(boats);
 else
    if boats.timeaverage==0
       tboats = boats0d_time_average(boats);
    else
       tboats = boats;
    end
 end

%-----------------------------------------------------------------------------------------
% Total diagnostics
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------  
% Fish mass
 
 diagnostics.total_fish = nansum(tboats.fish(:));

%-----------------------------------------------------------------------------------------  
% Intercept mass class 1 (log10)

 tdfish = squeeze(nansum(tboats.dfish,1)); % add through all groups
 diagnostics.total_intercept_mc1 = log10(tdfish(1));

%-----------------------------------------------------------------------------------------  
% Total non-reproducing slope
 fmass = tboats.fmass;
 tdfish = squeeze(nansum(tboats.dfish,1)); % add through all groups
 % Ignore mass classes above malpha group 2
 ibad = (fmass > tboats.malpha(2)) | isnan(tdfish);
 fmass(ibad) = [];
 tdfish(ibad) = [];
 [pp ss] = polyfit(log10(fmass),log10(tdfish),1);
 diagnostics.total_nr_slope_loglog = pp(1);
  
%----------------------------------------------------------------------------------------
% Total reproducing slope
  fmass  = tboats.fmass;
  tdfish = squeeze(nansum(tboats.dfish,1)); % add through all groups
  % Ignore mass classes below malpha group 2
  ibad = (fmass < tboats.malpha(2)) | isnan(tdfish);
  fmass(ibad)  = [];
  tdfish(ibad) = [];
  [pp ss]      = polyfit(log10(fmass),log10(tdfish),1);
  diagnostics.total_r_slope_loglog = pp(1);
 
%----------------------------------------------------------------------------------------
% Group diagnostics - broken down (not strictly necessary)
 if (1)

   %--------------------------------------------------------------------------------------
   % Fish mass
   for indg=1:tboats.nfish
     group_fish = squeeze(nansum(tboats.fish(indg,:)));
     diagnostics.(['group' num2str(indg) '_fish']) = group_fish;
   end

   %--------------------------------------------------------------------------------------
   % intercept mass class 1 (log10)
   for indg=1:tboats.nfish
     tdfish = squeeze(tboats.dfish(indg,:));
     diagnostics.(['group' num2str(indg) '_intercept_mc1']) = log10(tdfish(1));
   end

   %--------------------------------------------------------------------------------------
   % Non-reproducing slope
   for indg=1:tboats.nfish
     fmass = tboats.fmass;
     tdfish = squeeze(tboats.dfish(indg,:));
     % Ignore mass classes above malpha
     ibad = (fmass > tboats.malpha(indg)) | isnan(tdfish);
     fmass(ibad) = [];
     tdfish(ibad) = [];
     [pp ss] = polyfit(log10(fmass),log10(tdfish),1);
     diagnostics.(['group' num2str(indg) '_nr_slope_loglog']) = pp(1);
   end

   %--------------------------------------------------------------------------------------    
   % Reproducing slope (log10-log10)
   for indg=1:tboats.nfish
     fmass = tboats.fmass;
     tdfish = squeeze(tboats.dfish(indg,:));
     % Ignore mass classes below malpha
     ibad = (fmass < tboats.malpha(indg)) | isnan(tdfish) ;
     fmass(ibad) = [];
     tdfish(ibad) = [];
     [pp ss] = polyfit(log10(fmass),log10(tdfish),1);
     diagnostics.(['group' num2str(indg) '_r_slope_loglog']) = pp(1);
   end
         
   %-------------------------------------------------------------------------------------- 
   % rec_regime
   for indg=1:tboats.nfish
     rec_regime = squeeze(tboats.rec_regime(indg));
     % Gets rid of NaNs
     ibad = isnan(rec_regime);
     rec_regime(ibad) = [];
     diagnostics.(['group' num2str(indg) '_recregime']) = rec_regime;
   end

   %--------------------------------------------------------------------------------------   
   % flux_in_rep
   for indg=1:tboats.nfish
     flux_in_rep = squeeze(tboats.flux_in_rep(indg));
     % Gets rid of NaNs
     ibad = isnan(flux_in_rep);
     flux_in_rep(ibad) = [];
     diagnostics.(['group' num2str(indg) '_fluxinrep']) = flux_in_rep;
   end

   %--------------------------------------------------------------------------------------
   % flux_in_P
   for indg=1:tboats.nfish
     flux_in_P = squeeze(tboats.flux_in_P(indg));
     % Gets rid of NaNs
     ibad = isnan(flux_in_P);
     flux_in_P(ibad) = [];
     diagnostics.(['group' num2str(indg) '_fluxinP']) = flux_in_P;
   end
    
   %-------------------------------------------------------------------------------------- 
   % flux_in
   for indg=1:tboats.nfish
     flux_in = squeeze(tboats.flux_in(indg,1));
     % Gets rid of NaNs
     ibad = isnan(flux_in);
     flux_in(ibad) = [];
     diagnostics.(['group' num2str(indg) '_fluxin']) = flux_in;
   end
  
   %--------------------------------------------------------------------------------------
   % total flux_in_rep, flux_in_P, and flux_in 
   diagnostics.total_flux_in_rep = nansum(tboats.flux_in_rep(:));
   diagnostics.total_flux_in_P   = nansum(tboats.flux_in_P(:));
   diagnostics.total_flux_in     = nansum(tboats.flux_in(:,1));

   %-------------------------------------------------------------------------------------- 
   % ena_regime (energy available)
   for indg=1:tboats.nfish
     ena_regime = tboats.ena_regime(indg,:);
     % Gets rid of NaNs
     ibad = isnan(ena_regime);
     ena_regime(ibad) = [];
     diagnostics.(['group' num2str(indg) '_enaregime']) = ena_regime;
   end
  
   %--------------------------------------------------------------------------------------
   % ratio_Myers
   for indg=1:tboats.nfish
     ratio_Myers = tboats.ratio_Myers(indg,:);
     % Gets rid of NaNs   
     ibad = (isnan(ratio_Myers));
     ratio_Myers(ibad) = [];
     diagnostics.(['group' num2str(indg) '_ratioMyers']) = ratio_Myers;
   end

   %--------------------------------------------------------------------------------------
   % flux_in_spawners
   for indg=1:tboats.nfish
     flux_in_spawners = tboats.flux_in_spawners(indg,:);
     % Gets rid of NaNs   
     ibad = (isnan(flux_in_spawners));
     flux_in_spawners(ibad) = [];
     diagnostics.(['group' num2str(indg) '_flux_in_spawners']) = flux_in_spawners;
   end
   
   %--------------------------------------------------------------------------------------
   % number_spawners
   for indg=1:tboats.nfish
     number_spawners = tboats.number_spawners(indg,:);
     % Gets rid of NaNs   
     ibad = (isnan(number_spawners));
     number_spawners(ibad) = [];
     diagnostics.(['group' num2str(indg) '_number_spawners']) = number_spawners;
   end
  
   %--------------------------------------------------------------------------------------
   % Economic harvesting
   if (boats.parameters.idoecon==1)

     %------------------------------------------------------------------------------------
     % total harvest  
     diagnostics.total_harvest = nansum(tboats.harvest(:))*boats.parameters.spery;
     
     %------------------------------------------------------------------------------------ 
     % group harvest
     for indg=1:tboats.nfish
       group_harvest = squeeze(nansum(tboats.harvest(indg,:)))*boats.parameters.spery;
       diagnostics.(['group' num2str(indg) '_harvest']) = group_harvest;
     end
     
    end
    
end % if (1)

%----------------------------------------------------------------------------------------
% Save diagnostics to boats
 boats.diagnostics = diagnostics;

%----------------------------------------------------------------------------------------
% END OF SCRIPT
