function boats0d_plot_base(boats)
%-----------------------------------------------------------------------------------------
% boats0d_plot_base(boats)
% For the time dependent case plot fish time series and dfish spectra
% For the time dependent case plot harvest time series and dharvest spectra 
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% For the time dependent case plot fish time series and dfish spectra
 if size(boats.dfish,1)==boats.ntime

    fig1 = figure;

    subplot(2,1,1)
    plot(boats.year,squeeze(nansum(boats.fish,3)),'-','linewidth',3)
    title('fish biomass by group','fontsize',12);
    xlabel('time (years)','fontsize',12);
    ylabel('biomass (gwB m^{-2})','fontsize',12);
    xlim([boats.year(1) boats.year(end)]);

    subplot(2,1,2)
    plot(log10(boats.fmass),squeeze(log10(boats.dfish(end,:,:))),'.-','markersize',10)
    hold on
    plot(log10(boats.fmass),squeeze(log10(  nansum(boats.dfish(end,:,:),2))),'.-k','markersize',10)
	hold off
    title('equilibrium fish spectra by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('spectrum (log10 gwB m^{-2} g^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);

    en_input_vb = squeeze(boats.en_input_vb(end-1,:,:));
    en_input_vb(boats.mask_notexist) = NaN;
    en_input_P = squeeze(boats.en_input_P(end-1,:,:));
    en_input_P(boats.mask_notexist) = NaN;
    en_input = squeeze(boats.en_input(end-1,:,:));
    gamma = squeeze(boats.gamma(end-1,:,:));
   
    % want the min and max of the y axes in this figure to all be the same
    vec_all = [en_input_vb(:)' en_input_P(:)' en_input(:)' gamma(:)'];
    ylim_min = floor(log10(min(vec_all)));
    ylim_max = ceil(log10(max(vec_all)));

    figure
    subplot(2,2,1)
    
    plot(log10(boats.fmass),log10(en_input_vb),'.-','markersize',10)
    title('input energy VB by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('rate (log10 g s^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]); ylim([ylim_min ylim_max]);
    
    subplot(2,2,2)
    
    plot(log10(boats.fmass),log10(en_input_P),'.-','markersize',10)
    title('input energy P by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('rate (log10 g s^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]); ylim([ylim_min ylim_max]);

    subplot(2,2,3)
    
    plot(log10(boats.fmass),squeeze(log10(boats.en_input(end-1,:,:))),'.-','markersize',10)
    title('input energy by group','fontsize',12);
    xlabel('size (log10 mass)','fontsize',12);
    ylabel('rate (log10 g s^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]); ylim([ylim_min ylim_max]);

    subplot(2,2,4)
    
    plot(log10(boats.fmass),squeeze(log10(boats.gamma(end-1,:,:))),'.-','markersize',10)
    title('growth rate by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('rate (log10 g s^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]); ylim([ylim_min ylim_max]);

    figure
    
    flux_in = squeeze(boats.flux_in(end-1,:,:));
    flux_out = squeeze(boats.flux_out(end-1,:,:));
    flux_fish_growth = squeeze(boats.flux_out(end-1,:,:));
    mortality = boats.mortality(end-1,:,:);
   
    % want the min and max of the y axes in this figure to all be the same
    vec_all = [flux_in(:)' flux_out(:)' flux_fish_growth(:)' mortality(:)'];
    ylim_min = floor(log10(min(vec_all)));
    ylim_max = ceil(log10(max(vec_all)));
    
    subplot(2,3,1)
    
    plot(log10(boats.fmass),squeeze(log10(boats.flux_in(end-1,:,:))),'.-','markersize',10)
    title('flux in by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('flux in (log10 gwB m^{-2}s^{-1}g^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);

    subplot(2,3,2)
    
    plot(log10(boats.fmass),squeeze(log10(boats.flux_out(end-1,:,:))),'.-','markersize',10)
    title('flux out by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('flux out (log10 gwB m^{-2}s^{-1}g^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);

    subplot(2,3,3)
    
    plot(log10(boats.fmass),squeeze(log10(boats.flux_in(end-1,:,:)) - log10(boats.flux_out(end-1,:,:))),'.-','markersize',10)
    title('flux in - flux out by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('flux in - flux out (log10 gwB m^{-2}s^{-1}g^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);

    subplot(2,3,4)
    plot(log10(boats.fmass),squeeze(log10(boats.flux_fish_growth(end-1,:,:))),'.-','markersize',10)
    title('growth by size (term 2 of equation 1) by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('growth by size (log10 gwB m^{-2}s^{-1}g^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]); ylim([ylim_min ylim_max]);

    subplot(2,3,5)
    
    plot(log10(boats.fmass),squeeze(log10(boats.mortality(end-1,:,:))),'.-','markersize',10)
    title('natural mortality by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('natural mortality (log10 gwB m^{-2}s^{-1}g^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]); ylim([ylim_min ylim_max]);

    figure
    
    subplot(2,2,1)
    mortality = boats.mortality0 .* boats.parameters.minf_2d.^(boats.parameters.h_allo + boats.parameters.b_allo - 1) .* boats.parameters.fmass_2d.^(-boats.parameters.h_allo);
    mortality(boats.mask_notexist) = NaN;

    plot(log10(boats.fmass),log10(squeeze(mortality)),'.-','markersize',10)
    title('natural mortality rate by group','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('mortality rate (log10 s^{-1})','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);

    % plot fraction of total contribution to the flux_in_rep from each mass class for each group

    flux_in_rep = nansum( boats.rep_alloc_frac .* squeeze(boats.en_input(end-1,:,:).*boats.dfish(end-1,:,:)) .* ...
    boats.parameters.delfm_2d ./ boats.parameters.fmass_2d,2);
    
    flux_in_rep_part = boats.rep_alloc_frac .* squeeze(boats.en_input(end-1,:,:).*boats.dfish(end-1,:,:)) .* ...
    boats.parameters.delfm_2d ./ boats.parameters.fmass_2d;

    flux_in_rep_frac = flux_in_rep_part ./ repmat(flux_in_rep,[1 boats.nfmass]);

    subplot(2,2,2)    
    plot(log10(boats.fmass),flux_in_rep_frac,'.-','markersize',10)
    title('fraction of total energy to egg reproduction','fontsize',12);
    xlabel('size (log10 mass)','fontsize',12);
    ylabel('fraction','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);
    
    rep_alloc_frac = boats.rep_alloc_frac;
    rep_alloc_frac(boats.mask_notexist) = NaN;

    subplot(2,2,3)    
    plot(log10(boats.fmass),rep_alloc_frac,'.-','markersize',10)
    title('fraction of input energy to growth','fontsize',12);
    xlabel('size (log10 g)','fontsize',12);
    ylabel('fraction','fontsize',12);
    xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);  
   
 end

%-----------------------------------------------------------------------------------------
% For the time dependent case plot harvest time series and dharvest spectra
 
 % Economic harvesting
 if boats.parameters.idoecon==1
 
   if size(boats.dfish,1)==boats.ntime
     fig2 = figure;

     subplot(3,1,1)
     plot(boats.year(1:end-1),squeeze(nansum(boats.harvest(1:end-1,:,:),3)*boats.parameters.spery),'-','linewidth',3)
     title('harvest by group','fontsize',12);
     xlabel('time (years)','fontsize',12);
     ylabel('harvest (gwB m^{-2} y^{-1}','fontsize',12);
     xlim([boats.year(1) boats.year(end)]);

     subplot(3,1,2)
     plot(log10(boats.fmass),squeeze(log10(boats.dharvest(end-1,:,:))),'.-','markersize',10)
     hold on;
     plot(log10(boats.fmass),squeeze(log10(nansum(boats.dharvest(end-1,:,:),2))),'.-k','markersize',10)
     hold off;
     title('harvest spectra by group','fontsize',12);
     xlabel('size (log10 g)','fontsize',12);
     ylabel('harvest distribution (log10 gwB m^{-2}s^{-1}g^{-1})','fontsize',12);
     xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);
     ylim([-20 0]);

     subplot(3,1,3)
     plot(log10(boats.fmass),squeeze(boats.selectivity'),'.-','markersize',10)
     title('selectivity by group','fontsize',12);
     xlabel('size (log10 g)','fontsize',12);
     ylabel('selectivity','fontsize',12);
     xlim([log10(boats.fmass(1)) log10(boats.fmass(end))]);
     ylim([0 1.05]);
    
   end
 end 

%-----------------------------------------------------------------------------------------
% END OF FUNCTION