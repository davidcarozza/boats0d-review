function newboats = boats0d_time_average(boats,varargin)

%-----------------------------------------------------------------------------------------
% boats0d_time_average(boats,varargin)
% Loop through the variables and remove the time dimension by taking an average over
% timesteps for all variables with time as a dimension
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Define standard input
 A.timesteps = nan;	% time step indeces for averaging - use NaN to take the last timestep
 A.timesteps_frac = 0.4;

%-----------------------------------------------------------------------------------------
% Parse required variables, substituting defaults where necessary
 A = parse_pv_pairs(A, varargin);
%-----------------------------------------------------------------------------------------

 newboats = boats;

 ntime = boats.ntime;

%-----------------------------------------------------------------------------------------
% Select the indices for averaging
 if isnan(A.timesteps)
    istart = ntime-ceil(A.timesteps_frac*ntime);
    iend   = ntime;
 elseif ndims(A.timesteps)==1
    istart = A.timesteps(1);
    iend   = ntime;
 else
    istart = A.timesteps(1);
    iend   = A.timesteps(2);
 end
 
 allvar = fieldnames(boats);
 nvar = length(allvar);

 for indv=1:nvar
    newvar = boats.(allvar{indv});
    if isnumeric(newvar)
       if any(size(newvar)==ntime)
          tdim = find(size(newvar)==ntime);
          meanvar = generalmean(newvar,tdim,istart,iend);
          newboats.(allvar{indv}) = meanvar;
       end
    end
 end

 newboats.timeaverage = 1;

%-----------------------------------------------------------------------------------------
% END OF SCRIPT