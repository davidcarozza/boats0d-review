function boats = boats0d_save_restart(boats)

%-----------------------------------------------------------------------------------------
% boats0d_save_restart(boats)
% Define the standard input parameters
%-----------------------------------------------------------------------------------------

%-----------------------------------------------------------------------------------------
% Save dfish to restart
 filename = boats.parameters.sname_rest; 
 restart.dfish 	= squeeze(boats.dfish(end,:,:));

%-----------------------------------------------------------------------------------------
% Economic harvesting
% Save effort to restart
 if boats.parameters.idoecon == 1
    restart.effort = squeeze(boats.effort(end,:,:));
 end

%-----------------------------------------------------------------------------------------
% Display message and save restart 
 disp(['saving restart.... name: ' boats.parameters.sname_rest]);
 save(['restart_' filename],'restart');

%-----------------------------------------------------------------------------------------
% Save restart to boats
 boats.restart = restart;

%-----------------------------------------------------------------------------------------
% END OF SCRIPT