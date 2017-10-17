function [Ta_out MassFlow] = CoolingTower(P_w,options)
% COOLINGTOWER is a cooling tower 0D modelisation
% COOLINGTOWER(P_w,options) compute the thermodynamics states for a Cooling
% tower based on several inputs (given in OPTION) and based on a given 
% water power to dissipate P_w.
% It returns the main results. 
%
% INPUTS :
% P_W = Heat power output at the condenser [W]
% OPTIONS is a structure containing :
%   -options.Tw_out [°C]: Cooling water temperature at the condenser outlet
%   -options.Tw_in  [°C]: Cooling water temperature at the condenser inlet
%   -options.Tpinch [°C]: Minimum tempearture pinch between Tw_out and the
%                         outlet air temperature.
%   -options.Ta_in  [°C]: Atmospheric air temperature 
%   -options.Phi_atm [-]: Relative humidity of atmospheric air
%   -options.Phi_out [-]: Maximum relative humidity of air at the cooling 
%                         tower outlet.
%
% OUTPUT :
% Ta_out     [°C]: Air temperature at the outlet
% MassFlow [kg/s]: Vector containing the different massflow :
%   -massflow(1) : water massflow at the condenser
%   -massflow(2) : additionnal water massflow = water flow evaporated
%   -massflow(3) : air massflow at the cooling tower   


% Example of how to handle with options structure

if isfield(options,'Tpinch')
    Tpinch = options.Tpinch;
else
    Tpinch = 5; %[K]
end

% [...] your work


end