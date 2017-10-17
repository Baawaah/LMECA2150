function [ETA DATEN DATEX DATST DATGT MASSFLOW COMBUSTION] = CCGT2P(P_e,options,display)
% CCGT2P is a Combine cycle Gas Turbine with 2 pressure level
% CCGT2P(P_e,options,display) compute the thermodynamics states for a CCGT
% with 2 pressure level (cfr p164 english reference book) including
% combustion, exchanger and cycles. This is done based on several inputs
% (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input 
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [W]
% OPTIONS is a structure containing :
%   -options.T0       [°C] : Reference temperature
%   -options.T_ext    [K]  : External temperature
%   -options.T_STmax  [°C] : maximum temperature on ST cycle
%   -options.eta_mec  [-] : mecanic efficiency of shafts bearings
%   -options.comb is a structure containing combustion data : 
%       -comb.Tmax    [°C] : maximum combustion temperature
%       -comb.lambda  [-] : air excess
%       -comb.x       [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y       [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
%   -options.plow    [bar]: Low pressure
%   -options.x5       [-] : Vapor ratio [gaseous/liquid] (titre)
%   -option.eta_SiC   [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT   [-] : Isotrenpic efficiency for compression
% DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then 
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1)  : eta_STcyclen, cycle energy efficiency
%   -eta(2)  : eta_GTcyclen, cycle energy efficiency
%   -eta(3)  : eta_toten, overall energy efficiency
%   -eta(4)  : eta_STcyclex, cycle exegy efficiency
%   -eta(5)  : eta_GTcyclex, cycle exegy efficiency
%   -eta(6)  : eta_totex, overall exergie efficiency
%   -eta(7)  : eta_gen, Steam generator energy efficiency
%   -eta(8)  : eta_gex, Steam generator exergy efficiency
%   -eta(9)  : eta_combex, Combustion exergy efficiency
%   -eta(10) : eta_chemex, Chimney exergy efficiency (losses)
%   -eta(11) : eta_transex, Heat exchanger overall exergy efficiency
% MASSFLOW is a vector containing : 
%   -massflow(1) [kg/s]: massflow of steam at 3 (shortcut high pressure)
%   -massflow(2) [kg/s]: massflow of steam at 4 (HP)
%   -massflow(3) [kg/s]: massflow of steam at 5 (all steam flow)
%   -massflow(4) [kg/s]: air massflow at gas turbine inlet 
%   -massflow(5) [kg/s]: combustible massflow
% DATEN is a vector with : 
%   -daten(1) : perte_gen  [W]
%   -daten(2) : perte_mec  [W]
%   -daten(3) : perte_cond [W]
% DATEX is a vector with :
%   -datex(1) : perte_mec    [W]
%   -datex(2) : perte_totex  [W]
%   -datex(3) : perte_rotex  [W]
%   -datex(4) : perte_combex [W]
%   -datex(5) : perte_condex [W]
%   -datex(6) : perte_chemex [W]
%   -datex(7) : perte_transex[W]
% DATST is a matrix containing :
% datST = {T_1       , T_2       , ...       , T_5;  [°C]
%        p_1       , p_2       , ...       , p_5;  [bar]
%        h_1       , h_2       , ...       , h_5;  [kJ/kg]
%        s_1       , s_2       , ...       , s_5;  [kJ/kg/K]
%        e_1       , e_2       , ...       , e_5;  [kJ/kg]
%        x_1       , x_2       , ...       , x_5;};[-]
% DATGT is a matrix containing :
% datGT = {T_1g     , T_2g       , ...       , T_5g;  [°C]
%        p_1g       , p_2g       , ...       , p_5g;  [bar]
%        h_1g       , h_2g       , ...       , h_5g;  [kJ/kg]
%        s_1g       , s_2g       , ...       , s_5g;  [kJ/kg/K]
%        e_1g       , e_2g       , ...       , e_5g;  [kJ/kg]
%        x_1g       , x_2g       , ...       , x_5g;};[-]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combustible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas     [kJ/kg/K]
%   -combustion.fumTG  : is a vector of the exhaust gas composition :
%       -fumTG(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
%       -fumTG(1) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
%       -fumTG(1) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
%       -fumTG(1) = m_H2Of : massflow of H2O in exhaust gas [kg/s] 





% Example of how to handle with options structure

if isfield(options,'T_0')
    T_0 = options.T_0;
else
    T_0 = 288.15; %[K]
end

% [...] your work

end