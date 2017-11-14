function [ETA DATEN DATEX DAT MASSFLOW COMBUSTION] = GT(P_e,options,display)
% GT Gas turbine modelisation
% GT(P_e,options,display) compute the thermodynamics states for a Gas
% turbine based on several inputs (given in OPTION) and based on a given 
% electricity production P_e. It returns the main results. It can as well
% plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_E = electrical power output target [W]
% OPTIONS is a structure containing :
%   -options.k_mec [-] : Shaft losses 
%   -options.T_0   [K] : Reference temperature
%   -options.T_ext [K] : External temperature
%   -options.r     [-] : Comperssion ratio
%   -options.k_cc  [-] : Coefficient of pressure losses due to combustion
%                        chamber
%   -options.T_3   [K] : Temperature after combustion (before turbine)
%   -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for compression
%   -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for expansion
%DISPLAY = 1 or 0. If 1, then the code should plot graphics. If 0, then the
%          do not plot.
%
%OUPUTS : 
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_rotex, compressor-turbine exergy efficiency
%   -eta(6) : eta_combex, Combustion exergy efficiency
% DATEN is a vector with : 
%   -daten(1) : perte_mec [W]
%   -daten(2) : perte_ech [W]
% DATEX is a vector with :
%   -datex(1) : perte_mec [W]
%   -datex(2) : perte_rotex [W]
%   -datex(3) : perte_combex [W]
%   -datex(4) : perte_echex  [W]
% DAT is a matrix containing :
% dat = {T_1       , T_2       , T_3       , T_4; [°C]
%        p_1       , p_2       , p_3       , p_4; [bar]
%        h_1       , h_2       , h_3       , h_4; [kJ/kg]
%        s_1       , s_2       , s_3       , s_4; [kJ/kg/K]
%        e_1       , e_2       , e_3       , e_4;};[kJ/kg]
% MASSFLOW is a vector containing : 
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_c, combustible massflow [kg/s] 
%   -massflow(3) = m_f, exhaust gas massflow [kg/s]
% 
% COMBUSTION is a structure with :
%   -combustion.LHV    : the Lower Heat Value of the fuel [kJ/kg]
%   -combustion.e_c    : the combuistible exergie         [kJ/kg]
%   -combustion.lambda : the air excess                   [-]
%   -combustion.Cp_g   : heat capacity of exhaust gas     [kJ/kg/K]
%   -combustion.fumTG  : is a vector of the exhaust gas composition :
%       -fumTG(1) = m_O2f  : massflow of O2 in exhaust gas [kg/s]
%       -fumTG(1) = m_N2f  : massflow of N2 in exhaust gas [kg/s]
%       -fumTG(1) = m_CO2f : massflow of CO2 in exhaust gas [kg/s]
%       -fumTG(1) = m_H2Of : massflow of H2O in exhaust gas [kg/s] 





% Example of how to handle with options structure
if isfield(options,'T_0')
    data.T_0 = options.T_0;
else
    data.T_0 = 288.15; %[K]
end
if isfield(options,'k_mec')
    data.k_mec = options.k_mec;
else
    data.k_mec = 0.1; 
end
if isfield(options,'T_ext')
    data.T_ext = options.T_ext;
else
    data.T_ext = 288.15; %[K]
end
if isfield(options,'r')
    data.r = options.r;
else
    data.r = 10; 
end
if isfield(options,'k_cc')
    data.k_cc = options.k_cc;
else
    data.k_cc = 0.1; 
end
if isfield(options,'T3')
    data.T3 = options.T3;
else
    data.T3 = 1323.15; 
end
if isfield(options,'eta_PiC')
    data.eta_PiC = options.eta_PiC;
else
    data.eta_PiC = 0.9; 
end
if isfield(options,'eta_PiT')
    data.eta_PiT = options.eta_PiT;
else
    data.eta_PiT = 0.9; 
end
%% THE CODE
% INITIALISATION
MmNO2 = 28,0134;
MmO2  = 31,9988;
MmAr  = 39,948; %on prend le cas de Air = 79% N2 et 21% O2 ??
MmCO2 = 44,0095;
0.7808*MmNO2 + 0.2095*MmO2 + 0.00934*MmAr + 0.0004* MmCO2;
R = 8.3415*1000/28.66;

% 1 2 3 4 5
% p T h s v
%% STATE 1 AIR ENTRY
data.table(1,1) = 1; % [bar]
data.table(1,2) = data.T_ext;
data.table(1,3) = 0; % Etat de référence 
data.table(1,4) = 0; % Etat de référence
data.table(1,5) = 


end