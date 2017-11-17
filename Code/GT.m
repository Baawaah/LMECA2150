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


if nargin<3
    display = 1;
    if nargin<2
        options = struct();
        if nargin<1
            P_e = 250e6; % [W] Puissance Ã©nergÃ©tique de l'installation
        end
    end
end
% Example of how to handle with options structure
if isfield(options,'T_0')
    data.T_0 = options.T_0;
else
    data.T_0 = 288.15; %[K]
end

if isfield(options,'k_mec')
    data.k_mec = options.k_mec;
else
    data.k_mec = 0.015; 
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
    data.k_cc = 0.95; 
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
%% INITIALISATION
MmNO2 = 28.0134;
MmO2  = 31.9988;
MmAr  = 39.948; %on prend le cas de Air = 79% N2 et 21% O2 ??
MmCO2 = 44.0095;
MmCH4 = 16.04;
R_a = 8.314472*1000/(0.7808*MmNO2 + 0.2095*MmO2 + 0.00934*MmAr + 0.0004* MmCO2);
R_g = 8.314472*1000/(MmCH4);
Cp_a = (0.0004*janaf('c','CO2',300)+0.21*janaf('c','O2',300)+0.78*janaf('c','N2',300))*10^3;
Cv_a = (Cp_a-R_a);
gamma = Cp_a/Cv_a;

% http://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
%     298 - 1300         1300 - 6000
A_g1 = -0.703029;	A_g2 =  85.81217;
B_g1 =	108.4773;	B_g2 =  11.26467;
C_g1 = -42.52157;	C_g2 = -2.114146;
D_g1 =	5.862788;	D_g2 =  0.138190;
E_g1 =	0.678565;	E_g2 = -26.42221;
% F_g1 = -76.84376;	F_g2 = -153.5327;
% G_g1 =	158.7163;	G_g2 =  224.4143;
% H_g1 = -74.87310;	H_g2 = -74.87310;
Cp_g_fun = @(T) (A_g1 + B_g1*(T/1000) + C_g1*(T/1000)^2 + D_g1*(T/1000)^3 + E_g1/((T/1000)^2))/MmCH4*10^3;


%% FUEL
% CxHy - CH4
x = 1; y = 4;
ma1    = 34.39*( (x+y/4) / (3*x+y/4) );
lambda = 1;
%% Pe Power efficiency
P_e

%% STATE 1 AIR ENTRY
data.table(1,1) = 10^5; % [Pa]
data.table(1,2) = data.T_ext;
data.table(1,3) = 0; % Etat de référence 
data.table(1,4) = 0; % Etat de référence
data.table(1,5) = R_a*data.table(1,2)/data.table(1,1);

%% STATE 2 AFTER COMPRESSOR
data.table(2,1) = data.table(1,1)*data.r;
data.table(2,2) = data.table(1,2)*data.r^((gamma-1)/gamma);
data.table(2,3) = 1/data.eta_PiC * gamma/(gamma-1) * R_a * (data.table(2,2)-data.table(1,2));
data.table(2,4) = (1-data.eta_PiC)*Cp_a*log(data.table(2,2)/data.table(1,2));
data.table(2,5) = R_a*data.table(2,2)/data.table(2,1);

%% STATE 3 AFTER COMBUSTION [PRESSURE AND TEMPERATURE]
data.table(3,1) = data.k_cc*data.table(2,1);
data.table(3,2) = data.T3;


%% STATE 4 AFTER TURBINE
data.table(4,1) = 10^5;
data.table(4,2) = data.table(3,2)*(data.table(4,1)/data.table(3,1))^((gamma-1)/gamma);
data.table(4,3) = data.eta_TiC*R_


Cp_g32 = (Cp_g_fun(data.table(3,2))+Cp_g_fun(data.table(2,2)))/2;
Qcomb  = Cp_a*data.table(1,2)*( (1 + 1/(lambda*ma1)) * Cp_g32/Cp_a * data.table(3,2)/data.table(1,2) - data.r^( 1/data.eta_PiC * R_a/Cp_a) );
data.table(3,3) = -1;
 data.table(3,4) = 1;




%% Display
% Plot T-S
figure;
subplot(2,2,1);
hold on;
for i = 1 : length(data.table(:,1))
plot(data.table(i,4),data.table(i,2),'b*');
end
hold off;
% Plot P-v
subplot(2,2,2);
hold on;
for i = 1 : length(data.table(:,1))
plot(data.table(i,5),data.table(i,1),'b*');
end
hold off;
%%
format short;
disp('    p[bar]     T[K]      h[KJ]         s         v')
disp([data.table(:,1)/10^5 data.table(:,2) data.table(:,3)/10^3 data.table(:,4) data.table(:,5)])
end