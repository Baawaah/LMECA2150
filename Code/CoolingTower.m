function [DAT_WATER DAT_AIR MASSFLOW] = CoolingTower(P_w,options)
% COOLINGTOWER is a cooling tower 0D modelisation
% COOLINGTOWER(P_w,options) compute the thermodynamics states for a Cooling
% tower based on several inputs (given in OPTION) and based on a given 
% water power to dissipate P_w.
% It returns the main results. 
%
% INPUTS :
% P_W = Heat power output at the condenser [W]
% OPTIONS is a structure containing :
%   -options.Tcond  [degres C]: Temperature in the condenser
%   -options.Tpinch [degres C]: Minimum tempearture pinch between Tw_out and the
%                         condenser temperature.
%   -options.Tw_out [degres C]: Cooling water temperature at the condenser outlet
%   -options.Tw_in  [degres C]: Cooling water temperature at the condenser inlet
%   -options.Triver [degres C]: River temperature 
%   -options.Ta_in  [degres C]: Atmospheric air temperature 
%   -options.Ta_out [degres C]: Air outlet temperature of cooling tower 
%   -options.Phi_atm [-]: Relative humidity of atmospheric air
%   -options.Phi_out [-]: Maximum relative humidity of air at the cooling 
%                         tower outlet.
%
% OUTPUT :
% MassFlow [kg/s]: Vector containing the different massflow :
%   -massflow(1) : water massflow at the condenser
%   -massflow(2) : additionnal water massflow = water flow evaporated
%   -massflow(3) : air massflow at the cooling tower   
%
%  dat_water = [T_e1       , T_e2       , T_e3       , T_e4;  %[degres C]
%               h_e1       , h_e2       , h_e3       , h_e4;  %[kJ/kg]
%               m_e1       , m_e2       , m_e3       , m_e4]; %[kg/s]
% 
%  dat_air   = [Ta_in       , Ta_out  ;  %[degres C]
%               ha_in       , ha_out  ;  %[kJ/kg]
%               xa_in       , xa_out  ;  %[kg_water/kg_dry_air]
%               Phia_in     , Phia_out]; %[-] relative humidity
%  
%
% ADDITIONNAL INFORMATIONS
% Water points : 
%       1 : water outlet of cooling tower
%       2 : water just before condenser
%       3 : water just after  condenser
%       4 : water from the river (coming between 1 & 2)
%
% Air points :
%       a_in : air at the cooling tower inlet
%       a_out : air at the cooling tower outlet
%

% Example of how to handle with options structure

if nargin<2
    options=struct();
    if nargin<1
        P_w=200e6; %[W] %200MW_heat
    end
end

if isfield(options,'Tcond')
    Tcond = options.Tcond;
else
    Tpinch = 40;   %[degres C]
end

if isfield(options,'Tpinch')
    Tpinch = options.Tpinch;
else
    Tpinch = 10;   %[degres C] %NB : 5 a 10
end

if isfield(options,'Phi_atm')
    Phi_atm = options.Phi_atm;
else
    Phi_atm = 0.6; %[degres C]
end

if isfield(options,'Phi_out')
    Phi_out = options.Phi_out;
else
    Phi_out = 1;  %[-]
end

if isfield(options,'Ta_in')
    Ta_in = options.Ta_in;
else
    Ta_in = 20;   %[degres C]
end

if isfield(options,'Ta_out')
    Ta_out = options.Ta_out;
else
    Ta_out = 25;  %[degres C]
end

if isfield(options,'Triver')
    Triver = options.Triver;
else
    Triver = 15;  %[degres C]
end

if isfield(options,'Tw_out')
    Tw_out = options.Tw_out;
else
    Tw_out = 38;  %[degres C]
end

if isfield(options,'Tw_in')
    Tw_in = options.Tw_in;
else
    Tw_in = 30;   %[degres C]
end

% ETATS EAU CONNUS
%% entree condenseur : etat 2
T_e2 = Tw_in;                         %[degres C]
p_e2 = 1;                             %[bar]
h_e2 = XSteam('h_pT',p_e2,T_e2)*1e3;  %[J/kg]

%% sortie condenseur : etat 3
T_e3 = Tw_out;                       %[degres C]
p_e3 = 1;                            %[bar]
h_e3 = XSteam('h_pT',p_e3,T_e3)*1e3; %[J/kg]

%% debit d'eau au condenseur : etats 2 = etat 3
m_cond = P_w/(h_e3-h_e2);   %[kg/s]
m_e2 = m_cond;              %[kg/s]
m_e3 = m_cond;              %[kg/s]

%% arrivee eau : riviere, etat 4
T_e4 = Triver;                       %[degres C]
p_e4 = 1;                            %[bar]
h_e4 = XSteam('h_pT',p_e4,T_e4)*1e3; %[J/kg]

% ETATS AIR CONNUS
%% arrivee d'air
Phia_in = Phi_atm;                                                           %[-]
[~,xa_in,~,ha_in,~,~,~] = Psychrometrics('Tdb',Ta_in,'phi',Phia_in*100);     %[-] et [J/kg]

%% sortie d'air 
Phia_out = Phi_out;                                                          %[-]
[~,xa_out,~,ha_out,~,~,~] = Psychrometrics('Tdb',Ta_out,'phi',Phia_out*100); %[-] et [J/kg]

%% Calcul des inconnues restante et du debit d'eau au condenseur

function equations = unknowns(x)
    %x(1) = m_e4; 
    %x(2) = debit_m_air; 
    %x(3) = h_e1; 
    %conservation masse
	equations(1) = x(1) - x(2)*(xa_out-xa_in); 
    %conservation energie volume de controle cooling tower
	equations(2) = x(2)*(ha_out-ha_in)+(m_cond-x(1))*x(3)-m_cond*h_e3;
    %conservation energie volume de controle embranchement river
    equations(3) = x(1)*h_e4+(m_cond-x(1))*x(3)-m_cond*h_e2;
end
results = fsolve(@unknowns,[0 m_cond h_e2],optimoptions('fsolve','Display','off'));
m_e4        = results(1);  %[kg/s]
debit_m_air = results(2);  %[kg/s]
h_e1        = results(3);  %[J/kg]

%% sortie d'eau de la tour
p_e1 = 1;                             %[bar]
T_e1 = XSteam('T_ph',p_e1,h_e1*1e-3); %[degres C]
m_e1 = m_cond - m_e4;                 %[kg/s]

%% Resultats 
MASSFLOW = [m_cond m_e4 debit_m_air];

DAT_WATER = [T_e1      , T_e2       , T_e3       , T_e4;      %[degres C]
             h_e1/1e3  , h_e2/1e3   , h_e3/1e3   , h_e4/1e3;  %[kJ/kg]
             m_e1      , m_e2       , m_e3       , m_e4];     %[kg/s]
 
DAT_AIR   = [Ta_in     , Ta_out;     %[degres C]
             ha_in/1e3 , ha_out/1e3; %[kJ/kg]
             xa_in     , xa_out;     %[kg_water/kg_dry_air]
             Phia_in   , Phia_out];  %[-] relative humidity


end