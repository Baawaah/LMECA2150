function [ETA,DATEN,DATEX,DAT,MASSFLOW,COMBUSTION] = GT(P_e,options,display)
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
%   -options.r     [-] : Comperession ratio
%   -options.k_cc  [-] : Coefficient of pressure losses due to combustion
%                        chamber
%   -options.T_3   [K] : Temperature after combustion (before turbine)
%   -option.eta_PiC[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for compression
%   -option.eta_PiT[-] : Intern polytropic efficiency (Rendement
%                        polytropique interne) for expansion
%   -options.comb is a structure containing combustion data : 
%       -comb.Tmax    [°C] : maximum combustion temperature
%       -comb.lambda  [-] : air excess
%       -comb.x       [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y       [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
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
            P_e = 150e6; % [W] Puissance Ã©nergÃ©tique de l'installation
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
    data.k_cc = 0.90; 
end
if isfield(options,'T3')
    data.T3 = options.T3;
else
    data.T3 = 1050+273.15; 
end
if isfield(options,'eta_PiC')
    data.eta_PiC = options.eta_PiC;
else
    data.eta_PiC = 0.90; 
end
if isfield(options,'eta_PiT')
    data.eta_PiT = options.eta_PiT;
else
    data.eta_PiT = 0.90; 
end
% OPTION RAJOUTER
%   -options.comb is a structure containing combustion data : 
%       -comb.Tmax    [°C] : maximum combustion temperature
%       -comb.lambda  [-] : air excess
%       -comb.x       [-] : the ratio O_x/C. Example 0.05 in CH_1.2O_0.05
%       -comb.y       [-] : the ratio H_y/C. Example 1.2 in CH_1.2O_0.05
if isfield(options,'comb')
    data.T3 = options.comb.Tmax+273.15;
    lambda  = options.comb.lambda;
    x = options.comb.x;
    y = options.comb.y;    
else
    x = 1; 
    y = 4;
    lambda = 3.5;
end

%% THE CODE
%% INITIALISATION
%% AIR
MmNO2 = 28.0134;
MmO2  = 31.9988;
MmAr  = 39.948; %on prend le cas de Air = 79% N2 et 21% O2 ??
MmCO2 = 44.0095;
MmH2O = 18;
R_a = 8.314472*1000/(0.7808*MmNO2 + 0.2095*MmO2 + 0.00934*MmAr + 0.0004* MmCO2);
Cp_a = (0.0004*janaf('c','CO2',300)+0.21*janaf('c','O2',300)+0.78*janaf('c','N2',300))*10^3;
Cp_a_fun = @(T) (0.0004*janaf('c','CO2',T)+0.21*janaf('c','O2',T)+0.78*janaf('c','N2',T))*10^3;
Cv_a = (Cp_a-R_a);
gamma_a = Cp_a/Cv_a;
%h_a_fun = @(T) (0.0004*janaf('h','CO2',T)+0.21*janaf('h','O2',T)+0.78*janaf('h','N2',T))*10^3;

%% FUEL
% CxHy - CH4

ma1    = 34.39*( (x+y/4) / (3*x+y/4) );

MmCH4 = 16.04;
%HHV = 55.50*10^3; %[KJ]
LHV = 50.00*10^3; %[KJ/KG]
%H_c_298 = -74.87;
% Cp of the Fuel
% http://webbook.nist.gov/cgi/cbook.cgi?ID=C74828&Mask=1
%     298 - 1300         1300 - 6000
A_g1 = -0.703029;	%A_g2 =  85.81217;
B_g1 =	108.4773;	%B_g2 =  11.26467;
C_g1 = -42.52157;	%C_g2 = -2.114146;
D_g1 =	5.862788;	%D_g2 =  0.138190;
E_g1 =	0.678565;	%E_g2 = -26.42221;
%F_g1 = -76.84376;	%F_g2 = -153.5327;
%G_g1 =	158.7163;	%G_g2 =  224.4143;
%H_g1 = -74.87310;	%H_g2 = -74.87310;
Cp_c_fun = @(T) (A_g1 + B_g1*(T/1000) + C_g1*(T/1000)^2 + D_g1*(T/1000)^3 + E_g1/((T/1000)^2))/MmCH4*10^3;
%h_c_fun = @(T) (A_g1*(T/1000) + B_g1*((T/1000)^2)/2 + C_g1*((T/1000)^3)/3 + D_g1*((T/1000)^4)/4 - E_g1/(T/1000) + F_g1 - H_g1 ); 

e_c = LHV*1.041*10^3;

R_c = 8.314472*1000/(MmCH4);
%% Gas
R_g = 1/(1+lambda*ma1) *R_c + lambda*ma1/(1+lambda*ma1) *R_a;
Cp_g_fun = @(T) 1/(1+lambda*ma1)*Cp_c_fun(T)+ lambda*ma1/(1+lambda*ma1)*Cp_a_fun(T);
Cv_g_fun = @(T) Cp_g_fun(T)-R_g;
gamma_g_fun= @(T) Cp_g_fun(T)/Cv_g_fun(T);
%figure;
%fplot(gamma_g_fun,[300 1500])
%% STATE 1 AIR ENTRY
data.table(1,1) = 10^5; % [Pa]
data.table(1,2) = data.T_ext;
data.table(1,3) = 0; % Etat de référence 
data.table(1,4) = 0; % Etat de référence
data.table(1,5) = R_a*data.table(1,2)/data.table(1,1);
data.table(1,6) = data.table(1,3) - data.T_0*data.table(1,4); 
%% STATE 2 AFTER COMPRESSOR
data.table(2,1) = data.table(1,1)*data.r;
data.table(2,2) = data.table(1,2)*data.r^((gamma_a-1)/gamma_a);
WmC = 1/data.eta_PiC * gamma_a/(gamma_a-1) * R_a * (data.table(2,2)-data.table(1,2));
data.table(2,3) = WmC + data.table(1,3); 
data.table(2,4) = (1-data.eta_PiC)*Cp_a*log(data.table(2,2)/data.table(1,2) + data.table(1,4));
data.table(2,5) = R_a*data.table(2,2)/data.table(2,1);
data.table(2,6) = data.table(2,3) - data.T_0*data.table(2,4); 
%% STATE 3 AFTER COMBUSTION [PRESSURE AND TEMPERATURE]
data.table(3,1) = data.k_cc*data.table(2,1);
data.table(3,2) = data.T3;
    Cp_g32 = (Cp_g_fun(data.table(3,2))+Cp_g_fun(data.table(2,2)))/2;
    Qcomb  = Cp_a*data.table(1,2)*( (1 + 1/(lambda*ma1)) * Cp_g32/Cp_a * data.table(3,2)/data.table(1,2) - data.r^( 1/data.eta_PiC * R_a/Cp_a) );
data.table(3,3) = Qcomb + data.table(2,3);
data.table(3,4) = Cp_g32*log(data.table(3,2)/data.table(2,2)) + data.table(2,4);
data.table(3,5) = R_g*data.table(3,2)/data.table(3,1);
data.table(3,6) = data.table(3,3) - data.T_0*data.table(3,4); 
%% STATE 4 AFTER TURBINE
data.table(4,1) = 101325;
    T4_guess     = data.table(3,2)*(data.table(4,1)/data.table(3,1))^((gamma_a-1)/gamma_a);
    T4_opt_guess = data.table(3,2)*(data.table(4,1)/data.table(2,1))^((gamma_a-1)/gamma_a);
    T4_error = 1;
    while T4_error > 10^-6
        T4_old = T4_guess;
        gamma_guess = (gamma_g_fun(T4_guess)+gamma_g_fun(data.table(3,2)))/2;      
        T4_guess = data.table(3,2)*(data.table(4,1)/data.table(3,1))^((gamma_guess-1)/gamma_guess);
        T4_opt_guess = data.table(3,2)*(data.table(4,1)/data.table(2,1))^((gamma_guess-1)/gamma_guess);
        T4_error = abs(T4_guess-T4_old);
    end
data.table(4,2) = T4_guess;
    %Cp_g43 = (Cp_g_fun(data.table(4,2))+Cp_g_fun(data.table(3,2)))/2;
    %WmT = Cp_g43*(data.table(3,2)-data.table(4,2))
    WmT= 1/data.eta_PiT * gamma_g_fun(data.table(4,2))/(gamma_g_fun(data.table(4,2))-1) * R_g * (data.table(4,2)-data.table(3,2));
    WmT_opt= 1/data.eta_PiT * gamma_g_fun(T4_opt_guess)/(gamma_g_fun(T4_opt_guess)-1) * R_g * (T4_opt_guess-data.table(3,2));
data.table(4,3) = data.table(3,3) + WmT ;
data.table(4,4) = (1-data.eta_PiT)*Cp_g_fun(data.table(4,2))*log(data.table(3,2)/data.table(4,2)) + data.table(3,4);
data.table(4,5) = R_g*data.table(4,2)/data.table(4,1);
data.table(4,6) = data.table(4,3) - data.T_0*data.table(4,4); 
%% MASS FLOW
% MASSFLOW is a vector containing : 
%   -massflow(1) = m_a, air massflow [kg/s]
%   -massflow(2) = m_c, combustible massflow [kg/s] 
%   -massflow(3) = m_f, exhaust gas massflow [kg/s]
 %eta_mec = (1-data.k_mec);
 %Pfmec = P_e *(1/eta_mec - 1); 
 %Pm = P_e + Pfmec; 
 %mdot_g = abs(Pm/(WmT + lambda*ma1/(1+lambda*ma1)*WmC));
 mdot_g = P_e/( (1-data.k_mec)*abs(WmT) - (1+data.k_mec)*WmC*(lambda*ma1/(1+lambda*ma1)));
 mdot_a =  mdot_g * lambda*ma1/(1+lambda*ma1);
 mdot_c = mdot_g * 1/(1+lambda*ma1);
 Pm = abs(mdot_g*WmT) - abs(mdot_a*WmC);
 Pfmec = data.k_mec * ( abs(mdot_g*WmT) + abs(mdot_a*WmC) ) ;
 Wm = Pm/mdot_a;
 eta_mec = 1- Pfmec / Pm;
 
 MASSFLOW(1) = mdot_a;
 MASSFLOW(2) = mdot_c;
 MASSFLOW(3) = mdot_g;
%% LOSSES
%OUPUTS : 
% ETA is a vector with :
%   -eta(1) : eta_cyclen, cycle energy efficiency
%   -eta(2) : eta_toten, overall energy efficiency
%   -eta(3) : eta_cyclex, cycle exegy efficiency
%   -eta(4) : eta_totex, overall exergie efficiency
%   -eta(5) : eta_rotex, compressor-turbine exergy efficiency
%   -eta(6) : eta_combex, Combustion exergy efficiency
ETA(1) = Wm/Qcomb;
ETA(2) = ETA(1)*(1-data.k_mec);
ETA(3) = Pm/(mdot_g*data.table(3,6) - mdot_a*data.table(2,6) );
ETA(4) = P_e/(mdot_c*e_c);
ETA(5) = (mdot_g*(data.table(3,3)-data.table(4,3))-mdot_a*(data.table(2,3)-data.table(1,3)))/(mdot_g*(data.table(3,6)-data.table(4,6))-mdot_a*(data.table(2,6)-data.table(1,6)));
ETA(6) = (mdot_g*data.table(3,6)-mdot_a*data.table(2,6))/(mdot_c*e_c);
%% COMBUSTION
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
COMBUSTION.LHV      = LHV;
COMBUSTION.e_c      = e_c;
COMBUSTION.lambda   = lambda;
COMBUSTION.Cp_g     = Cp_g_fun(data.T3);
CH4_mole = mdot_c/ MmCH4; 
COMBUSTION.fumTG(1) = mdot_a*0.2095*(lambda-1)/lambda;
COMBUSTION.fumTG(2) = mdot_a*0.7808;
COMBUSTION.fumTG(3) = CH4_mole*MmCO2+mdot_a*0.0004;
COMBUSTION.fumTG(3) = CH4_mole*MmH2O*2;
%% DATA ENER
% DATEN is a vector with : 
%   -daten(1) : perte_mec [W]
%   -daten(2) : perte_ech [W]
DATEN(1) = Pfmec;
DATEN(2) = abs(WmT_opt)-abs(WmT);

%% DATA EXER
% DATEX is a vector with :
%   -datex(1) : perte_mec [W]
%   -datex(2) : perte_rotex [W]
%   -datex(3) : perte_combex [W]
%   -datex(4) : perte_echex  [W]
DATEX(1) = 0;
DATEX(2) = (data.table(2,6)) * MASSFLOW(1) - (data.table(4,6)-data.table(3,6)) * MASSFLOW(3);
DATEX(3) = (data.table(3,6)-data.table(2,6)) * MASSFLOW(2);
DATEX(4) = 0;
%% DATA OVERALL
% DAT is a matrix containing :
% dat = {T_1       , T_2       , T_3       , T_4; [°C]
%        p_1       , p_2       , p_3       , p_4; [bar]
%        h_1       , h_2       , h_3       , h_4; [kJ/kg]
%        s_1       , s_2       , s_3       , s_4; [kJ/kg/K]
%        e_1       , e_2       , e_3       , e_4;};[kJ/kg]
TEMP = data.table(:,2) - 273.15;
PRES = data.table(:,1)/10^5;
HENT = data.table(:,3)/10^3;
ENTR = data.table(:,4)/10^3;
EXER = data.table(:,6)/10^3;
DAT = [TEMP'; PRES'; HENT';ENTR';EXER'];
%% Display
if display == 1
    % Plot T-S
    gamma_g = gamma_g_fun(data.T3);
    %curve1 = @(S) exp( (gamma_g)/( gamma_g-1) * (S) * 1/R_g  ) +data.table(2,2)
    curve1 = @(S) exp( ( S*10^3 - data.table(2,4)  )/((1-data.eta_PiC)*Cp_g32) ) * data.table(2,2); 
    curve2 = @(S) exp( ( S*10^3 - data.table(2,4)  )/Cp_g32 ) * data.table(2,2);
    curve3 = @(S) exp( -( S*10^3 - data.table(3,4)  )/((1-data.eta_PiT)*Cp_g_fun(data.table(4,2))) ) * data.table(3,2);
    curve4 = @(S) exp( ( S*10^3 - data.table(4,4)  )/Cp_g32 ) * data.table(4,2);
    figure;
    hold on;
    for i = 1 : length(data.table(:,1))
        plot(data.table(i,4)/10^3,data.table(i,2),'b^');
    end
    %XC1  = linspace(0,500,100);
    %XC1Y = curve1(XC1);
    %plot(XC1/10^3,XC1Y)
    fplot(curve1,[data.table(1,4)/10^3,data.table(2,4)/10^3])
    fplot(curve2,[data.table(2,4)/10^3,data.table(3,4)/10^3])
    fplot(curve3,[data.table(3,4)/10^3,data.table(4,4)/10^3])
    fplot(curve4,[data.table(4,4)/10^3,data.table(1,4)/10^3])
    %axis([0,1.5,0,inf])
    hold off;
    grid on;
    grid minor;
    title('Gas Turbine T-S');
    xlabel('Entropy s[kJ/kgK]');
    ylabel('Temperature T [K]');
    
    figure;
    sum_energ = abs(diff(data.table(:,3)))
    label_en = {'Compressor','Combustion Chamber','Turbine'}
    pie(sum_energ,label_en)
    
    figure;
    sum_exerg = abs(diff(data.table(:,6)))
    sum_exerg(1) = sum_exerg(1) * MASSFLOW(1)
    sum_exerg(2) = sum_exerg(2) * MASSFLOW(3)
    sum_exerg(3) = sum_exerg(2) * MASSFLOW(3)
    label_ex = {'Compressor','Combustion Chamber','Turbine',}
    pie(sum_energ,label_en)
    
%     figure;
%     hold on;
%     for i = 1 : length(data.table(:,1))
%         plot(data.table(i,4),data.table(i,3),'b^');
%     end
%     title('h-s');
%     hold off;

    %%
    format short g;
    disp('    p[bar]     T[K]      h[MJ]     s[KJ/K]        v      e[KJ]');
    disp([data.table(:,1)/10^5 data.table(:,2) data.table(:,3)/10^3 data.table(:,4)/10^3 data.table(:,5) data.table(:,6)/10^3]);
    disp('  WmC[KJ]   WmT[KJ]   Wm[KJ]    mdot_g     mdot_c   mdot_a');
    disp([WmC/10^3 WmT/10^3 Wm/10^3 mdot_g mdot_c mdot_a]);
    disp('eta_cyclen eta_toten eta_cyclex eta_totex eta_rotex eta_combex eta_mec');
    disp([ETA(1) ETA(2) ETA(3) ETA(4) ETA(5) ETA(6) eta_mec]);
end
end