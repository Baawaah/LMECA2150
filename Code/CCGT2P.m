function [ETA MASSFLOW FIG] = CCGT2P(P_e,options,display)
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
% P_E  = electrical power output target for CC turbine [W]
% P_EG = electrical power output target for gas turbine [W]
% OPTIONS is a structure containing :
%   -options.T0       [°C]  : Reference temperature
%   -options.T_ext    [K]   : External temperature
%   -options.T_STmax  [°C]  : maximum temperature on ST cycle
%   -options.eta_mec  [-]   : mecanic efficiency of shafts bearings
%   -options.plow     [bar] : Low pressure
%   -options.x5       [-]   : Vapor ratio [gaseous/liquid] (titre)
%   -options.eta_SiC  [-]   : Isotrenpic efficiency for compression
%   -options.eta_SiT  [-]   : Isotrenpic efficiency for expansion
%   -options.GT    [struct] : options for Gas turbine (see GT function) 
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
%   -massflow(1) [kg/s]: water massflow at high pressure turbine inlet
%   -massflow(2) [kg/s]: water massflow at low pressure turbine inlet
%   -massflow(3) [kg/s]: air massflow at gas turbine inlet 
%   -massflow(4) [kg/s]: combustible massflow
%% INPUT HANDLER
if nargin<3, display = 1; if nargin<2, options = struct(); if nargin<1, P_e = 160e6; end, end, end

if isfield(options,'T_0')       , T_0 = options.T_0;          else T_0 = 15;         end %   -options.T0       [°C] : Reference temperature
if isfield(options,'T_ext')     , T_ext = options.T_ext;      else T_ext = 26;       end %   -options.T_ext    [K]  : Condenseur cold outlet temperature
if isfield(options,'TpinchCond'), TpinchCond = options.T_ext; else TpinchCond = 4;   end %   -options.TpinchCond[°C] : Temperature pinch at condenser 
if isfield(options,'T_STMax')   , T_STMax = options.T_STMax;  else T_STMax = 525;    end %   -options.T_STmax  [°C] : maximum temperature on ST cycle
if isfield(options,'eta_mec')   , eta_mec = options.eta_mec;  else eta_mec = 0.9;    end %   -options.eta_mec  [-]  : mecanic efficiency of shafts bearings
if isfield(options,'plow')      , plow = options.plow;        else plow = 5.8;       end %   -options.plow    [bar] : Low pressure
if isfield(options,'phig')      , phig = options.phig;        else phig = 78;        end %   -options.phig    [bar] : High pressure
if isfield(options,'x5')        , x5 = options.x5;            else x5 = 1;           end %   -options.x5       [-]  : Vapor ratio [gaseous/liquid] (titre)
if isfield(options,'eta_SiC')   , eta_SiC = options.eta_SiC;  else eta_SiC = 0.9;    end %   -option.eta_SiC    [-] : Isotrenpic efficiency for compression
if isfield(options,'eta_SiT')   , eta_SiT = options.eta_SiT;  else eta_SiT(1) = 0.9;
                                                                   eta_SiT(2) = 0.9; end %   -option.eta_SiT   [-]  : Isotrenpic efficiency for compression
if isfield(options,'GT')        , optionsGT = options.GT;     else optionsGT.comb.Tmax = 1150;
                                                                   optionsGT.comb.x = 1; 
                                                                   optionsGT.comb.y = 4;
                                                                   optionsGT.comb.lambda = 3; end 
                                                                                         %   -options.GT    [struct] : options for Gas turbine (see GT function)    
%% GAS TURBINE COMPUTATION
% REQUIRE GT.m
% GT(P_e,options,display) compute the thermodynamics states for a Gas
% turbine based on several inputs (given in OPTION) and based on a given 
% electricity production P_E_G. It returns the main results. It can as well
% plots graphs if input argument DISPLAY = true (<=> DISPLAY=1)
P_E_G = 150*10^6;
[GT_ETA GT_DATEN GT_DATEX GT_DAT GT_MASSFLOW GT_COMBUSTION] = GT(P_E_G,optionsGT,0);
%% INITIALISATION OF VALUE
% Cycle Table Formalisme
% [ 1 2 3 4 5 6 7 ]
% [ p T h s v e x ]
TpinchChim = 10;               % Heat exchanger of the exit GT
T_pinch_EVA_LP = 139; 
Pump_comp = [250 13];          % Compression Ratio
Cond_loss = 0.1;               % Condensor Pressure Loss
h0_273 = XSteam('h_pT',1,T_0); % Enthalpy0 at 273
s0_273 = XSteam('s_pT',1,T_0); % Entropy0 at 273
k_hp = 0.5;                    % Ratio of HP LP water fraction

%% STATE AFTER HEATING
table(3,1) = phig; %[bar] 
table(3,2) = min( T_STMax,GT_DAT(1,4)-TpinchChim); %[C]
table(3,3) = XSteam('h_pT',table(3,1),table(3,2));
table(3,4) = XSteam('s_pT',table(3,1),table(3,2));
table(3,5) = XSteam('v_pT',table(3,1),table(3,2));
table(3,6) = table(3,3)-h0_273 - (T_0+273.15)*(table(3,4)-s0_273); 
%% STATE AFTER HP Turbine
table(4,1) = plow;
    h4s = XSteam('h_ps',table(4,1),table(3,4));
table(4,3) = eta_SiT(1)*h4s +  ( 1-eta_SiT(1) ) * table(3,3); 
table(4,2) = XSteam('T_ph',table(4,1),table(4,3));
table(4,4) = XSteam('s_ph',table(4,1),table(4,3));
table(4,5) = XSteam('v_ph',table(4,1),table(4,3));
table(4,6) = table(4,3)-h0_273 - (T_0+273.15)*(table(4,4)-s0_273); 
WmT_HP = table(4,3)-table(3,3); 
%% STATE AFTER LP Turbine
table(5,2) = T_ext+TpinchCond+3;
table(5,1) = XSteam('psat_T',table(5,2));
table(5,7) = x5;
table(5,3) = XSteam('h_px',table(5,1),table(5,7));
table(5,4) = XSteam('s_ph',table(5,1),table(5,3));
table(5,5) = XSteam('v_ph',table(5,1),table(5,3));
table(5,6) = table(5,3)-h0_273 - (T_0+273.15)*(table(5,4)-s0_273); 
WmT_LP = table(5,3)-table(4,3);
%% STATE AFTER CONDENSOR
table(1,1) = table(5,1)* (1-Cond_loss);
table(1,2) = T_ext+TpinchCond;
table(1,3) = XSteam('h_pT',table(1,1),table(1,2));
table(1,4) = XSteam('s_pT',table(1,1),table(1,2));
table(1,5) = XSteam('v_pT',table(1,2),table(1,2));
table(1,6) = table(1,3)-h0_273 - (T_0+273.15)*(table(1,4)-s0_273); 
%% STATE AFTER PUMP
table(2,1) = table(1,1)*Pump_comp(1); 
    h2s = XSteam('h_ps',table(2,1),table(1,4));
table(2,3) = 1/eta_SiC * h2s + (1-1/eta_SiC)*table(1,3);    
table(2,2) = XSteam('T_ph',table(2,1),table(2,3));
table(2,4) = XSteam('s_ph',table(2,1),table(2,3));
table(2,5) = XSteam('v_ph',table(2,1),table(2,3));
table(2,6) = table(2,3)-h0_273 - (T_0+273.15)*(table(2,4)-s0_273); 
WmC = table(2,3) - table(1,3)
%% COMPUTATION OF THE MASSFLOW
Wm = k_hp*abs(WmT_HP)+abs(WmT_LP) - WmC;
mdot_w = P_e/(Wm*10^3 /(eta_mec));
mdot_w_HP = k_hp*mdot_w; 
mdot_w_LP = (1-k_hp)*mdot_w;
%% Chimney Part
% 6 - EXIT ECONOMIZER
table(6,1) = table(2,1);
table(6,2) = XSteam('Tsat_p',plow);
table(6,3) = XSteam('H_pt',table(6,1),table(6,2));
table(6,4) = XSteam('s_pt',table(6,1),table(6,2));
table(6,5) = XSteam('v_pt',table(6,1),table(6,2));
table(6,6) = table(6,3)-h0_273 - (T_0+273.15)*(table(6,4)-s0_273);
% 7 - EXIT LP BALLON LIQUID AND PUMP
table(7,1) = phig;
    h7s = XSteam('h_ps',table(7,1),table(6,4));
table(7,3) = 1/eta_SiC * h7s + (1-1/eta_SiC)*table(6,3);
table(7,2) = XSteam('T_ph',table(7,1),table(7,3));
table(7,4) = XSteam('s_pt',table(7,1),table(7,2));
table(7,5) = XSteam('v_pt',table(7,1),table(7,2));
table(7,6) = table(7,3)-h0_273 - (T_0+273.15)*(table(7,4)-s0_273);
% 8 - EXIT LP BALLON VAPOR
table(8,1) = plow;
table(8,2) = XSteam('Tsat_p',plow);
table(8,3) = XSteam('HV_p',table(8,1));
table(8,4) = XSteam('sV_p',table(8,1));
table(8,5) = XSteam('vV_p',table(8,1));
table(8,6) = table(8,3)-h0_273 - (T_0+273.15)*(table(8,4)-s0_273);
% 9 - EXIT ECO HP
table(9,1) = phig;
table(9,2) = XSteam('Tsat_p',phig);
table(9,3) = XSteam('HL_p',table(9,1));
table(9,4) = XSteam('sL_p',table(9,1));
table(9,5) = XSteam('vL_p',table(9,1));
table(9,6) = table(9,3)-h0_273 - (T_0+273.15)*(table(9,4)-s0_273);
% 10 - EXIT HP BALLON VAPOR
table(10,1) = phig
table(10,2) = XSteam('Tsat_p',phig);
table(10,3) = XSteam('HV_p',table(10,1));
table(10,4) = XSteam('sV_p',table(10,1));
table(10,5) = XSteam('vV_p',table(10,1));
table(10,6) = table(10,3)-h0_273 - (T_0+273.15)*(table(10,4)-s0_273);
%% EXHAUST GT GAS
Cp_g   = GT_COMBUSTION.Cp_g;
mdot_g = GT_MASSFLOW(3);
T_pinch_EVA_HP = 5;
T_pinch_EVA_LP = 10; 
% ENTRY OF THE CHIMNEY
table(11,1) = GT_DAT(2,4);
table(11,2) = GT_DAT(1,4);
table(11,3) = GT_DAT(3,4);
table(11,4) = GT_DAT(4,4);
table(11,6) = GT_DAT(5,4);
% AFTER SUP HP
table(12,1) = table(11,1);
table(12,3) = table(11,3) - mdot_w_HP/mdot_g  * (table(3,3)-table(10,3));
table(12,2) = table(11,2) + (table(12,3)-table(11,3))/Cp_g;
table(12,4) = Cp_g*log((table(12,2)+273.15)/(table(11,2)+273.15)) + table(11,4);
table(12,6) = table(12,3) - (T_0+273.15)*table(12,4); 
% AFTER EVA HP
table(13,1) = table(12,1);
table(13,2) = XSteam('Tsat_p',phig) + T_pinch_EVA_HP;
table(13,3) = table(12,3) + Cp_g*(table(13,2)-table(12,2));
table(13,4) = Cp_g*log((table(13,2)+273.15)/(table(12,2)+273.15)) + table(12,4);
table(13,6) = table(13,3) - (T_0+273.15)*table(13,4); 
% AFTER ECO HP SUP LP
table(14,1) = table(13,1);
table(14,2) = (table(7,2)+table(8,2))/2 + T_pinch_EVA_HP;
table(14,3) = table(13,3) + Cp_g*(table(14,2)-table(13,2));
table(14,4) = Cp_g*log((table(14,2)+273.15)/(table(13,2)+273.15)) + table(13,4);
table(14,6) = table(14,3) - (T_0+273.15)*table(14,4); 
% AFTER EVA LP
table(15,1) = table(14,1);
table(15,2) = table(6,2) + T_pinch_EVA_HP;
table(15,3) = table(14,3) + Cp_g*(table(15,2)-table(14,2));
table(15,4) = Cp_g*log((table(15,2)+273.15)/(table(14,2)+273.15)) + table(14,4);
table(15,6) = table(15,3) - (T_0+273.15)*table(15,4); 
% AFTER ECO LP - EXIT CHIMNEY 
table(16,1) = table(15,1);
table(16,3) = table(15,3) - mdot_w/mdot_g*(table(6,2)-table(2,2));
table(16,2) = table(15,2) + (table(16,3)-table(15,3))/Cp_g;
table(16,4) = Cp_g*log((table(16,2)+273.15)/(table(15,2)+273.15)) + table(15,4);
table(16,6) = table(16,3) - (T_0+273.15)*table(16,4); 
%% CYCLE HP AND LP MASSFLOW
h_EVA_HP = XSteam('HV_p',phig)-XSteam('HL_p',phig);
h_EVA_LP = XSteam('HV_p',plow)-XSteam('HL_p',plow);
mdot_v_HP = (table(12,3)-table(13,3))/(h_EVA_HP) * mdot_g;
mdot_v_LP = (table(14,3)-table(15,3))/(h_EVA_LP) * mdot_g;

%% ETA Efficiency
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

%ST Cyclen
%P_chm = mdot_g * (table_chim(1,2)-table_chim(6,2));
P_ST  = mdot_w*Wm;
ETA(1)= P_ST/(GT_COMBUSTION.LHV * GT_MASSFLOW(2));%P_chm;
%GT Cyclen
ETA(2)= GT_ETA(1);
%CC Toten
Qech = mdot_g*GT_DAT(3,4) - GT_MASSFLOW(1)* GT_DAT(3,1);
epsilon_ech = Qech/(GT_COMBUSTION.LHV * GT_MASSFLOW(2));
ETA(3)= ETA(1) + ETA(2) - ETA(1)*ETA(2) - epsilon_ech*ETA(2);
% ST Cyclex
ETA(4) = P_ST/( (k_hp*mdot_w)*(table(3,6)-table(2,6)) + ((1-k_hp)*mdot_w)*(table(4,6)-table(2,6)));
% GT Cyclex
ETA(5) = GT_ETA(3);
% ETA TOTEX
ETA(6) = (P_e+P_E_G)/(GT_MASSFLOW(2)*GT_COMBUSTION.e_c);
% ETA GEN
% ETA(7) = ( (k_hp*mdot_w)*(table(3,3)-table(2,3)) + ((1-k_hp)*mdot_w)*(table(4,3)-table(2,3)) )/(GT_MASSFLOW(2)*GT_COMBUSTION.LHV) 
% ETA GEX
%ETA(8) = 
%% OUTPUT MASSFLOW
% MASSFLOW is a vector containing : 
MASSFLOW(1) = mdot_w_HP;      %   -massflow(1) [kg/s]: water massflow at high pressure turbine inlet
MASSFLOW(2) = mdot_w_LP;      %   -massflow(2) [kg/s]: water massflow at low pressure turbine inlet
MASSFLOW(3) = mdot_w ;        %   -massflow(3) [kg/s]: air massflow at gas turbine inlet 
MASSFLOW(4) = GT_MASSFLOW(2); %   -massflow(4) [kg/s]: combustible massflow

%% DISPLAY
if display == 1, visibility = 'on'; , else visibility = 'off'; , end
%% FIG 1 - CYCLE T-S DIAGRAM
    FIG(1) = figure('visible',visibility);
    hold on;
     T=linspace(0,375,1000);  
     for i= 1 : length(T)
          S1(i)=XSteam('sL_T',T(i));
          S2(i)=XSteam('sV_T',T(i));
     end   
     S=[S1,S2];
     T=[T,T];
    plot(S,T,'b:');
    for i = 1 : length(table(:,1))
        plot(table(i,4),table(i,2),'k^');
        text(table(i,4)+0.1,table(i,2)+0.1,sprintf('%d',i));
        %if i <= length(table(:,1))-1, plot([table(i,4) table(i+1,4)], [table(i,2) table(i+1,2)]), end
        if i <= 4 && i ~= 2, plot([table(i,4) table(i+1,4)], [table(i,2) table(i+1,2)],'b'), end
    end
    %plot([table(length(table(:,1)),4) table(1,4)], [table(length(table(:,1)),2) table(1,2)])
    plot([table(5,4) table(1,4)], [table(5,2) table(1,2)])
    % CHIMNEE
    plot([table(2,4) table(6,4)], [table(2,2) table(6,2)],'r--')
    plot([table(6,4) table(7,4)], [table(6,2) table(7,2)],'g--')
    plot([table(6,4) table(8,4)], [table(6,2) table(8,2)],'g--')
    plot([table(8,4) table(4,4)], [table(8,2) table(4,2)],'g--')
    plot([table(7,4) table(9,4)], [table(7,2) table(9,2)],'r--')
    plot([table(9,4) table(10,4)], [table(9,2) table(10,2)],'r--')
    plot([table(10,4) table(3,4)], [table(10,2) table(3,2)],'r--')
    
    grid on;
    grid minor;
    hold off;
    title('CCGT2P T-S Diagram ');
    xlabel('Entropy [kJ/C°]');
    ylabel('Temperature [C°]');
    %visibility = 'off';
%%  FIG 2 - WATER REPARTITION  
    FIG(2)= figure('visible',visibility);
    F2labels = {'T-WMF';'H- WMF';'LP-WMF';'HP-V-WMF';'LP-V-WMF'};
    bar([1:5],[mdot_w mdot_w_HP mdot_w_LP mdot_v_HP mdot_v_LP ]);
    set(gca,'xticklabel',F2labels)
    colormap parula
%%  FIG 3 - ENERGETIC POWER
    FIG(3)= figure('visible',visibility);
    LOSS_MECA = (P_e+P_E_G)/1e3*(1-eta_mec);
    LOSS_COND = abs((table(5,3)-table(1,3)))*mdot_w;
    LOSS_CHIM = abs((table(11,3)-table(16,3)))*GT_MASSFLOW(3);
    ENPO = [P_e/1e3 P_E_G/1e3 LOSS_MECA LOSS_COND LOSS_CHIM];
    F3labels = {'Effective Power ST';'Effective Power GT';'Mecanic Losses';'Condensor Losses';'Chimnee Losses'};
    pie(ENPO,F3labels);
    colormap summer;
%%  FIG 4 - EXERGY FLOW
    FIG(3)= figure('visible',visibility);
    LOSS_COMB = GT_DATEX(3);
    EXFL = [P_e/1e3 P_E_G/1e3 LOSS_MECA LOSS_COMB];
    pie(EXFL);
    colormap winter;
    
    
    format short g;
    disp('        STATE       p[bar]         T[C]         h[MJ]     s[KJ/C]            v        e[KJ]');
    disp([linspace(1,length(table(:,1)),length(table(:,1)))' table(:,1) table(:,2) (table(:,3)/10^3) (table(:,4)) table(:,5) table(:,6) ]);
    %     format short g;
%     disp('Tgas[C]       hgas       p[bar]     T[C]      h[MJ]     s[J/K]        v      x[KJ]');
%     disp([table_chim(:,1) table_chim(:,2) table_chim(:,3) table_chim(:,4) table_chim(:,5)/10^3  ]);

    
    disp('  WmC[KJ]    WmT[KJ]    Wm[KJ]    mdot_g    mdot_w    mdot_w_HP  ');
    disp([WmT_HP WmT_LP  Wm mdot_g  mdot_w_HP]);

end
