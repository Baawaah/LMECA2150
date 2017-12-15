function [ETA MASSFLOW FIG] = CCGT3P(P_e,options,display)
% CCGT3P is a Combine cycle Gas Turbine with 2 pressure level
% CCGT3P(P_e,options,display) compute the thermodynamics states for a CCGT
% with 3 pressure level (cfr p166 english reference book) including
% combustion, exchanger and cycles. This is done based on several inputs
% (given in OPTION) and based on a given electricity production P_e.
% It returns the main results. It can as well plots graphs if input 
% argument DISPLAY = true (<=> DISPLAY=1)
%
% INPUTS (some inputs can be dependent on others => only one of these 2 can
%         be activated)
% P_EG = electrical power output target for gas turbine [W]
% OPTIONS is a structure containing :
%   -options.T0       [°C] : Reference temperature
%   -options.T_ext    [K]  : External temperature
%   -options.T_STmax  [°C] : maximum temperature on ST cycle
%   -options.eta_mec  [-] : mecanic efficiency of shafts bearings
%   -options.pdrum   [bar]: Drum pressure
%   -options.pmid    [bar]: Intermediary pressure level
%   -options.x7       [-] : Vapor ratio [gaseous/liquid] (titre)
%   -option.eta_SiC   [-] : Isotrenpic efficiency for compression
%   -option.eta_SiT   [-] : Isotrenpic efficiency for compression
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
%   -massflow(2) [kg/s]: water massflow at medium pressure turbine inlet
%   -massflow(3) [kg/s]: water massflow at low pressure turbine inlet
%   -massflow(4) [kg/s]: air massflow at gas turbine inlet 
%   -massflow(5) [kg/s]: combustible massflow
%% INPUT HANDLER
if nargin<3, display = 1; if nargin<2, options = struct(); if nargin<1, P_e = 153.8e6; end, end, end
if isfield(options,'T_0')    , T_0 = options.T_0;                else T_0 = 15;                 end %   -options.T0       [°C] : Reference temperature
if isfield(options,'T_ext')  , T_river = options.T_ext;          else T_river = 15;             end %   -options.T_ext    [K]  : External temperature
if isfield(options,'T_STMax'), T_STMax = options.T_STMax;        else T_STMax = 525;            end %   -options.T_STmax  [°C] : maximum temperature on ST cycle
if isfield(options,'eta_mec'), eta_mec = options.eta_mec;        else eta_mec = 0.9;            end %   -options.eta_mec  [-] : mecanic efficiency of shafts bearings
if isfield(options,'pdrum')  , pdrum = options.pdrum;            else pdrum =3.6;               end %   -options.pdrum   [bar]: Drum pressure
if isfield(options,'plow')   , plow = options.plow;              else plow =3.6;                end %   -options.plow    [bar]: Low pressure
if isfield(options,'pmid')   , pmid = options.pmid;              else pmid =27.3;               end %   -options.pmid    [bar]: Intermediary pressure level
if isfield(options,'phig')   , phig = options.phig;              else phig = 122.8;             end %   -options.phig    [bar]: Intermediary pressure level
if isfield(options,'x7')     , x7   = option.x7;                 else x7 = 1;                   end %   -options.x7       [-] : Vapor ratio [gaseous/liquid] (titre)
if isfield(options,'eta_SiC'), eta_SiC = options.eta_SiC;        else eta_SiC = 0.9;            end %   -option.eta_SiC   [-] : Isotrenpic efficiency for compression
if isfield(options,'eta_SiT'), eta_SiT = options.eta_SiT;        else eta_SiT(1) = 0.9;
                                                                      eta_SiT(2) = 0.9;  
                                                                      eta_SiT(3) = 0.9;         end %   -option.eta_SiT   [-] : Isotrenpic efficiency for compression
if isfield(options,'GT')     , optionsGT = options.GT;           else optionsGT.T3 = 1050+273.15;      end %-options.GT    [struct] : options for Gas turbine (see GT function) 
%% GTurbine
P_E_G = 283.7e6;
[GT_ETA GT_DATEN GT_DATEX GT_DAT GT_MASSFLOW GT_COMBUSTION] = GT(P_E_G,optionsGT,0);
mdot_g = GT_MASSFLOW(3);
Cp_g   = GT_COMBUSTION.Cp_g;
e_c    = GT_COMBUSTION.e_c;
%% INITIALISATION OF A FEW VALUE
T_pinch = 10;
%T_river = 15;
Cond_loss = 0.1;

h0_273 = XSteam('h_pT',1,T_0);
s0_273 = XSteam('s_pT',1,T_0);
%% AT THE ENTRANCE OF HP
% p T h s v x
table(3,1) = phig; %[bar] 
table(3,2) = min( T_STMax,GT_DAT(1,4)-T_pinch); %[C]
table(3,3) = XSteam('h_pT',table(3,1),table(3,2));
table(3,4) = XSteam('s_pT',table(3,1),table(3,2));
table(3,5) = XSteam('v_pT',table(3,1),table(3,2));
table(3,6) = table(3,3)-h0_273 - (T_0+273.15)*(table(3,4)-s0_273); 
%% AT THE EXIT OF THE TURBINE HP 
table(4,1) = pmid;
    h4s = XSteam('h_ps',table(4,1),table(3,4));
table(4,3) = eta_SiT(1)*h4s +  ( 1-eta_SiT(1) ) * table(3,3); 
table(4,2) = XSteam('T_ph',table(4,1),table(4,3));
table(4,4) = XSteam('s_ph',table(4,1),table(4,3));
table(4,5) = XSteam('v_ph',table(4,1),table(4,3));
table(4,6) = table(4,3)-h0_273 - (T_0+273.15)*(table(4,4)-s0_273); 
WmT_HP = table(4,3)-table(3,3);
%% AT THE ENTRANCE OF THE TURBINE IP
table(5,1) = pmid; %[bar] 
table(5,2) = min( T_STMax,GT_DAT(1,4)-T_pinch); %[C]
table(5,3) = XSteam('h_pT',table(5,1),table(5,2));
table(5,4) = XSteam('s_pT',table(5,1),table(5,2));
table(5,5) = XSteam('v_pT',table(5,1),table(5,2));
table(5,6) = table(5,3)-h0_273 - (T_0+273.15)*(table(5,4)-s0_273); 
%% AT THE EXIT OF THE TURBINE IP
table(6,1) = plow;
    h6s = XSteam('h_ps',table(6,1),table(5,4));
table(6,3) = eta_SiT(1)*h4s +  ( 1-eta_SiT(1) ) * table(5,3); 
table(6,2) = XSteam('T_ph',table(6,1),table(6,3));
table(6,4) = XSteam('s_ph',table(6,1),table(6,3));
table(6,5) = XSteam('v_ph',table(6,1),table(6,3));
table(6,6) = table(6,3)-h0_273 - (T_0+273.15)*(table(6,4)-s0_273); 
WmT_IP = table(6,3)-table(5,3);
%% AT THE EXIT OF THE TURBINE LP
table(7,2) = T_river+T_pinch+3;
table(7,1) = XSteam('psat_T',table(7,2));
table(7,7) = x7;
table(7,3) = XSteam('h_px',table(7,1),table(7,7));
table(7,4) = XSteam('s_ph',table(7,1),table(7,3));
table(7,5) = XSteam('v_ph',table(7,1),table(7,3));
table(7,6) = table(7,3)-h0_273 - (T_0+273.15)*(table(7,4)-s0_273); 
WmT_LP = table(7,3)-table(6,3);
%% AT THE EXIT OF THE CONDENSOR
table(1,1) = table(7,1)* (1-Cond_loss);
table(1,2) = T_river+T_pinch;
table(1,3) = XSteam('h_pT',table(1,1),table(1,2));
table(1,4) = XSteam('s_pT',table(1,1),table(1,2));
table(1,5) = XSteam('v_pT',table(1,2),table(1,2));
table(1,6) = table(1,3)-h0_273 - (T_0+273.15)*(table(1,4)-s0_273); 
%% AT THE EXIT OF THE MAIN PUMP
table(2,1) = pdrum; 
    h2s = XSteam('h_ps',table(2,1),table(1,4));
table(2,3) = 1/eta_SiC * h2s + (1-1/eta_SiC)*table(1,3);    
table(2,2) = XSteam('T_ph',table(2,1),table(2,3));
table(2,4) = XSteam('s_ph',table(2,1),table(2,3));
table(2,5) = XSteam('v_ph',table(2,1),table(2,3));
table(2,6) = table(2,3)-h0_273 - (T_0+273.15)*(table(2,4)-s0_273); 
WmC = table(2,3) - table(1,3);
%% Massflow
k_hp = 0.5;
k_ip = 0.4;
if abs(k_hp+k_ip) > 1, disp('TURBINE FRACTION WRONG') ,end  % 0 < k_hp+k_ip < 1   

Wm = k_hp*abs(WmT_HP)+(k_ip+k_hp)*abs(WmT_IP)+abs(WmT_LP) - WmC;
mdot_w = P_e/(Wm*10^3 /(eta_mec));
mdot_w_HP = k_hp*mdot_w; 
mdot_w_IP = k_ip*mdot_w;
mdot_w_LP = (1-k_hp-k_ip)*mdot_w;
%% THE CHIMNEY (OH NO ! NOT THE CHIMNEY !)
% 8 AFTER ECO LP
table(8,1) = table(2,1);
table(8,2) = XSteam('Tsat_p',pdrum);
table(8,3) = XSteam('hL_p',table(8,1));
table(8,4) = XSteam('sL_p',table(8,1));
table(8,5) = XSteam('vL_p',table(8,1));
table(8,6) = table(8,3)-h0_273 - (T_0+273.15)*(table(8,4)-s0_273); 
% 9 AFTER LP BALLON LIQUID + PUMP
table(9,1) = pmid;
    h9s = XSteam('h_ps',table(9,1),table(8,4));
table(9,3) = 1/eta_SiC * h9s + (1-1/eta_SiC)*table(8,3);
table(9,2) = XSteam('T_ph',table(9,1),table(9,3));
table(9,4) = XSteam('s_pt',table(9,1),table(9,2));
table(9,5) = XSteam('v_pt',table(9,1),table(9,2));
table(9,6) = table(9,3)-h0_273 - (T_0+273.15)*(table(9,4)-s0_273);
WmC2 = abs(table(9,3) - XSteam('HL_p',plow));
% 10 AFTER LP BALLON VAPOR
table(10,1) = plow;
table(10,2) = XSteam('Tsat_p',plow);
table(10,3) = XSteam('hV_p',table(10,1));
table(10,4) = XSteam('sV_p',table(10,1));
table(10,5) = XSteam('vV_p',table(10,1));
table(10,6) = table(10,3)-h0_273 - (T_0+273.15)*(table(10,4)-s0_273);
% 11 AFTER SUP LP
table(11,1) = plow;
table(11,2) = XSteam('Tsat_p',pmid);
table(11,3) = XSteam('h_pt',table(11,1),table(11,2));
table(11,4) = XSteam('s_pt',table(11,1),table(11,2));
table(11,5) = XSteam('v_pt',table(11,1),table(11,2));
table(11,6) = table(10,3)-h0_273 - (T_0+273.15)*(table(10,4)-s0_273);
% 12 AFTER ECO IP
table(12,1) = pmid;
table(12,2) = XSteam('Tsat_p',pmid);
table(12,3) = XSteam('hL_p',table(12,1));
table(12,4) = XSteam('sL_p',table(12,1));
table(12,5) = XSteam('vL_p',table(12,1));
table(12,6) = table(12,3)-h0_273 - (T_0+273.15)*(table(12,4)-s0_273);
% 13 AFTER IP BALLON LIQUID + Pump
table(13,1) = phig;
    h13s = XSteam('h_ps',table(13,1),table(12,4));
table(13,3) = 1/eta_SiC * h13s + (1-1/eta_SiC)*table(12,3);
table(13,2) = XSteam('T_ph',table(13,1),table(13,3));
table(13,4) = XSteam('s_pt',table(13,1),table(13,2));
table(13,5) = XSteam('v_pt',table(13,1),table(13,2));
table(13,6) = table(13,3)-h0_273 - (T_0+273.15)*(table(13,4)-s0_273);
WmC3 = abs(table(13,3) - XSteam('HL_p',plow));
% 14 AFTER ECO HP
table(14,1) = phig;
table(14,2) = XSteam('Tsat_p',phig);
table(14,3) = XSteam('hL_p',table(14,1));
table(14,4) = XSteam('sL_p',table(14,1));
table(14,5) = XSteam('vL_p',table(14,1));
table(14,6) = table(14,3)-h0_273 - (T_0+273.15)*(table(14,4)-s0_273);
% 15 AFTER IP BALLON Vapor
table(15,1) = pmid;
table(15,2) = XSteam('Tsat_p',pmid);
table(15,3) = XSteam('hV_p',table(15,1));
table(15,4) = XSteam('sV_p',table(15,1));
table(15,5) = XSteam('vV_p',table(15,1));
table(15,6) = table(15,3)-h0_273 - (T_0+273.15)*(table(15,4)-s0_273);
% 16 AFTER HP BALLOON VAPOR 
table(16,1) = phig;
table(16,2) = XSteam('Tsat_p',phig);
table(16,3) = XSteam('hV_p',table(16,1));
table(16,4) = XSteam('sV_p',table(16,1));
table(16,5) = XSteam('vV_p',table(16,1));
table(16,6) = table(16,3)-h0_273 - (T_0+273.15)*(table(16,4)-s0_273);
%% EXHAUST GAS FLOW
T_pinch_EVA_HP = 10;
T_pinch_EVA_IP = 10;
T_pinch_EVA_LP = 20;
table(17,1) = GT_DAT(2,4);
table(17,2) = GT_DAT(1,4);
table(17,3) = GT_DAT(3,4);
table(17,4) = GT_DAT(4,4);
table(17,6) = GT_DAT(5,4);
% AFTER SUP HP REHEAT IP
table(18,1) = table(17,1);
table(18,3) = table(17,3) - mdot_w_HP/mdot_g*(table(3,3)-table(16,3)) - (mdot_w_IP+mdot_w_LP)/mdot_g*(table(5,3)-table(4,3));
table(18,2) = table(17,2) + (table(18,3)-table(17,3))/Cp_g;
table(18,4) = Cp_g*log((table(18,2)+273.15)/(table(17,2)+273.15)) + table(17,4);
table(18,6) = table(18,3) - (T_0+273.15)*table(18,4); 
% AFTER EVA HP
table(19,1) = table(18,1);
table(19,2) = XSteam('Tsat_p',phig)+T_pinch_EVA_HP;
table(19,3) = table(18,3) + Cp_g*(table(19,2)-table(18,2));
table(19,4) = Cp_g*log((table(19,2)+273.15)/(table(18,2)+273.15)) + table(18,4);
table(19,6) = table(19,3) - (T_0+273.15)*table(19,4); 
% AFTER ECO HP
table(20,1) = table(19,1);
table(20,3) = table(19,3) - mdot_w_LP/mdot_g*(table(6,3)-table(11,3))... 
                          - mdot_w_IP/mdot_g*(table(5,3)-table(4,3))...
                          - mdot_w_HP/mdot_g*(table(14,3)-table(13,3));
table(20,2) = table(19,2) + (table(20,3)-table(19,3))/Cp_g;                      
table(20,4) = Cp_g*log((table(20,2)+273.15)/(table(19,2)+273.15)) + table(19,4);
table(20,6) = table(20,3) - (T_0+273.15)*table(20,4); 
% AFTER EVA IP
table(21,1) = table(20,1);
table(21,2) = XSteam('Tsat_p',pmid)+T_pinch_EVA_IP;
table(21,3) = table(20,3) + Cp_g*(table(21,2)-table(20,2));
table(21,4) = Cp_g*log((table(21,2)+273.15)/(table(20,2)+273.15)) + table(20,4);
table(21,6) = table(21,3) - (T_0+273.15)*table(21,4); 
% AFTER SUP LP ECO IP
table(22,1) = table(21,1);
table(22,3) = table(21,3) - mdot_w_LP/mdot_g*(table(11,3)-table(10,3))... 
                          - (mdot_w_HP+mdot_w_IP)/mdot_g*(table(12,3)-table(9,3));
table(22,2) = table(21,2) + (table(22,3)-table(21,3))/Cp_g;                      
table(22,4) = Cp_g*log((table(22,2)+273.15)/(table(21,2)+273.15)) + table(21,4);
table(22,6) = table(22,3) - (T_0+273.15)*table(22,4); 
% AFTER EVA LP
table(23,1) = table(22,1);
table(23,2) = XSteam('Tsat_p',plow)+T_pinch_EVA_LP;
table(23,3) = table(22,3) + Cp_g*(table(23,2)-table(22,2));
table(23,4) = Cp_g*log((table(23,2)+273.15)/(table(22,2)+273.15)) + table(22,4);
table(23,6) = table(23,3) - (T_0+273.15)*table(23,4); 
% AFTER ECO LP
table(24,1) = table(23,1);
table(24,3) = table(23,3) - mdot_w/mdot_g*(table(8,3)-table(2,3));
table(24,2) = table(23,2) + (table(24,3)-table(23,3))/Cp_g;                      
table(24,4) = Cp_g*log((table(24,2)+273.15)/(table(23,2)+273.15)) + table(23,4);
table(24,6) = table(24,3) - (T_0+273.15)*table(24,4); 
%% CYCLE HP AND LP MASSFLOW
h_EVA_HP = XSteam('HV_p',phig)-XSteam('HL_p',phig);
h_EVA_IP = XSteam('HV_p',pmid)-XSteam('HL_p',pmid);
h_EVA_LP = XSteam('HV_p',plow)-XSteam('HL_p',plow);
mdot_v_HP = (table(18,3)-table(19,3))/(h_EVA_HP) * mdot_g;
mdot_v_IP = (table(20,3)-table(21,3))/(h_EVA_IP) * mdot_g;
mdot_v_LP = (table(22,3)-table(23,3))/(h_EVA_LP) * mdot_g;
%% Output ETA
%OUPUTS : 
% ETA is a vector with :
%   -eta(1)  : eta_STcyclen, cycle energy efficiency
P_ST  = mdot_w*Wm;
ETA(1)= P_ST/(GT_COMBUSTION.LHV * GT_MASSFLOW(2));
%   -eta(2)  : eta_GTcyclen, cycle energy efficiency
ETA(2)= GT_ETA(1);
%   -eta(3)  : eta_toten, overall energy efficiency
Qech = mdot_g*GT_DAT(3,4) - GT_MASSFLOW(1)* GT_DAT(3,1);
epsilon_ech = Qech/(GT_COMBUSTION.LHV * GT_MASSFLOW(2));
ETA(3)= ETA(1) + ETA(2) - ETA(1)*ETA(2) - epsilon_ech*ETA(2);
%   -eta(4)  : eta_STcyclex, cycle exegy efficiency
ETA(4) = P_ST/( (mdot_w_HP)*(table(3,6)-table(2,6)) + (mdot_w_IP)*(table(4,6)-table(2,6)) + (mdot_w_LP)*(table(6,6)-table(2,6)) );
%   -eta(5)  : eta_GTcyclex, cycle exegy efficiency
ETA(5) = GT_ETA(3);
%   -eta(6)  : eta_totex, overall exergie efficiency
ETA(6) = (P_e+P_E_G)/1e3/(GT_MASSFLOW(2)*GT_COMBUSTION.e_c);
%   -eta(7)  : eta_gen, Steam generator energy efficiency
ETA(7) = ( mdot_w_HP*(table(3,3)-table(2,3)) + mdot_w_IP*(table(4,3)-table(2,3))+ mdot_w_LP*(table(6,3)-table(2,3)))/(GT_MASSFLOW(2)*GT_COMBUSTION.LHV);
%   -eta(8)  : eta_gex, Steam generator exergy efficiency
ETA(8) = ( mdot_w_HP*(table(3,6)-table(2,6)) + mdot_w_IP*(table(4,6)-table(2,6))+ mdot_w_LP*(table(6,6)-table(2,6)))/(GT_MASSFLOW(2)*e_c);
%   -eta(9)  : eta_combex, Combustion exergy efficiency
ETA(9) = GT_ETA(6);
%   -eta(10) : eta_chemex, Chimney exergy efficiency (losses)
ETA(10) = (GT_DAT(5,4)-table(24,6))/(GT_DAT(5,4));
%   -eta(11) : eta_transex, Heat exchanger overall exergy efficiency
ETA(11) = ( mdot_w_HP*(table(3,6)-table(2,6)) + mdot_w_IP*(table(4,6)-table(2,6))+ mdot_w_LP*(table(6,6)-table(2,6)))/(GT_MASSFLOW(3)*(GT_DAT(5,4)-table(24,6)));

%% Output Massflow
% MASSFLOW is a vector containing : 
%   -massflow(1) [kg/s]: massflow of steam at 3 (shortcut high pressure)
%   -massflow(2) [kg/s]: massflow of steam at 6 (IP)
%   -massflow(3) [kg/s]: massflow of steam at 7 (all steam flow)
%   -massflow(4) [kg/s]: air massflow at gas turbine inlet 
%   -massflow(5) [kg/s]: combustible massflow
MASSFLOW(1) = mdot_w_HP;
MASSFLOW(2) = mdot_w_IP;
MASSFLOW(3) = mdot_w_LP;
MASSFLOW(4) = GT_MASSFLOW(1);
MASSFLOW(5) = GT_MASSFLOW(2);

%% DISPLAY
if display == 1, visibility = 'on'; , else visibility = 'off'; , end
%% FIG 1 Cycle T-S Diagram
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
        if i <= 6 && i ~= 2, plot([table(i,4) table(i+1,4)], [table(i,2) table(i+1,2)],'b'), end
    end
    plot([table(7,4) table(1,4)], [table(7,2) table(1,2)],'b')
    
    plot([table(2,4) table(8,4)], [table(2,2) table(8,2)],'b')
    plot([table(8,4) table(10,4)], [table(8,2) table(10,2)],'r--')
        plot([table(10,4) table(11,4)], [table(10,2) table(11,2)],'r--')
        plot([table(11,4) table(6,4)], [table(11,2) table(6,2)],'r--')
    plot([table(8,4) table(9,4)], [table(8,2) table(9,2)],'g--')
        plot([table(9,4) table(12,4)], [table(9,2) table(12,2)],'g--')
        plot([table(12,4) table(13,4)], [table(12,2) table(13,2)],'g--')
        plot([table(13,4) table(15,4)], [table(13,2) table(15,2)],'g--')
        plot([table(15,4) table(4,4)], [table(15,2) table(4,2)],'g--')
    plot([table(12,4) table(14,4)], [table(12,2) table(14,2)],'y--') 
        plot([table(14,4) table(16,4)], [table(14,2) table(16,2)],'y--')
        plot([table(16,4) table(3,4)], [table(16,2) table(3,2)],'y--') 
    
    plot([table(17,4) table(18,4)], [table(17,2) table(18,2)],'k--')
    plot([table(18,4) table(19,4)], [table(18,2) table(19,2)],'k--') 
    plot([table(19,4) table(20,4)], [table(19,2) table(20,2)],'k--') 
    plot([table(20,4) table(21,4)], [table(20,2) table(21,2)],'k--') 
    plot([table(21,4) table(22,4)], [table(21,2) table(22,2)],'k--') 
    plot([table(22,4) table(23,4)], [table(22,2) table(23,2)],'k--') 
    plot([table(23,4) table(24,4)], [table(23,2) table(24,2)],'k--') 
    
    h = zeros(4, 1);
    h(1) = plot(NaN,NaN,'--r');
    h(2) = plot(NaN,NaN,'--g');
    h(3) = plot(NaN,NaN,'--y');
    h(4) = plot(NaN,NaN,'--k');
    legend(h,'HP Evolution Part','IP Part Evolution','LP Part Evolution','Exhaust Gas Evolution','Location','southeast');    
        
    grid on;
    grid minor;
    hold off;
    title('CCGT3P T-S Diagram ');
    xlabel('Entropy [kJ/C°]')
    ylabel('Temperature [C°]')
%% FIG 2 Cycle H-S Diagram
    FIG(2) = figure('visible',visibility);
    hold on;
     T=linspace(0,375,1000);  
     for i= 1 : length(T)
          S1(i)=XSteam('sL_T',T(i));
          S2(i)=XSteam('sV_T',T(i));
          H1(i)=XSteam('HL_T',T(i));
          H2(i)=XSteam('HV_T',T(i));
     end   
     S=[S1,S2];
     H=[H1,H2];
    plot(S,H,'b:');
    for i = 1 : 16
        plot(table(i,4),table(i,3),'k^');
        text(table(i,4)+0.1,table(i,3)+0.1,sprintf('%d',i));
        if i <= 6 && i ~= 2, plot([table(i,4) table(i+1,4)], [table(i,3) table(i+1,3)],'b'), end
    end
    plot([table(7,4) table(1,4)], [table(7,3) table(1,3)],'b')
    
    plot([table(2,4) table(8,4)], [table(2,3) table(8,3)],'b')
    plot([table(8,4) table(10,4)], [table(8,3) table(10,3)],'r--')
        plot([table(10,4) table(11,4)], [table(10,3) table(11,3)],'r--')
        plot([table(11,4) table(6,4)], [table(11,3) table(6,3)],'r--')
    plot([table(8,4) table(9,4)], [table(8,3) table(9,3)],'g--')
        plot([table(9,4) table(12,4)], [table(9,3) table(12,3)],'g--')
        plot([table(12,4) table(13,4)], [table(12,3) table(13,3)],'g--')
        plot([table(13,4) table(15,4)], [table(13,3) table(15,3)],'g--')
        plot([table(15,4) table(4,4)], [table(15,3) table(4,3)],'g--')
    plot([table(12,4) table(14,4)], [table(12,3) table(14,3)],'y--') 
        plot([table(14,4) table(16,4)], [table(14,3) table(16,3)],'y--')
        plot([table(16,4) table(3,4)], [table(16,3) table(3,3)],'y--') 
        
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'--r');
    h(2) = plot(NaN,NaN,'--g');
    h(3) = plot(NaN,NaN,'--y');
    legend(h,'HP Evolution Part','IP Part Evolution','LP Part Evolution','Location','southeast');    
        
    grid on;
    grid minor;
    hold off;
    title('CCGT3P H-S Diagram ');
    xlabel('Entropy [kJ/C°kg]')
    ylabel('Enthalpy [kJ/kg]')    
%%  FIG 3 - WATER REPARTITION  
    FIG(3)= figure('visible',visibility);
    F2labels = {'T-WMF';'HP-WMF';'IP-WMF';'LP-WMF';'HP-V-WMF';'IP-V-WMF';'LP-V-WMF'};
    bar(1:7,[mdot_w mdot_w_HP mdot_w_IP mdot_w_LP mdot_v_HP mdot_v_IP mdot_v_LP ]);
    set(gca,'xticklabel',F2labels)
    title('Water Massflow Repartion')
    colormap parula    
%%  FIG 4 - ENERGETIC POWER
    FIG(4)= figure('visible',visibility);
    LOSS_MECA   = (P_e)/1e3*(1-eta_mec);
    LOSS_COND   = abs((table(7,3)-table(1,3)))*mdot_w;
    LOSS_CHIM   = table(24,3)*GT_MASSFLOW(3);
    LOSS_MECAGT = GT_DATEN(1)/1e3;
    ENPO = [P_e/1e3 P_E_G/1e3 LOSS_MECA LOSS_MECAGT LOSS_COND LOSS_CHIM];
    F3labels = {'Effective Power ST';'Effective Power GT';'Mecanic Losses';'Mecanic Losses GT';'Condensor Losses';'Chimney Losses'};
    pie(ENPO,F3labels);
    STREN1 = sprintf('Effective Power ST %.1f [MW]' ,P_e/1e6);
    STREN2 = sprintf('Effective Power GT %.1f [MW]' ,P_E_G/1e6);
    STREN3 = sprintf('Mecanic Losses %.1f [MW]'     ,LOSS_MECA/1e3);
    STREN4 = sprintf('Mecanic Losses GT %.1f [MW]'  ,LOSS_MECAGT/1e3);
    STREN5 = sprintf('Condensor Losses %.1f [MW]'   ,LOSS_COND/1e3);
    STREN6 = sprintf('Chimney Losses %.1f [MW]'     ,LOSS_CHIM/1e3);
    legend(STREN1,STREN2,STREN3,STREN4,STREN5,STREN6,'Location','southeastoutside')
    title('Primary Energy Power')
    colormap summer;
%%  FIG 5 - EXERGY FLOW
    FIG(5)= figure('visible',visibility);
    LOSS_EX_COMB = GT_DATEX(3)/1e3;
    LOSS_EX_CHEM = table(24,6)*GT_MASSFLOW(3);
    LOSS_EX_COND = mdot_w*(table(7,6)-table(1,6));
    LOSS_EX_ROTO = (mdot_w * WmC + mdot_w_HP * WmC2 - mdot_w*(table(2,6)-table(1,6))...
                                                    + mdot_w_HP*(table(3,6)-table(4,6))...
                                                    + (mdot_w_HP+mdot_w_IP)*(table(6,6)-table(5,6))...
                                                    + mdot_w*(table(6,6)-table(7,6)));
    LOSS_EX_TRAN = GT_MASSFLOW(2)*e_c...
                  -P_e/1e3-P_E_G/1e3...
                  -LOSS_MECA-LOSS_MECAGT-LOSS_EX_CHEM-LOSS_EX_COND;
    
    EXFL = [P_e/1e3 P_E_G/1e3 LOSS_MECA LOSS_EX_COMB LOSS_EX_CHEM LOSS_EX_TRAN LOSS_EX_ROTO LOSS_EX_COND];
    F4labels = {'Effective Power ST';'Effective Power GT';'Mecanic Losses';'Combustion Losses';'Chimney Losses';'Transfert Losses';'Rotoric Losses';'Condensor Losses'};
    pie(EXFL,F4labels);
    STREX1 = sprintf('Effective Power ST %.1f [MW]' ,P_e/1e6);
    STREX2 = sprintf('Effective Power GT %.1f [MW]' ,P_E_G/1e6);
    STREX3 = sprintf('Mecanic Losses %.1f [MW]'     ,(LOSS_MECA+LOSS_MECAGT)/1e3);
    STREX4 = sprintf('Combustion Losses %.1f [MW]'  ,LOSS_EX_COMB/1e3);
    STREX5 = sprintf('Chimney Losses %.1f [MW]'     ,LOSS_EX_CHEM/1e3);
    STREX6 = sprintf('Transfer Losses %.1f [MW]'    ,LOSS_EX_TRAN/1e3);
    STREX7 = sprintf('Rotoric Losses %.1f [MW]'     ,LOSS_EX_ROTO/1e3);
    STREX8 = sprintf('Condensor Losses %.1f [MW]'   ,LOSS_EX_COND/1e3);
    
    legend(STREX1,STREX2,STREX3,STREX4,STREX5,STREX6,STREX7,STREX8,'Location','southeastoutside')
    title('Primary Exergy Flux')
    colormap winter;    

%% DUE POINT
    T_due_point  = XSteam('Tsat_p', GT_COMBUSTION.fumTG(4)/GT_MASSFLOW(3)); 
if display == 1
    
format short g;
    disp('        STATE       p[bar]         T[C]         h[KJ]     s[KJ/C]            v        e[KJ]');
    disp([linspace(1,length(table(:,1)),length(table(:,1)))' table(:,1) table(:,2) (table(:,3)) (table(:,4)) table(:,5) table(:,6) ]);
    disp('Chimney Due point');
    disp(T_due_point); 
end

end