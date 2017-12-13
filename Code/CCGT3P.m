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
if isfield(options,'T_ext')  , optionsGT.T_ext = options.T_ext;  else optionsGT.T_ext = 288.15; end %   -options.T_ext    [K]  : External temperature
if isfield(options,'T_STMax'), T_STMax = options.T_STMax;        else T_STMax = 525;            end %   -options.T_STmax  [°C] : maximum temperature on ST cycle
if isfield(options,'eta_mec'), eta_mec = options.eta_mec;        else eta_mec = 0.9;            end %   -options.eta_mec  [-] : mecanic efficiency of shafts bearings
if isfield(options,'pdrum')  , pdrum = options.pdrum;            else pdrum =3.6;               end %   -options.pdrum   [bar]: Drum pressure
if isfield(options,'plow')   , plow = options.plow;              else plow =3.6;                end %   -options.plow    [bar]: Low pressure
if isfield(options,'pmid')   , pmid = options.pmid;              else pmid =27.3;               end %   -options.pmid    [bar]: Intermediary pressure level
if isfield(options,'phig')   , phig = options.phig;              else phig = 122.8;             end %   -options.phig    [bar]: Intermediary pressure level
if isfield(options,'x7')     , x7   = option.x7;                 else x7 = 1;                   end %   -options.x7       [-] : Vapor ratio [gaseous/liquid] (titre)
if isfield(options,'eta_SiC'), eta_SiC = options.eta_SiC;        else eta_SiC = 0.9;            end %   -option.eta_SiC   [-] : Isotrenpic efficiency for compression
if isfield(options,'eta_SiT'), eta_SiT = options.eta_SiT;        else eta_SiT(1) = 0.9;
                                                                      eta_SiT(2) = 0.9;         end %   -option.eta_SiT   [-] : Isotrenpic efficiency for compression
if isfield(options,'GT')     , optionsGT = options.GT;           else optionsGT.comb.Tmax = 1150;
                                                                      optionsGT.comb.x = 1; 
                                                                      optionsGT.comb.y = 4;
                                                                      optionsGT.comb.lambda = 3;end %-options.GT    [struct] : options for Gas turbine (see GT function) 
%% GTurbine
P_E_G = 283.7e6;
[GT_ETA GT_DATEN GT_DATEX GT_DAT GT_MASSFLOW GT_COMBUSTION] = GT(P_E_G,optionsGT,0);
%% INITIALISATION OF A FEW VALUE
T_pinch = 10;
T_river = 15;
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
table(2,1) = plow; 
    h2s = XSteam('h_ps',table(2,1),table(1,4));
table(2,3) = 1/eta_SiC * h2s + (1-1/eta_SiC)*table(1,3);    
table(2,2) = XSteam('T_ph',table(2,1),table(2,3));
table(2,4) = XSteam('s_ph',table(2,1),table(2,3));
table(2,5) = XSteam('v_ph',table(2,1),table(2,3));
table(2,6) = table(2,3)-h0_273 - (T_0+273.15)*(table(2,4)-s0_273); 
WmC = table(2,3) - table(1,3);
%% Mass flow
k_hp = 0.8;
k_ip = 0.1;
if abs(k_hp+k_ip) > 1
    disp('TURBINE FRACTION WRONG')
end    
% 0 < k_hp+k_ip < 1
Wm = k_hp*abs(WmT_HP)+(k_ip+k_hp)*abs(WmT_IP)+abs(WmT_LP) - WmC;
mdot_w = P_e/(Wm*10^3 /(eta_mec));
mdot_w_HP = k_hp*mdot_w; % First Guess
mdot_w_IP = k_ip*mdot_w;
mdot_w_LP = (1-k_hp-k_ip)*mdot_w;

%% COMBUSTION
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
    plot(S,T);
    for i = 1 : length(table(:,1))
        plot(table(i,4),table(i,2),'b^');
        text(table(i,4)+0.1,table(i,2)+0.1,sprintf('%d',i));
        if i <= length(table(:,1))-1, plot([table(i,4) table(i+1,4)], [table(i,2) table(i+1,2)]), end
    end
    plot([table(length(table(:,1)),4) table(1,4)], [table(length(table(:,1)),2) table(1,2)])
    grid on;
    grid minor;
    hold off;
    title('CCGT3P T-S Diagram ');
    xlabel('Entropy [kJ/C°]')
    ylabel('Temperature [C°]')
%%    
%     format short g;
%     disp('       p[bar]         T[C]         h[MJ]     s[KJ/C]            v        e[KJ]');
%     disp([table(:,1) table(:,2) (table(:,3)/10^3) (table(:,4)) table(:,5) table(:,6) ]);
%     
% %     format short;
% %     disp('Tgas[C]       hgas       p[bar]     T[C]      h[MJ]     s[J/K]        v      x[KJ]');
% %     disp([table_chim(:,1) table_chim(:,2) table_chim(:,3) table_chim(:,4) table_chim(:,5)/10^3  ]);
%     
%     disp('  WmC[KJ]    WmT[KJ]    Wm[KJ]    mdot_g    mdot_w    mdot_w_HP  ');
%     disp([WmT_HP WmT_LP  Wm mdot_g mdot_w mdot_w_HP]);

end