function [ data ] = init_data()
%% =============================
% Master Table
% Je réserve:
% 1,2 Pompe Alimentaire In-Out 
% 3,4 Steam Generator   In-Out 
% 5,6 Turbine           In-Out 
% 7,8 Condensor         In-Out
%                     p   T v h s
 data.Table(1,:) = [0.050 0 0 0 0]; 
 data.Table(1,2) = XSteam('Tsat_p',data.Table(1,1));
 data.Table(1,3) = XSteam('vL_T',data.Table(1,2));
 data.Table(1,4) = XSteam('hL_T',data.Table(1,2));
 data.Table(1,5) = XSteam('s_ph',data.Table(1,1),data.Table(1,4));

% Rankin
data.rankin_PA_ratio =  1.0353*10^3;  % Compresion Ratio of the Pump% P2/P1
data.rankin_SG_T_out = 525         ;  % Maximum Temperature out of the Steam Genetator
data.rankin_TU_ratio = 0.0012575   ;  % Expansion Ratio of the Turbine
data.rankin_SG_pdrop = 0.1         ;  % Pressure drop in the Steam Generator
data.rankin_CO_pdrop = 0.14653     ;  % Pressure drop in the condensor

% Soutirage
% data.Table
%     .Sout_Table
%     .Sout_deb
%     .Sout_pitch
%     .Sout_n
%     .Sout_p_extract
data.Sout_pitch     = 5;
data.Sout_n         = 5;
data.Sout_extract = [ 0.2 0.3 0.5 0.7 0.8 ];

end

