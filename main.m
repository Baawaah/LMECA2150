%% =============================
% 
%
%
%% =============================
%Ratkins Cycle
%
% Parameter
PA_ratio =  1.0353*10^3;  % Compresion Ratio of the Pump% P2/P1
SG_T_out = 525         ;  % Maximum Temperature out of the Steam Genetator
TU_ratio = 0.0012575   ;  % Expansion Ratio of the Turbine
SG_pdrop = 0.1         ;  % Pressure drop in the Steam Generator
CO_pdrop = 0.14653     ;  % Pressure drop in the condensor
%% =============================
% Master Table
% Je réserve:
% 1,2 Pompe Alimentaire In-Out
% 3,4 Steam Generator   In-Out
% 4,5 Turbine           In-Out
% 5,6 Condensor         In-Out

%        p T v h
Table = [];

 %Table(1,:) = [0.0425 30 0 0];
 Table(1,:) = [0.05 0 0 0]; 
 Table(1,2) = XSteam('Tsat_p',Table(1,1));
 Table(1,3) = XSteam('vL_T',Table(1,2));
 Table(1,4) = XSteam('hL_T',Table(1,2));
%% =============================
% Pump
 Table(2,1) = Table(1,1) * PA_ratio; 
 Table(2,4) = Table(1,4)+ (Table(2,1)-Table(1,1))*Table(1,3)*10^2;
 Table(2,2) = XSteam('T_pH',Table(2,1),Table(2,4));
 Table(2,3) = XSteam('vL_T',Table(2,2));
% Connection P-SG 
 Table(3,:) = Table(2,:);
 
% Steam Generator
 Table(4,1) = Table(3,1)*(1-SG_pdrop);
 Table(4,2) = SG_T_out;
 Table(4,3) = XSteam('v_pT',Table(4,1),Table(4,2));
 Table(4,4) = Table(3,4)+ XSteam('h_pT',Table(4,1),Table(4,2)) - XSteam('h_pT',Table(3,1),Table(3,2)); 
 
% Connection SG-TU
 Table(5,:) = Table(4,:);
 
% Turbine 
 Table(6,1) = Table(5,1) * TU_ratio;
 Table(6,2) = XSteam('Tsat_p',Table(6,1));
 Table(6,3) = XSteam('vV_p',Table(6,1));
 Table(6,4) = Table(5,4) + XSteam('hV_p',Table(6,1)) - XSteam('h_pT',Table(5,1),Table(5,2));
% Connection TU-CO 
 Table(7,:) = Table(6,:);
  
% Condensor
 Table(8,1) = XSteam('psat_T',Table(7,2))*(1-CO_pdrop);
 Table(8,2) = XSteam('Tsat_p',Table(8,1));
 Table(8,3) = XSteam('vL_p',Table(8,1));
 Table(8,4) = Table(7,4) + XSteam('hL_T',Table(8,2)) - XSteam('hV_T',Table(7,2));

 

 Table
