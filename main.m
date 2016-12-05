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
% 5,6 Turbine           In-Out 
% 7,8 Condensor         In-Out

%        p T v h s
Table = [];

 %Table(1,:) = [0.0425 30 0 0 0];
 Table(1,:) = [0.05 0 0 0 0]; 
 Table(1,2) = XSteam('Tsat_p',Table(1,1));
 Table(1,3) = XSteam('vL_T',Table(1,2));
 Table(1,4) = XSteam('hL_T',Table(1,2));
 Table(1,5) = XSteam('s_ph',Table(1,1),Table(1,4));
%% =============================
% Pump
 Table(2,1) = Table(1,1) * PA_ratio; 
 Table(2,4) = Table(1,4)+ (Table(2,1)-Table(1,1))*Table(1,3)*10^2;
 Table(2,2) = XSteam('T_ph',Table(2,1),Table(2,4));
 Table(2,3) = XSteam('vL_T',Table(2,2));
 Table(2,5) = XSteam('s_ph',Table(2,1),Table(2,4));
 
% Connection P-SG 
 Table(3,:) = Table(2,:);
 
% Steam Generator
 Table(4,1) = Table(3,1)*(1-SG_pdrop);
 Table(4,2) = SG_T_out;
 Table(4,3) = XSteam('v_pT',Table(4,1),Table(4,2));
 Table(4,4) = Table(3,4)+ XSteam('h_pT',Table(4,1),Table(4,2)) - XSteam('h_pT',Table(3,1),Table(3,2)); 
 Table(4,5) = XSteam('s_ph',Table(4,1),Table(4,4));
 
% Connection SG-TU
 Table(5,:) = Table(4,:);
 
% Turbine 
 Table(6,1) = Table(5,1) * TU_ratio;
 Table(6,2) = XSteam('Tsat_p',Table(6,1));
 Table(6,3) = XSteam('vV_p',Table(6,1));
 Table(6,4) = Table(5,4) + XSteam('hV_p',Table(6,1)) - XSteam('h_pT',Table(5,1),Table(5,2));
 Table(6,5) = XSteam('s_ph',Table(6,1),Table(6,4));
 
% Connection TU-CO 
 Table(7,:) = Table(6,:);
  
% Condensor
 Table(8,1) = XSteam('psat_T',Table(7,2))*(1-CO_pdrop);
 Table(8,2) = XSteam('Tsat_p',Table(8,1));
 Table(8,3) = XSteam('vL_p',Table(8,1));
 Table(8,4) = Table(7,4) + XSteam('hL_T',Table(8,2)) - XSteam('hV_T',Table(7,2));
 Table(8,5) = XSteam('s_ph',Table(8,1),Table(8,4));

 

 Table
 
%%  
% Soutirage

% Approximation T(p) = To-Ti/ Po-Pi * p + (Ti*Po - T0*Pi)/(Po-Pi)

TB_sout     = zeros(N_soutirage*4+4,5);
pitch = 5;
% TB_sout est le tableau des états  du feedwater se trouvant dans le 
% soutirage qui sera construit en parallele avec la matrice de résolution 
% des débits

% Connection
TB_sout(1,:) = Table(8,:);

% Turbine Extraction
p_extract = [ 0.2 0.3 0.5 0.7 0.8 ] * Table(5,1);
N_soutirage = length(p_extract);
% Entropy
s_i = @(i) Table(5,5) + (Table(6,5)-Table(5,5))*(N_soutirage + 1 - i)/(N_soutirage + 1);
% Remplissage de TB_sout
% p T v h s
% Chaque niveau possède 2 pts waterfeed and bleeding
% Reheat Block
reheat_sout = 10;
% Just pour chauffer un peu l'eau du soutirage ici 5 degree
%  Input Bleeding
TB_sout(3,1) = p_extract(1);
TB_sout(3,2) = XSteam('Tsat_p',TB_sout(3,1));
TB_sout(3,4) = XSteam('hL_p',TB_sout(3,1));
% Output Bleeding
TB_sout(4,1) = TB_sout(3,1);
TB_sout(4,2) = Table(8,2) + reheat_sout;
TB_sout(4,4) = XSteam('h_pT',TB_sout(4,1),TB_sout(4,2));
%  Output Feedwater
TB_sout(2,1) = TB_sout(1,1);
TB_sout(2,2) = TB_sout(4,2)+5;
TB_sout(2,4) = XSteam('h_pT',TB_sout(2,1),TB_sout(2,2));

for i = 1 : N_soutirage
    %  Input Bleeding
    TB_sout(i*4+3,1) = p_extract(i);
    TB_sout(i*4+3,5) = s_i(i);
    TB_sout(i*4+3,2) = XSteam('T_ps',TB_sout(i*4+3,1),TB_sout(i*4+3,5));
    TB_sout(i*4+3,4) = XSteam('h_ps',TB_sout(i*4+3,1),TB_sout(i*4+3,5));
    % Output Bleeding
    TB_sout(i*4+4,1) = TB_sout(i*4+3,1); 
    TB_sout(i*4+4,2) = XSteam('Tsat_p',TB_sout(i*4+4,1));
    TB_sout(i*4+4,4) = XSteam('hL_T',TB_sout(i*4+4,2));
    TB_sout(i*4+4,5) = XSteam('sL_p',TB_sout(i*4+4,1));
    %  Input Feedwater
    TB_sout(i*4+1,1) = TB_sout(1,1);
    TB_sout(i*4+1,2) = TB_sout(i*4-2,2);
    TB_sout(i*4+1,4) = XSteam('h_pT',TB_sout(i*4+1,1),TB_sout(i*4+1,2));
    % Output Feedwater
    TB_sout(i*4+2,1) = TB_sout(i*4+1,1);
    TB_sout(i*4+2,2) = TB_sout(i*4+4,2) + pitch; 
    TB_sout(i*4+2,4) = XSteam('h_pT',TB_sout(i*4+2,1),TB_sout(i*4+2,2));
    
end
% Calcule des debits
TB_sout_A = zeros(N_soutirage, N_soutirage);
TB_sout_B = zeros(N_soutirage, 1);
for i = 1 : N_soutirage
    for j = 1 : N_soutirage
        H_tfo = TB_sout(i*4+2,4);
        H_tfi = TB_sout(i*4+1,4);
        TB_sout_A(i,j) = ( H_tfo - H_tfi ) + TB_sout_A(i,j);
    end
    H_tbo = TB_sout(i*4+4,4);
    H_tbi = TB_sout(i*4+3,4);
    TB_sout_A(i,i) = TB_sout_A(i,i) - (H_tbi - H_tbo); 
    for j = i+1 : N_soutirage
        TB_sout_A(i,j) = TB_sout_A(i,j) - (H_tbo - H_tbi); 
    end
    TB_sout_B(i) = H_tfo - H_tfi;
end
% Vecteur des débits
Debit_sout = (-1)*(TB_sout_A\TB_sout_B);
% on met à jour le tableau



    
 