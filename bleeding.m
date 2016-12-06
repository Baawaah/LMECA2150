function [data_output] = bleeding(data)
%%  
% TB_sout est le tableau des états  du feedwater se trouvant dans le 
% soutirage qui sera construit en parallele avec la matrice de résolution 
% des débits
% Turbine Extraction data.Sout_n
p_extract = data.Sout_extract * data.Table(5,1);
data.Sout_Table = zeros(data.Sout_n*4+4,5);
% Connection
data.Sout_Table(1,:) = data.Table(8,:);
% Entropy
s_i = @(i) data.Table(5,5) + (data.Table(6,5)-data.Table(5,5))*(data.Sout_n + 1 - i)/(data.Sout_n + 1);
% Remplissage de TB_sout
% p T v h s
% Chaque niveau possède 2 pts waterfeed and bleeding
% Reheat Block
reheat_sout = 10;
% Just pour chauffer un peu l'eau du soutirage ici 5 degree
%  Input Bleeding
data.Sout_Table(3,1) = p_extract(1);
data.Sout_Table(3,2) = XSteam('Tsat_p',data.Sout_Table(3,1));
data.Sout_Table(3,4) = XSteam('hL_p',data.Sout_Table(3,1));
% Output Bleeding
data.Sout_Table(4,1) = data.Sout_Table(3,1);
data.Sout_Table(4,2) = data.Table(8,2) + reheat_sout;
data.Sout_Table(4,4) = XSteam('h_pT',data.Sout_Table(4,1),data.Sout_Table(4,2));
%  Output Feedwater
data.Sout_Table(2,1) = data.Sout_Table(1,1);
data.Sout_Table(2,2) = data.Sout_Table(4,2)+5;
data.Sout_Table(2,4) = XSteam('h_pT',data.Sout_Table(2,1),data.Sout_Table(2,2));

for i = 1 : data.Sout_n
    %  Input Bleeding
    data.Sout_Table(i*4+3,1) = p_extract(i);
    data.Sout_Table(i*4+3,5) = s_i(i);
    data.Sout_Table(i*4+3,2) = XSteam('T_ps',data.Sout_Table(i*4+3,1),data.Sout_Table(i*4+3,5));
    data.Sout_Table(i*4+3,4) = XSteam('h_ps',data.Sout_Table(i*4+3,1),data.Sout_Table(i*4+3,5));
    % Output Bleeding
    data.Sout_Table(i*4+4,1) = data.Sout_Table(i*4+3,1); 
    data.Sout_Table(i*4+4,2) = XSteam('Tsat_p',data.Sout_Table(i*4+4,1));
    data.Sout_Table(i*4+4,4) = XSteam('hL_T',data.Sout_Table(i*4+4,2));
    data.Sout_Table(i*4+4,5) = XSteam('sL_p',data.Sout_Table(i*4+4,1));
    %  Input Feedwater
    data.Sout_Table(i*4+1,1) = data.Sout_Table(1,1);
    data.Sout_Table(i*4+1,2) = data.Sout_Table(i*4-2,2);
    data.Sout_Table(i*4+1,4) = XSteam('h_pT',data.Sout_Table(i*4+1,1),data.Sout_Table(i*4+1,2));
    % Output Feedwater
    data.Sout_Table(i*4+2,1) = data.Sout_Table(i*4+1,1);
    data.Sout_Table(i*4+2,2) = data.Sout_Table(i*4+4,2) + data.Sout_pitch; 
    data.Sout_Table(i*4+2,4) = XSteam('h_pT',data.Sout_Table(i*4+2,1),data.Sout_Table(i*4+2,2));
    
end
% Calcule des debits
TB_sout_A = zeros(data.Sout_n,data.Sout_n);
TB_sout_B = zeros(data.Sout_n, 1);
for i = 1 : data.Sout_n
    for j = 1 : data.Sout_n
        H_tfo = data.Sout_Table(i*4+2,4);
        H_tfi = data.Sout_Table(i*4+1,4);
        TB_sout_A(i,j) = ( H_tfo - H_tfi ) + TB_sout_A(i,j);
    end
    H_tbo = data.Sout_Table(i*4+4,4);
    H_tbi = data.Sout_Table(i*4+3,4);
    TB_sout_A(i,i) = TB_sout_A(i,i) - (H_tbi - H_tbo); 
    for j = i+1 : data.Sout_n
        TB_sout_A(i,j) = TB_sout_A(i,j) - (H_tbo - H_tbi); 
    end
    TB_sout_B(i) = H_tfo - H_tfi;
end
% on met à jour le tableau
data.Sout_deb = (-1)*(TB_sout_A\TB_sout_B);
data_output = data;
end

