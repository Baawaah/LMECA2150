function [ data_output ] = rankin(data)
%% =============================
% Pump
 data.Table(2,1) = data.Table(1,1) * data.rankin_PA_ratio; 
 data.Table(2,4) = data.Table(1,4)+ (data.Table(2,1)-data.Table(1,1))*data.Table(1,3)*10^2;
 data.Table(2,2) = XSteam('T_ph',data.Table(2,1),data.Table(2,4));
 data.Table(2,3) = XSteam('vL_T',data.Table(2,2));
 data.Table(2,5) = XSteam('s_ph',data.Table(2,1),data.Table(2,4));
 
% Connection P-SG 
 data.Table(3,:) = data.Table(2,:);
 
% Steam Generator
 data.Table(4,1) = data.Table(3,1)*(1-data.rankin_SG_pdrop);
 data.Table(4,2) = data.rankin_SG_T_out;
 data.Table(4,3) = XSteam('v_pT',data.Table(4,1),data.Table(4,2));
 data.Table(4,4) = data.Table(3,4)+ XSteam('h_pT',data.Table(4,1),data.Table(4,2)) - XSteam('h_pT',data.Table(3,1),data.Table(3,2)); 
 data.Table(4,5) = XSteam('s_ph',data.Table(4,1),data.Table(4,4));
 
% Connection SG-TU
 data.Table(5,:) = data.Table(4,:);
 
% Turbine 
 data.Table(6,1) = data.Table(5,1) * data.rankin_TU_ratio;
 data.Table(6,2) = XSteam('Tsat_p',data.Table(6,1));
 data.Table(6,3) = XSteam('vV_p',data.Table(6,1));
 data.Table(6,4) = data.Table(5,4) + XSteam('hV_p',data.Table(6,1)) - XSteam('h_pT',data.Table(5,1),data.Table(5,2));
 data.Table(6,5) = XSteam('s_ph',data.Table(6,1),data.Table(6,4));
 
% Connection TU-CO 
 data.Table(7,:) = data.Table(6,:);
  
% Condensor
 data.Table(8,1) = XSteam('psat_T',data.Table(7,2))*(1-data.rankin_CO_pdrop);
 data.Table(8,2) = XSteam('Tsat_p',data.Table(8,1));
 data.Table(8,3) = XSteam('vL_p',data.Table(8,1));
 data.Table(8,4) = data.Table(7,4) + XSteam('hL_T',data.Table(8,2)) - XSteam('hV_T',data.Table(7,2));
 data.Table(8,5) = XSteam('s_ph',data.Table(8,1),data.Table(8,4));
 
 data_output = data;

end

