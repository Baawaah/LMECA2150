function [data] = ST_Condensor(data)
%% Sortie de condenseur (liquide saturé : x=0) :
data.result(1).T = data.T_river+data.deltaT_river;
data.result(1).p = XSteam('psat_T',data.result(1).T)
data.result(1).h = XSteam('hL_p', data.result(1).p);
data.result(1).s = XSteam('sL_p', data.result(1).p);
data.result(1).x = 0;
data.result(1).ex = exergy(data.result(1).h, data.h_ref, data.result(1).s, data.s_ref, data.T_ref);
end

