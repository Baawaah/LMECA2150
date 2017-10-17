function [result]=Rankine_cycle(T_river, deltaT_river, p_turb, T_turb, eta_is_turb, eta_pomp, eta_meca, P_eff)
clc
close all
if nargin==0
    T_river=15+273.15;
    deltaT_river=15;
    p_turb=200;
    T_turb=525;
    eta_is_turb=0.88;
    eta_pomp=1;
    eta_meca=0.9;
    P_eff=400;
    
    data.T_river=15+273.15;
    data.deltaT_river=15;
    data.p_turb=200;
    data.T_turb=525;
    data.eta_is_turb=0.88;
    data.eta_pomp=1;
    data.eta_meca=0.9;
    data.P_eff=400;
end

%% etat de reference
data.T_ref=25;
data.p_ref=1;
data.h_ref=XSteam('h_pT',data.p_ref,data.T_ref);
data.s_ref=XSteam('s_pT',data.p_ref,data.T_ref);

%% Entrée turbine :
% result(3).p = p_turb;
% result(3).T = T_turb;
% result(3).h = XSteam('h_pT',result(3).p,result(3).T);
% result(3).s = XSteam('s_pT',result(3).p,result(3).T);
% result(3).x = XSteam('x_ph',result(3).p,result(3).h);
% result(3).ex = exergy(result(3).h, h_ref, result(3).s, s_ref, T_ref);
data = ST_Init(data);
data = ST_Turbine(data);
data = ST_Condensor(data);
data = ST_SteamGen(data);



%% Plot 
figure
hold on;
title('Rankine-Hirn Cycle')
xlabel('s[kJ/kgK]')
ylabel('T[K]')
plot(data.result(1).s,data.result(1).T,'.','MarkerSize',15)
plot(data.result(20).s,data.result(20).T,'.','MarkerSize',15)
plot(data.result(21).s,data.result(21).T,'.','MarkerSize',15)
plot(data.result(22).s,data.result(22).T,'.','MarkerSize',15)
plot(data.result(3).s,data.result(3).T,'.','MarkerSize',15)
plot(data.result(4).s,data.result(4).T,'.','MarkerSize',15)

grid off

end