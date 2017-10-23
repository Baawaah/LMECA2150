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
    data.Turb_comp=7;
    data.SG_ploss=0.1;
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
subplot(2,2,1);

hold on;
title('Rankine-Hirn Cycle')
xlabel('s[kJ/kgK]')
ylabel('T[K]')
plot(data.result(1).s,data.result(1).T,'.','MarkerSize',15)
text(data.result(1).s,data.result(1).T,'1')
plot(data.result(20).s,data.result(20).T,'.','MarkerSize',15)
text(data.result(20).s,data.result(20).T,'20')
plot(data.result(21).s,data.result(21).T,'.','MarkerSize',15)
text(data.result(21).s,data.result(21).T,'21')
plot(data.result(22).s,data.result(22).T,'.','MarkerSize',15)
text(data.result(22).s,data.result(22).T,'22')
plot(data.result(3).s,data.result(3).T,'.','MarkerSize',15)
text(data.result(3).s,data.result(3).T,'3')
plot(data.result(4).s,data.result(4).T,'.','MarkerSize',15)
text(data.result(4).s,data.result(4).T,'4')

plot([data.result(1).s data.result(20).s], [data.result(1).T data.result(20).T])
plot([data.result(20).s data.result(21).s], [data.result(20).T data.result(21).T] )
plot([data.result(21).s data.result(22).s], [data.result(21).T data.result(22).T] )
plot([data.result(22).s data.result(3).s], [data.result(22).T data.result(3).T] )
plot([data.result(3).s data.result(4).s], [data.result(3).T data.result(4).T] )
plot([data.result(4).s data.result(1).s], [data.result(4).T data.result(1).T] )
hold off;

subplot(2,2,2);
hold on;
title('Rankine-Hirn Cycle')
xlabel('v[m^3/kg]')
ylabel('p[Pa]')

data.result(1).v
data.result(1).p

plot(data.result(1).v,data.result(1).p,'.','MarkerSize',15)
text(data.result(1).v,data.result(1).p,'1')
plot(data.result(20).v,data.result(20).p,'.','MarkerSize',15)
text(data.result(20).v,data.result(20).p,'20')
plot(data.result(21).v,data.result(21).p,'.','MarkerSize',15)
text(data.result(21).v,data.result(21).p,'21')
plot(data.result(22).v,data.result(22).p,'.','MarkerSize',15)
text(data.result(22).v,data.result(22).p,'22')
plot(data.result(3).v,data.result(3).p,'.','MarkerSize',15)
text(data.result(3).v,data.result(3).p,'3')
plot(data.result(4).v,data.result(4).p,'.','MarkerSize',15)
text(data.result(4).v,data.result(4).p,'4')

plot([data.result(1).v data.result(20).v], [data.result(1).p data.result(20).p] )
plot([data.result(20).v data.result(21).v], [data.result(20).p data.result(21).p] )
plot([data.result(21).v data.result(22).v], [data.result(21).p data.result(22).p] )
plot([data.result(22).v data.result(3).v], [data.result(22).p data.result(3).p] )
plot([data.result(3).v data.result(4).v], [data.result(3).p data.result(4).p] )
plot([data.result(4).v data.result(1).v], [data.result(4).p data.result(1).p] )
hold off;
grid off;

subplot(2,2,3);
labels = {'pump','steam generator','turbine','condensor'};
sum_energ = [data.result(1).h data.result(20).h+data.result(21).h+data.result(22).h data.result(3).h data.result(4).h];
pie(sum_energ,labels);

subplot(2,2,4);
labels = {'pump','steam generator','turbine','condensor'};
sum_exerg = [data.result(1).ex data.result(20).ex+data.result(21).ex+data.result(22).ex data.result(3).ex data.result(4).ex];
pie(sum_exerg,labels);
end