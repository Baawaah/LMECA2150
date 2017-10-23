function [data] = ST_SteamGen(data)
data.v_eau = 1/1000; %(volume massique eau)
data.result(20).p = data.result(3).p*(1+data.SG_ploss);
data.result(20).h = data.result(1).h+data.v_eau*(data.result(20).p-data.result(1).p)*data.eta_pomp;
data.result(20).T = XSteam('T_ph',data.result(20).p,data.result(20).h);
data.result(20).s = XSteam('s_ph',data.result(20).p,data.result(20).h);
data.result(20).v = XSteam('v_pT',data.result(20).p,data.result(20).T);
data.result(20).x = XSteam('x_ph',data.result(20).p,data.result(20).h);
data.result(20).ex = exergy(data.result(20).h, data.h_ref,data.result(20).s,data.s_ref,data.T_ref);

data.result(21).p = data.result(3).p*(1+data.SG_ploss/2);
data.result(21).h = XSteam('hL_p', data.result(21).p);
data.result(21).T = XSteam('Tsat_p',data.result(21).p);
data.result(21).s = XSteam('sL_p',data.result(21).p);
data.result(21).v = XSteam('vL_p',data.result(21).p);
data.result(21).x = 0;
data.result(21).ex = exergy(data.result(21).h, data.h_ref, data.result(21).s, data.s_ref, data.T_ref);

data.result(22).p = data.result(3).p;
data.result(22).h = XSteam('hV_p',data.result(22).p);
data.result(22).T = XSteam('Tsat_p',data.result(22).p);
data.result(22).s = XSteam('sV_p',data.result(22).p);
data.result(22).v = XSteam('vV_p',data.result(22).p);
data.result(22).x = 1;
data.result(22).ex = exergy(data.result(22).h, data.h_ref, data.result(22).s, data.s_ref, data.T_ref);
end

