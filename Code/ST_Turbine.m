function [data] = ST_Turbine(data)
data.result(4).p = data.result(3).p;
s_4_is = data.result(3).s;
h_4_is = XSteam('h_ps', data.result(4).p, s_4_is);
data.result(4).h = data.result(3).h + data.eta_is_turb*(h_4_is-data.result(3).h);
data.result(4).T = XSteam('T_ph', data.result(4).p, data.result(4).h);
data.result(4).x = XSteam('x_ph', data.result(4).p, data.result(4).h);
data.result(4).s = XSteam('s_ph', data.result(4).p, data.result(4).h);
data.result(4).ex = exergy(data.result(4).h, data.h_ref, data.result(4).s, data.s_ref, data.T_ref);
end

