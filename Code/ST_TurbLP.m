function [data] = ST_TurbLP(data)
data.result(4).p = data.result(32).p/data.TurbLP_comp;
s_4_is = data.result(32).s;
h_4_is = XSteam('h_ps', data.result(32).p, s_4_is);
%data.result(4).T = data.T_river + data.Cond_pinch;
%data.result(4).p = XSteam('Psat_T',data.T_river + data.Cond_pinch);
data.result(4).h = data.result(32).h + data.eta_is_turb*(h_4_is-data.result(32).h);
data.result(4).T = XSteam('T_ph',data.result(4).p,data.result(4).h);
data.result(4).x = XSteam('x_ph', data.result(4).p, data.result(4).h);
data.result(4).s = XSteam('s_ph', data.result(4).p, data.result(4).h);
data.result(4).v = XSteam('v_ph',data.result(4).p,data.result(4).h);
data.result(4).ex = exergy(data.result(4).h, data.h_ref, data.result(4).s, data.s_ref, data.T_ref);
data.result(32).p
data.result(4).p
end
