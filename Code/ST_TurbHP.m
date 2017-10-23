function [data] = ST_TurbHP(data)
data.result(31).p = data.result(3).p/data.TurbHP_comp;
s_31_is = data.result(3).s;
h_31_is = XSteam('h_ps', data.result(31).p, s_31_is);
data.result(31).h = data.result(3).h + data.eta_is_turb*(h_31_is-data.result(3).h);
data.result(31).T = XSteam('T_ph', data.result(31).p, data.result(31).h);
data.result(31).x = XSteam('x_ph', data.result(31).p, data.result(31).h);
data.result(31).s = XSteam('s_ph', data.result(31).p, data.result(31).h);
data.result(31).v = XSteam('v_pT',data.result(31).p,data.result(31).T);
data.result(31).ex = exergy(data.result(31).h, data.h_ref, data.result(31).s, data.s_ref, data.T_ref);
end

