function [data] = ST_TurbIP(data)
data.result(32).p = data.result(31).p/data.TurbIP_comp;
s_32_is = data.result(31).s;
h_32_is = XSteam('h_ps', data.result(32).p, s_32_is);
data.result(32).h = data.result(31).h + data.eta_is_turb*(h_32_is-data.result(31).h);
data.result(32).T = XSteam('T_ph', data.result(32).p, data.result(32).h);
data.result(32).x = XSteam('x_ph', data.result(32).p, data.result(32).h);
data.result(32).s = XSteam('s_ph', data.result(32).p, data.result(32).h);
data.result(32).v = XSteam('v_ph',data.result(32).p,data.result(32).h);
data.result(32).ex = exergy(data.result(32).h, data.h_ref, data.result(32).s, data.s_ref, data.T_ref);
end
