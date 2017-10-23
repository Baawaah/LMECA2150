function [data] = ST_Init(data)
data.result(3).p = data.p_turb;
data.result(3).T = data.T_turb;
data.result(3).h = XSteam('h_pT',data.result(3).p,data.result(3).T);
data.result(3).s = XSteam('s_pT',data.result(3).p,data.result(3).T);
data.result(3).v = XSteam('v_pT',data.result(3).p,data.result(3).T);
data.result(3).x = XSteam('x_ph',data.result(3).p,data.result(3).h);
data.result(3).ex = exergy(data.result(3).h,data.h_ref,data.result(3).s,data.s_ref,data.T_ref);
end

