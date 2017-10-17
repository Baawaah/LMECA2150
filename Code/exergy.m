function [ex] = exergy(h,h_ref,s,s_ref,T_ref)
        ex = (h-h_ref+T_ref*(s-s_ref));
end

