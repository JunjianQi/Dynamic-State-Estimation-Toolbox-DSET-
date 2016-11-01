% Jacobian of the state transition function
% Copyright (C) 2016 Junjian Qi
% Aug. 8, 2014
function df = power_system_df_dx(x,para)
mac_con = para{1};
n_mac = size(mac_con,1);
n_s = size(x,1);
mac_tra_idx = para{7};
df = zeros(n_s,n_s);
n_tra = size(mac_tra_idx,1);
F = power_system_f_tra1(x,para);

for i=1:2*n_mac+2*n_tra
    x_tmp = x;
    Delta_x = abs(x_tmp(i))*1e-4;
    x_tmp(i) = x_tmp(i) + Delta_x;
    F1 = power_system_f_tra1(x_tmp,para);
    df(:,i) = (F1-F)/Delta_x;
end

