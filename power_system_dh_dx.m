% Jacobian of the measurement model function
% Copyright (C) 2016 Junjian Qi
% Aug. 7, 2014
function dY = power_system_dh_dx(x,para)
mac_con = para{1};
n_mac = size(mac_con,1);
mac_tra_idx = para{7};
mac_em_idx = para{8};
Coutput = para{10};
n_s = 4*size(mac_tra_idx,1) + 2*size(mac_em_idx,1);
dY = zeros(4*size(Coutput,1),n_s);
Y = power_system_h_tra1(x,para);
n_tra = size(mac_tra_idx,1);

for i=1:2*n_mac+2*n_tra
    x_tmp = x;
    Delta_x = abs(x_tmp(i))*1e-4;
    x_tmp(i) = x_tmp(i) + Delta_x;
    Y1 = power_system_h_tra1(x_tmp,para);
    dY(:,i) = (Y1-Y)/Delta_x;
end

