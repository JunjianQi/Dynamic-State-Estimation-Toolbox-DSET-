% Dynamical model function of power systems used for CKF
% Copyright (C) 2016 Junjian Qi
% Apr. 15, 2015
function x_n_all = power_system_f_tra1(x,para)
mac_con = para{1};
n_mac = size(mac_con,1);
xd = mac_con(:,6);
xdp = mac_con(:,7);
Td0p = mac_con(:,9);
xq = mac_con(:,11);
xqp = mac_con(:,12);
Tq0p = mac_con(:,14);
H2 = 2*mac_con(:,16);
wR = para{2};
Pm = para{3};
Efd = para{4};
dt = para{5};
Y_gprf = para{6};
mac_tra_idx = para{7};
mac_em_idx = para{8};
s_pos = para{9};
xcont = para{11};
mac_pot = 100*ones(n_mac,1)./mac_con(:,3);

n_col = size(x,2);
xaug = zeros(4*n_mac,n_col);
xaug(s_pos,:) = x;
xaug(setdiff(1:4*n_mac,s_pos),:) = repmat(xcont,1,n_col);

x_n_all = zeros(max(size(s_pos)),n_col);
for col=1:n_col
    dlt = xaug(1:n_mac,col);
    omg = xaug(n_mac+1:2*n_mac,col);
    eqp = xaug(2*n_mac+1:3*n_mac,col);
    edp = xaug(3*n_mac+1:end,col);

    [iqg,idg,eq,ed,Te] = calTemp(dlt,eqp,edp,mac_con,Y_gprf,mac_pot);
    ddlt = omg - wR;
    domg = wR*(Pm - Te - mac_con(:,17).*(omg-wR)/wR)./H2;
    deqp = zeros(n_mac,1);
    dedp = zeros(n_mac,1);
    deqp(mac_em_idx) = 0;
    dedp(mac_em_idx) = 0;
    deqp(mac_tra_idx) = (Efd(mac_tra_idx) - eqp(mac_tra_idx) - ...
        (xd(mac_tra_idx)-xdp(mac_tra_idx)).*idg(mac_tra_idx))./Td0p(mac_tra_idx);
    dedp(mac_tra_idx) = (-edp(mac_tra_idx) + (xq(mac_tra_idx) - ...
        xqp(mac_tra_idx)).*iqg(mac_tra_idx))./Tq0p(mac_tra_idx);
    dstate = [ddlt; domg; deqp; dedp];
    x_n1 = xaug(1:4*n_mac,col) + dt*dstate;

    dlt1 = x_n1(1:n_mac);
    omg1 = x_n1(n_mac+1:2*n_mac);
    eqp1 = x_n1(2*n_mac+1:3*n_mac);
    edp1 = x_n1(3*n_mac+1:4*n_mac);
    [iqg1,idg1,eq1,ed1,Te1] = calTemp(dlt1,eqp1,edp1,mac_con,Y_gprf,mac_pot);

    ddlt1 = omg1 - wR;
    domg1 = wR*(Pm - Te1 - mac_con(:,17).*(omg1-wR)/wR)./H2;
    deqp1 = zeros(n_mac,1);
    dedp1 = zeros(n_mac,1);
    deqp1(mac_em_idx) = 0;
    dedp1(mac_em_idx) = 0;
    deqp1(mac_tra_idx) = (Efd(mac_tra_idx) - eqp1(mac_tra_idx) - ...
        (xd(mac_tra_idx)-xdp(mac_tra_idx)).*idg1(mac_tra_idx))./Td0p(mac_tra_idx);
    dedp1(mac_tra_idx) = (-edp1(mac_tra_idx) + (xq(mac_tra_idx) - ...
        xqp(mac_tra_idx)).*iqg1(mac_tra_idx))./Tq0p(mac_tra_idx);

    dstate1 = [ddlt1; domg1; deqp1; dedp1];
    x_n = xaug(1:4*n_mac,col) + dt/2*(dstate + dstate1);
    x_n = x_n(s_pos);
    x_n_all(:,col) = x_n;
end



function [iqg,idg,eq,ed,Te] = calTemp(dlt,eqp,edp,mac_con,Y_gprf,mac_pot)
psi_re = sin(dlt).*edp + cos(dlt).*eqp;
psi_im = -cos(dlt).*edp + sin(dlt).*eqp;
psi = psi_re + 1i*psi_im;
iR = real(Y_gprf*psi);
iI = imag(Y_gprf*psi);
iq = iI.*sin(dlt) + iR.*cos(dlt);
id = iR.*sin(dlt) - iI.*cos(dlt);
idg = id.*mac_pot;
iqg = iq.*mac_pot;
eq = eqp - mac_con(:,7).*idg;
ed = edp + mac_con(:,7).*iqg;
Te = (ed.*id + eq.*iq).*mac_pot;

