% Measurement model function for transient generator model
% Copyright (C) 2016 Junjian Qi
% Apr. 29, 2015
function Y = power_system_h_tra(x,para)
%% para
mac_con = para{1};
n_mac = size(mac_con,1);
mac_pot = 100*ones(n_mac,1)./mac_con(:,3);
Y_gprf = para{6};
Coutput = para{10};
%% states
dlt = x(1:n_mac);
eqp = x(2*n_mac+1:3*n_mac);
edp = x(3*n_mac+1:4*n_mac);
%% measurement
psi_re = sin(dlt).*edp + cos(dlt).*eqp;
psi_im = -cos(dlt).*edp + sin(dlt).*eqp;
psi = psi_re + 1i*psi_im;
iR = real(Y_gprf*psi);
iI = imag(Y_gprf*psi);
iq = iI.*sin(dlt) + iR.*cos(dlt);
id = iR.*sin(dlt) - iI.*cos(dlt);
% added for dealing with base of generators
idg = id.*mac_pot;
iqg = iq.*mac_pot;
eq = eqp - mac_con(:,7).*idg;
ed = edp + mac_con(:,7).*iqg;
eR = ed.*sin(dlt)+eq.*cos(dlt);
eI = eq.*sin(dlt)-ed.*cos(dlt);
Y = [eR(Coutput); eI(Coutput); iR(Coutput); iI(Coutput)]; % meas: eR eI iR iI

