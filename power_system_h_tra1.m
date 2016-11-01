% Measurement model function for transient generator model used for CKF
% Copyright (C) 2016 Junjian Qi
% Apr. 15, 2015

function Y_all = power_system_h_tra1(x,para)
%% para
mac_con = para{1};
n_mac = size(mac_con,1);
Y_gprf = para{6};
mac_pot = 100*ones(n_mac,1)./mac_con(:,3);
Coutput = para{10};
s_pos = para{9};
xcont = para{11};

n_col = size(x,2);
xaug = zeros(4*n_mac,n_col);
xaug(s_pos,:) = x;
xaug(setdiff(1:4*n_mac,s_pos),:) = repmat(xcont,1,n_col);

Y_all = zeros(size(Coutput,1)*4,n_col);
%% states
for col=1:n_col
    dlt = xaug(1:n_mac,col);
    eqp = xaug(2*n_mac+1:3*n_mac,col);
    edp = xaug(3*n_mac+1:4*n_mac,col);
    %% measurement
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
    eR = ed.*sin(dlt)+eq.*cos(dlt);
    eI = eq.*sin(dlt)-ed.*cos(dlt);
    Y = [eR(Coutput); eI(Coutput); iR(Coutput); iI(Coutput)]; % meas: eR eI iR iI
    Y_all(:,col) = Y;
end

