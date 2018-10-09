% Dynamic State Estimation (DSE) by using 
% EKF, UKF-GPS, SR-UKF, UKF-schol, UKF-kappa, UKF-modified, and UKF-DeltaQ
% Copyright (C) 2016 Junjian Qi
% May 21, 2015
% Oct. 28, 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   References:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [1] Junjian Qi, Kai Sun, Jianhui Wang, and Hui Liu,
% "Dynamic State Estimation for Multi-Machine Power System by 
% Unscented Kalman Filter with Enhanced Numerical Stability,"
% IEEE Trans. Smart Grid, in press. DOI: 10.1109/TSG.2016.2580584 
% [2] Junjian Qi, Kai Sun, and Wei Kang, 
% "Optimal PMU Placement for Power System Dynamic State Estimation 
% by Using Empirical Observability Gramian,"
% IEEE Trans. Power Systems, vol. 30, no. 4, pp. 2041-2054, Jul. 2015.
% [3] Junjian Qi, Kai Sun, and Wei Kang, 
% "Adaptive optimal PMU placement based on empirical observability gramian," 
% 10th IFAC Symposium on Nonlinear Control Systems (NOLCOS), Monterey, CA USA, Aug. 2016."
% [4] Kai Sun, Junjian Qi, and Wei Kang, 
% "Power system observability and dynamic state estimation for stability 
% monitoring using synchrophasor measurements," Control Engineering Practice, 
% vol. 53, pp. 160-172, Aug. 2016. DOI: 10.1016/j.conengprac.2016.01.013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test_UKFGPS_SRUKF
clear global
clear
close all
clc
%% import data files
system = 3; % 3 for 3-machine system; 48 for 48-machine system
if system == 3
    load('x0_3m.mat');  % post-fault states (the fault is a three-phase fault at one bus of the line with highest power flow)
    load('x00_3m.mat'); % pre-fault states
    load('para_pf2_3m.mat');
    mac_tra_idx = para_pf2{7};
    mac_em_idx = para_pf2{8};
    n_mac = 3;
    n_pmu = 1;
    sensorpos_opt_det = 3;
elseif system == 48
    load('x0_48m.mat');  % post-fault states (the fault is a three-phase fault at one bus of the line with highest power flow)
    load('x00_48m.mat'); % pre-fault states
    load('para_pf2_48m.mat');
    mac_tra_idx = para_pf2{7};
    mac_em_idx = para_pf2{8};
    n_mac = 48;
    load('Placement_48m.mat');
    opt_det = x1_basic;
    n_pmu = 24;
    sensorpos_opt_det = find(opt_det(n_pmu,:)~=0);    
end
%% para
%% main part
index_record = zeros(7,8,1);
time_record_ukf_gps_record = [];
time_r = [];
idx = 0;
for i=n_pmu
    idx = idx + 1;
    n_tra = size(mac_tra_idx,1);
    n_s = 4*n_tra + 2*size(mac_em_idx,1);
    noisedltbd = 0.5*pi/180;
    noisewbd = 1e-3*2*pi*60;
    noiseeqpbd = 1e-3;
    noiseedpbd = 1e-3;
    P = [noisewbd^2*eye(n_mac) zeros(n_mac) zeros(n_mac,n_tra) zeros(n_mac,n_tra); ...
          zeros(n_mac) noisedltbd^2*eye(n_mac) zeros(n_mac,n_tra) zeros(n_mac,n_tra);...
          zeros(n_tra,n_mac) zeros(n_tra,n_mac) noiseeqpbd^2*eye(n_tra) zeros(n_tra);...
          zeros(n_tra,n_mac) zeros(n_tra,n_mac) zeros(n_tra) noiseedpbd^2*eye(n_tra)];    
    [index,time_record_ukf_gps,needTime,time] = testcase(n_mac,n_s,x00,x0,sensorpos_opt_det,para_pf2,P);
    index_record(:,:,idx) = index;
    time_r = [time_r; time];
    time_record_ukf_gps_record = [time_record_ukf_gps_record; time_record_ukf_gps];
    save(strcat('.\index_KF_',num2str(n_mac),'m.mat'), ...
        'index_record','index','time_r','time_record_ukf_gps','needTime','time');
    % In 'index_record.mat', each row corresponds to EKF, UKF-schol, UKF-GPS,
    % SR-UKF, UKF-kappa, UKF-modified, UKF-DeltaQ
end

%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%
function [index,time_record_ukf_gps,needTime,time] = testcase(n_mac,n_s,x00,x0,sensorpos_opt_det,para,P)
%% simulate
t0 = 0;
Tfinal = 10;
sensorfreq = 120;
tsequence = t0:1/sensorfreq:Tfinal;
%% slove ODE
para{5} = 1/sensorfreq;
s_pos = para{9};
states_nonoise = zeros(size(tsequence,2),4*n_mac);
states_nonoise(1,:) = x0';
for i=2:size(tsequence,2)
    states_nonoise(i,:) = power_system_f_tra(states_nonoise(i-1,:)',para)';
end
states_nonoise = states_nonoise(:,s_pos);
Q = diag((0.1*max(abs(diff(states_nonoise)))).^2);
states = zeros(size(tsequence,2),4*n_mac);
states(1,:) = x0';
for i=2:size(tsequence,2)
    states(i,:) = power_system_f_tra(states(i-1,:)',para)';
    processNoise = gauss_rnd(zeros(1,n_s)', Q)';
    states(i,s_pos) = states(i,s_pos) + processNoise;
end
%% set initial value
M = x00;
%% dynamic state estimation
para{5} = 1/60;
[index,time_record_ukf_gps,needTime,time] = ekf_ukf(n_mac,sensorpos_opt_det,states,M,P,Q,para,tsequence,Tfinal);


%%%%%%%%%%%%%%%% ekf_ukf %%%%%%%%%%%%%%%%%%
function [index,time_record_ukf_gps,needTime,time] = ekf_ukf(n_mac,sensorpos,states,M0,P0,Q,para,tsequence,Tfinal)
%% measurements
sensorfreq = 1/para{5};
ratio = (size(tsequence,2)-1)/sensorfreq/Tfinal;
tsequence = tsequence(1:ratio:end);
n_sensor = size(sensorpos,2);
%% create measurements with noise
f_func = @power_system_f_tra1;
h_func = @power_system_h_tra1;
para{10} = sensorpos';
states = states(1:ratio:end,:);
Y_real = zeros(4*n_sensor,size(states,1)); % meas: eR eI iR iI
for i=1:size(states,1)
    Y_real(:,i) = power_system_h_tra(states(i,:)',para);
end
R = diag(0.01^2*ones(1,4*n_sensor));
Y = Y_real;
for i=1:4*n_sensor % meas: eR eI iR iI
    Y(i,:) = Y(i,:) + sqrt(R(i,i)).*randn(1,size(Y,2));
end
%% estimate with UKF-GPS
s_pos = para{9};
n_s = size(s_pos,1);
U_MM_gps = zeros(n_s,size(Y,2));
U_MM_gps(:,1) = M0(s_pos);
U_PP_gps = zeros(n_s,n_s,size(Y,2));
M = M0(s_pos);
para{11} = M0(setdiff(1:4*n_mac,s_pos));
P = P0;
tPD = 0;
flag_p = 0;
num_solve = 0;
iteration_r = 0;
norm_r = 0;
allconverged = 1;
needTime = [];
tstart1 = tic;
for k=2:size(Y,2)
   [M,P] = ukf_predict1(M,P,f_func,Q,para);
   [~,p] = chol(P);
   if p
       num_solve = num_solve + 1;
       needTime = [needTime; tsequence(k) 1];
       flag_p = 1;
       tstart2 = tic;
       P_tmp = P;
       [P,normF,iterations,converged] = nearPD_matlab(P_tmp);
       if converged == 0
           allconverged = 0;
       end
       iteration_r = iteration_r + iterations;
       norm_r = norm_r + normF;
       tPD = tPD + toc(tstart2);
   end
   [M,P] = ukf_update1(M,P,Y(:,k),h_func,R,para);
   [~,p] = chol(P);
   if p
       num_solve = num_solve + 1;
       needTime = [needTime; tsequence(k) 2];
       flag_p = 1;
       P_tmp = P;
       tstart2 = tic;
       [P,normF,iterations,converged] = nearPD_matlab(P_tmp);
       if converged == 0
           allconverged = 0;
       end
       iteration_r = iteration_r + iterations;
       norm_r = norm_r + normF;
       tPD = tPD + toc(tstart2);
   end
   U_MM_gps(:,k) = M;
   U_PP_gps(:,:,k) = P;
end
time_ukfgps = toc(tstart1);

time_record_ukf_gps = [num_solve tPD iteration_r/num_solve norm_r/num_solve allconverged flag_p];

%% estimate with SR-UKF
tstart2 = tic;
s_pos = para{9};
n_s = size(s_pos,1);
U_MM_sq = zeros(n_s,size(Y,2));
U_MM_sq(:,1) = M0(s_pos);
U_PP_sq = zeros(n_s,n_s,size(Y,2));
M_sq = M0(s_pos);
para{11} = M0(setdiff(1:4*n_mac,s_pos));
P_sq = chol(P0);
for k=2:size(Y,2)
   [M_sq,P_sq] = ukf_predict1_sr(M_sq,P_sq,f_func,Q,para);
   [M_sq,P_sq] = ukf_update1_sr(M_sq,P_sq,Y(:,k),h_func,R,para);
   U_MM_sq(:,k) = M_sq;
   U_PP_sq(:,:,k) = P_sq;
end
time_srukf = toc(tstart2);
%% estimtate with UKF - (UKF-schol)
tstart3 = tic;
s_pos = para{9};
n_s = size(s_pos,1);
U_MM = zeros(n_s,size(Y,2));
U_MM(:,1) = M0(s_pos);
U_PP = zeros(n_s,n_s,size(Y,2));
M = M0(s_pos);
para{11} = M0(setdiff(1:4*n_mac,s_pos));
P = P0;
for k=2:size(Y,2)
   [M,P] = ukf_predict1(M,P,f_func,Q,para);
   [M,P] = ukf_update1(M,P,Y(:,k),h_func,R,para);
   U_MM(:,k) = M;
   U_PP(:,:,k) = P;
end
time_ukf = toc(tstart3);
%% estimtate with UKF1 - (UKF-kappa)
tstart4 = tic;
s_pos = para{9};
n_s = size(s_pos,1);
U_MM1 = zeros(n_s,size(Y,2));
U_MM1(:,1) = M0(s_pos);
U_PP1 = zeros(n_s,n_s,size(Y,2));
M = M0(s_pos);
para{11} = M0(setdiff(1:4*n_mac,s_pos));
P = P0;
alpha = 1;
beta = 0;
kappa = 0;
stable_ukf1 = 1;
for k=2:size(Y,2)
   [M,P] = ukf_predict1(M,P,f_func,Q,para,alpha,beta,kappa);
   [~,p] = chol(P);
   if p
       stable_ukf1 = 0;
   end   
   [M,P] = ukf_update1(M,P,Y(:,k),h_func,R,para,alpha,beta,kappa);
   [~,p] = chol(P);
   if p
       stable_ukf1 = 0;
   end   
   U_MM1(:,k) = M;
   U_PP1(:,:,k) = P;
end
time_ukf1 = toc(tstart4);
%% estimtate with UKF2 - (UKF-modified)
tstart5 = tic;
s_pos = para{9};
n_s = size(s_pos,1);
U_MM2 = zeros(n_s,size(Y,2));
U_MM2(:,1) = M0(s_pos);
U_PP2 = zeros(n_s,n_s,size(Y,2));
M = M0(s_pos);
para{11} = M0(setdiff(1:4*n_mac,s_pos));
P = P0;
stable_ukf2 = 1;
for k=2:size(Y,2)
   [M,P] = ukf_predict1_modified(M,P,f_func,Q,para);
   [~,p] = chol(P);
   if p
       stable_ukf2 = 0;
   end
   [M,P] = ukf_update1_modified(M,P,Y(:,k),h_func,R,para);
   [~,p] = chol(P);
   if p
       stable_ukf2 = 0;
   end  
   U_MM2(:,k) = M;
   U_PP2(:,:,k) = P;
end
time_ukf2 = toc(tstart5);
%% estimtate with UKF3 - (UKF-DeltaQ)
tstart6 = tic;
s_pos = para{9};
n_s = size(s_pos,1);
U_MM3 = zeros(n_s,size(Y,2));
U_MM3(:,1) = M0(s_pos);
U_PP3 = zeros(n_s,n_s,size(Y,2));
M = M0(s_pos);
para{11} = M0(setdiff(1:4*n_mac,s_pos));
P = P0;
stable_ukf3 = 1;
for k=2:size(Y,2)
   [M,P] = ukf_predict1(M,P,f_func,Q,para);
   P = P + 0.005^2*eye(n_s);
   [~,p] = chol(P);
   if p
       stable_ukf3 = 0;
   end
   [M,P] = ukf_update1(M,P,Y(:,k),h_func,R,para);
   [~,p] = chol(P);
   if p
       stable_ukf3 = 0;
   end
   U_MM3(:,k) = M;
   U_PP3(:,:,k) = P;
end
time_ukf3 = toc(tstart6);
%% estimate with EKF
tstart7 = tic;
df_dx_func = @power_system_df_dx;
dh_dx_func = @power_system_dh_dx;
M = M0(s_pos);
P = P0;
E_MM = zeros(n_s,size(Y,2));
E_MM(:,1) = M0(s_pos);
E_PP = zeros(n_s,n_s,size(Y,2));
for k=2:size(Y,2)
   [M,P] = ekf_predict1(M,P,df_dx_func,Q,f_func,[],para);
   [M,P] = ekf_update1(M,P,Y(:,k),dh_dx_func,R,h_func,[],para);
   E_MM(:,k)   = M;
   E_PP(:,:,k) = P;
end
time_ekf = toc(tstart7);

time = [time_ekf time_ukf time_ukfgps time_srukf time_ukf1 time_ukf2 time_ukf3];

index1 = error(states,E_MM,para,tsequence,s_pos,n_mac);
index2 = error(states,U_MM,para,tsequence,s_pos,n_mac);
index3 = error(states,U_MM_gps,para,tsequence,s_pos,n_mac);
index4 = error(states,U_MM_sq,para,tsequence,s_pos,n_mac);
index5 = error(states,U_MM1,para,tsequence,s_pos,n_mac);
index6 = error(states,U_MM2,para,tsequence,s_pos,n_mac);
index7 = error(states,U_MM3,para,tsequence,s_pos,n_mac);
index = [index1; index2; index3; index4; index5; index6; index7];

save(strcat('.\KF_',num2str(n_mac),'m.mat'),...
    'states','U_MM_gps','U_MM_sq','U_MM','U_MM1','U_MM2','U_MM3','E_MM','tsequence',...
    's_pos','n_mac','para','stable_ukf1','stable_ukf2');


%%%%%%%%%%%%%%%%%%%% error %%%%%%%%%%%%%%%%%%%%%%
function index = error(states,U_MM,para,tsequence,s_pos,n_mac)
U_MM1 = U_MM';
states = states(:,s_pos);
dltX = U_MM1 - states;
index_delta = 0;
index_omega = 0;
index_eqp = 0;
index_edp = 0;
mac_tra_idx = para{7};
for i=1:size(tsequence,2)
    index_delta = index_delta + 1/n_mac*sum((dltX(i,1:n_mac)).^2);
    index_omega = index_omega + 1/n_mac*sum((dltX(i,n_mac+1:2*n_mac)).^2);
    index_eqp = index_eqp + 1/size(mac_tra_idx,1)*sum((dltX(i,2*n_mac+1:2*n_mac+size(mac_tra_idx,1))).^2);
    index_edp = index_edp + 1/size(mac_tra_idx,1)*sum((dltX(i,2*n_mac+size(mac_tra_idx,1)+1:end)).^2);    
end
index_delta = sqrt(index_delta /(size(tsequence,2)));
index_omega = sqrt(index_omega /(size(tsequence,2)));
index_eqp = sqrt(index_eqp /(size(tsequence,2)));
index_edp = sqrt(index_edp /(size(tsequence,2)));
%% number of convergent states
index_num_dlt = zeros(1,n_mac);
index_num_omg = zeros(1,n_mac);
index_num_eqp = zeros(1,size(mac_tra_idx,1));
index_num_edp = zeros(1,size(mac_tra_idx,1));
t_4 = find(tsequence==9);
for i=1:n_mac
    if all(abs(dltX(t_4:end,i)) < abs(states(t_4:end,i))*0.01)
        index_num_dlt(1,i) = 1;
    end
    if all(abs(dltX(t_4:end,n_mac+i)) < abs(states(t_4:end,n_mac+i))*0.01)
        index_num_omg(1,i) = 1;
    end
end
for i=1:size(mac_tra_idx,1)
    if all(abs(dltX(t_4:end,2*n_mac+i)) < abs(states(t_4:end,2*n_mac+i))*0.01)
        index_num_eqp(1,i) = 1;
    end
    if all(abs(dltX(t_4:end,2*n_mac+size(mac_tra_idx,1)+i)) < abs(states(t_4:end,2*n_mac+size(mac_tra_idx,1)+i))*0.01)
        index_num_edp(1,i) = 1;
    end
end
index = [index_delta index_omega index_eqp index_edp sum(index_num_dlt) ...
    sum(index_num_omg) sum(index_num_eqp) sum(index_num_edp)];

