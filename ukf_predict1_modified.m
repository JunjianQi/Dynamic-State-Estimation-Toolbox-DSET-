% Modified "ukf_predict1" in ekfukf toolbox to get the modified UKF method
% The modified UKF can be found in 
% Junjian Qi, Kai Sun, Jianhui Wang, and Hui Liu,
% "Dynamic State Estimation for Multi-Machine Power System by 
% Unscented Kalman Filter with Enhanced Numerical Stability,"
% IEEE Trans. Smart Grid, in press. DOI: 10.1109/TSG.2016.2580584 
% and the references therein
% Copyright (C) 2016 Junjian Qi
% 2015

%UKF_PREDICT1  Nonaugmented (Additive) UKF prediction step
%
% Syntax:
%   [M,P] = UKF_PREDICT1(M,P,f,Q,f_param,alpha,beta,kappa,mat)
%
% In:
%   M - Nx1 mean state estimate of previous step
%   P - NxN state covariance of previous step
%   f - Dynamic model function as a matrix A defining
%       linear function a(x) = A*x, inline function,
%       function handle or name of function in
%       form a(x,param)                   (optional, default eye())
%   Q - Process noise of discrete model   (optional, default zero)
%   f_param - Parameters of f               (optional, default empty)
%   alpha - Transformation parameter      (optional)
%   beta  - Transformation parameter      (optional)
%   kappa - Transformation parameter      (optional)
%   mat   - If 1 uses matrix form         (optional, default 0)
%
% Out:
%   M - Updated state mean
%   P - Updated state covariance
%
% Description:
%   Perform additive form Unscented Kalman Filter prediction step.
%
%   Function a should be such that it can be given
%   DxN matrix of N sigma Dx1 points and it returns 
%   the corresponding predictions for each sigma
%   point. 
%
% See also:
%   UKF_UPDATE1, UKF_PREDICT2, UKF_UPDATE2, UKF_PREDICT3, UKF_UPDATE3,
%   UT_TRANSFORM, UT_WEIGHTS, UT_MWEIGHTS, UT_SIGMAS

% Copyright (C) 2003-2006 Simo S�rkk�
%
% $Id: ukf_predict1.m 483 2010-10-18 08:54:19Z jmjharti $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function [M,P,D] = ukf_predict1_modified(M,P,f,Q,f_param,alpha,beta,kappa,mat)

  %
  % Check which arguments are there
  %
  if nargin < 2
    error('Too few arguments');
  end
  if nargin < 3
    f = [];
  end
  if nargin < 4
    Q = [];
  end
  if nargin < 5
    f_param = [];
  end
  if nargin < 6
    alpha = [];
  end
  if nargin < 7
    beta = [];
  end
  if nargin < 8
    kappa = [];
  end
  if nargin < 9
    mat = [];
  end

  %
  % Apply defaults
  %
  if isempty(f)
    f = eye(size(M,1));
  end
  if isempty(Q)
    Q = zeros(size(M,1));
  end
  if isempty(mat)
    mat = 0;
  end

  %
  % Do transform
  % and add process noise
  %
  
  tr_param = {alpha beta kappa mat};
  [M,P,D] = ut_transform_modified(M,P,f,f_param,tr_param);
  P = P + Q;

