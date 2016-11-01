% Modified "ut_sigmas" in ekfukf toolbox for the square-root UKF method
% The square-root UKF can be found in 
% Junjian Qi, Kai Sun, Jianhui Wang, and Hui Liu,
% "Dynamic State Estimation for Multi-Machine Power System by 
% Unscented Kalman Filter with Enhanced Numerical Stability,"
% IEEE Trans. Smart Grid, in press. DOI: 10.1109/TSG.2016.2580584 
% and the references therein
% Copyright (C) 2016 Junjian Qi
% 2015

%UT_SIGMAS - Generate Sigma Points for Unscented Transformation
%
% Syntax:
%   X = ut_sigmas(M,P,c);
%
% In:
%   M - Initial state mean (Nx1 column vector)
%   P - Initial state covariance
%   c - Parameter returned by UT_WEIGHTS
%
% Out:
%   X - Matrix where 2N+1 sigma points are as columns
%
% Description:
%   Generates sigma points and associated weights for Gaussian
%   initial distribution N(M,P). For default values of parameters
%   alpha, beta and kappa see UT_WEIGHTS.
%
% See also UT_WEIGHTS UT_TRANSFORM UT_SIGMAS
%

% Copyright (C) 2006 Simo Särkkä
%
% $Id: ut_sigmas.m 109 2007-09-04 08:32:58Z jmjharti $
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function X = ut_sigmas_sr(M,P,c);
  A = P';
  X = [zeros(size(M)) A -A];
  X = sqrt(c)*X + repmat(M,1,size(X,2));
