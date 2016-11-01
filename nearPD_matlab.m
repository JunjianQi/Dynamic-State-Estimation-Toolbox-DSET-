% nearPD_matlab
% Author: Junjian Qi
function [X,normF,iter,converged] = nearPD_matlab(x,corr,keepDiag,do2eigen,doSym,doDykstra,onlyvalues,eigtol,convtol,posdtol,maxit)
if nargin < 11
    [corr,keepDiag,do2eigen,doSym,doDykstra,onlyvalues,eigtol,convtol,posdtol,maxit] = defaultPara(nargin);
end
if ~issymmetric(x)
    x = (x+x')/2;
end
n = size(x,2);
if keepDiag
    diagX0 = diag(x);
end
if doDykstra
    D_S = x;
    D_S(:,:) = 0;
end
X = x;
iter = 0;
converged = 0;
while iter < maxit && ~converged
    Y = X;
    if doDykstra
        R = Y - D_S;
        [Q,d] = mexeig(R);
    else
        [Q,d] = mexeig(Y);
    end
    p = d > eigtol * max(d);
    Q = Q(:,p);
    X = Q.*repmat(d(p)',n,1)*Q';
    if doDykstra
        D_S = X - R;
    end
    if doSym
        X = (X + X')/2;
    end
    if corr
        X(1:n+1:n*n) = 1;
    elseif keepDiag
        X(1:n+1:n*n) = diagX0;
    end
    conv = norm(Y-X,'fro')/norm(Y,'fro');
    iter = iter + 1;
    converged = conv <= convtol;
end

if (~converged)
    disp('not converged!');
end
if do2eigen || onlyvalues
    [Q,d] = mexeig(X);
    Eps = posdtol*abs(max(d));
    if min(d) < Eps
        d(d<Eps) = Eps;
        if ~onlyvalues
            diagX = diag(X);
            X = Q*(diag(d)*Q');
            D = sqrt(max(Eps,diagX)./diag(X));
            X = diag(D)*X.*repmat(D,1,n);
        end
        if onlyvalues
            X = diag(d);
        end
        if corr
            X(1:n+1:n*n) = 1;
        elseif keepDiag
            X(1:n+1:n*n) = diagX0;
        end
    end
    normF = norm(x-X,'fro');
else
    normF = norm(x-X,'fro');
end
X = (X + X')/2;

%%% function 
function [corr,keepDiag,do2eigen,doSym,doDykstra,onlyvalues,eigtol,convtol,posdtol,maxit] = defaultPara(nargin)
if (nargin<2)
    corr = 0; keepDiag = 0; do2eigen = 1; doSym = 0; doDykstra = 1; onlyvalues = 0;
    eigtol = 1e-6; convtol = 1e-7; posdtol = 1e-7; maxit = 100;
end
% if (nargin<3)
%     keepDiag = 0; do2eigen = 1; doSym = 0; doDykstra = 1; onlyvalues = 0;
%     eigtol = 1e-6; convtol = 1e-7; posdtol = 1e-8; maxit = 100;
% end
% if (nargin<4)
%     do2eigen = 1; doSym = 0; doDykstra = 1; onlyvalues = 0;
%     eigtol = 1e-6; convtol = 1e-7; posdtol = 1e-8; maxit = 100;
% end
% if (nargin<5)
%     doSym = 0; doDykstra = 1; onlyvalues = 0; eigtol = 1e-6; 
%     convtol = 1e-7; posdtol = 1e-8; maxit = 100;
% end
% if (nargin<6)
%     doDykstra = 1; onlyvalues = 0; eigtol = 1e-6; convtol = 1e-7; 
%     posdtol = 1e-8; maxit = 100;
% end
% if (nargin<7)
%     onlyvalues = 0; eigtol = 1e-6; convtol = 1e-7; posdtol = 1e-8; 
%     maxit = 100;
% end
% if (nargin<8)
%     eigtol = 1e-6; convtol = 1e-7; posdtol = 1e-8; maxit = 100;
% end
% if (nargin<9)
%     convtol = 1e-7; posdtol = 1e-8; maxit = 100;
% end
% if (nargin<10)
%     posdtol = 1e-8; maxit = 100;
% end
% if (nargin<11)
%     maxit = 100;
% end

