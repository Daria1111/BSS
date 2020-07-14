function [X,E,errs,iter,times] = lrr2(y,A,lambda,K,N,opts)

% Solve the Low-Rank Representation minimization problem by M-ADMM
%
% min_{X,E} ||X||_*+lambda*loss(E), s.t. A=BX+E
% loss(E) = ||E||_1 or 0.5*||E||_F^2 or ||E||_{2,1}
%
% ---------------------------------------------
% Input:
%       A       -    d*na matrix
%       B       -    d*nb matrix
%       lambda  -    >0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.loss       -   'l1': loss(E) = ||E||_1 
%                               'l2': loss(E) = 0.5*||E||_F^2
%                               'l21' (default): loss(E) = ||E||_{2,1}
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.DEBUG      -   0 or 1
%
% Output:
%       X       -    nb*na matrix
%       E       -    d*na matrix
%       obj     -    objective function value
%       err     -    residual
%       iter    -    number of iterations
%
% version 1.0 - 18/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 

tol = 1e-4; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
DEBUG = 0;
loss = 'l21';

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'loss');        loss = opts.loss;            end
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end


[d,na] = size(y);
[a,nb] = size(A);

X = zeros(2*K,N);
E = zeros(d,na);
J = X;

Y1 = E;
Y2 = X;

Ablock = A(1:(a/2),1:(nb/2));
AtA = Ablock'*Ablock;
Aty = A'*y;
I = eye(nb/2);
invAtAIblock = (AtA+I)\I;
invAtAI = blkdiag(invAtAIblock,invAtAIblock);

iter = 0;
errs = [];
times = [];
tic;
start = tic;
for iter = 1 : max_iter
    Xk = X;
    Ek = E;
    Jk = J;
    % first super block {J,E}
    
    [J,nuclearnormJ] = prox_nuclear(X+Y2/mu,1/mu);
    
    
    E = prox_l21(y-A*X(:)+Y1/mu,lambda/mu);
    
        % second  super block {X}
    
    X0 = reshape(A'*(Y1/mu-E)+Aty,[2*K,N])-Y2/mu+J;
    X = reshape(invAtAI * X0(:),[2*K,N]);
    dY1 = y-A*X(:)-E;
    dY2 = X-J;
    
    chgX = max(max(abs(Xk-X)));
    chgE = max(max(abs(Ek-E)));
    chgJ = max(max(abs(Jk-J)));
    chg = max([chgX chgE chgJ max(abs(dY1(:))) max(abs(dY2(:)))]);
  
    
    if chg < tol
        break;
    end 
    Y1 = Y1 + mu*dY1;
    Y2 = Y2 + mu*dY2;
    mu = min(rho*mu,max_mu);    
    errs(iter) = norm(y-A*X(:),'fro');
    times(iter) = toc(start);
end

