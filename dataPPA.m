clc;
%{
m = pim(1:500);
K = 10;
L = 1000;
N = 500;
C = randn(L,N);
h = randn(K,1);
B = randn(L,K);
X = h * reshape(m,[1,500]);

A0 = circulant(C(:,1)) * B;
for iter = 2 : N
    elem = circulant(C(:,iter)) * B;
    A0 = sparse(cat(2,A0,elem));
end
%}
y = A(X,A0);
b = y + 8;
delta = norm(y-b);
[nr,nc] = size(X);
m1 = 1000;
m2 = 0;
m3 = 0;
amap = @A;
atmap = @At;
[X1,iter,ttime,sd,runhist]=dualPPA(amap, atmap,A0,b,delta,nr,nc,m1,m2,m3)
%[X1,iter,ttime,sd,runhist]=primPPA(amap, atmap,A0,b,delta,nr,nc,m1,m2,m3)
