clc;

s1 = pim(1:500);
s2 = u1(1:500);
K = 10;
L = 2000;
N = 500;
C = randn(L,N);
h1 = randn(K,1);
h2 = randn(K,1);
h3 = randn(K,1);
h4 = randn(K,1);
B = randn(L,K);
A00 = circulant(C(:,1)) * B;
for iter = 2 : N
    elem = circulant(C(:,iter)) * B;
    A00 = cat(2,A00,elem);
end
X11 = h1 * reshape(s1,[1,N]);
X12 = h2 * reshape(s1, [1,N]);
X1 = cat(1,X11,X12);
X21 = h3 * reshape(s2, [1,N]);
X22 = h4 * reshape(s2, [1,N]);
X2 = cat(1,X21,X22);
X = X1 + X2;
A0 = blkdiag(A00,A00);

y = A(X,A0);
b = y + 0.01;
delta = norm(y-b);
[nr,nc] = size(X);
m1 = 4000;
m2 = 0;
m3 = 0;
amap = @A;
atmap = @At;
[X1,iter,ttime,sd,runhist]=dualPPA(amap, atmap,A0,b,delta,nr,nc,m1,m2,m3)