clc;

K = 2;
L = 2;
N = 20;
m1 = pim(1:N);
m2 = u1(1:N);
C = randn(L,N);
B = randn(L,K);
h1 = randn(K,1);
h2 = randn(K,1);
h3 = randn(K,1);
h4 = randn(K,1);
A = circulant(C(:,1)) * B;
for iter = 2 : N
    elem = circulant(C(:,iter)) * B;
    A = cat(2,A,elem);
end
X11 = h1 * reshape(m1,[1,N]);
X12 = h2 * reshape(m1,[1,N]);
X1 = cat(1,X11,X12);
X21 = h3 * reshape(m2,[1,N]);
X22 = h4 * reshape(m2,[1,N]);
X2 = cat(1,X21,X22);
A = blkdiag(A,A);
y1 = A * X1(:);
y2 = A * X2(:);

y = y1 + y2 + 1000 * reshape(nf(1:2*L),[2*L,1]);



[Xk,E, errs, iter, times] = lrr2(y,A,lambda,K,N);
[u,s,v] = svds(Xk);


min(acos((v(:,1)' * m1')/(norm(v(:,1)) * norm(m1)))^2,acos((v(:,1)' * -m1')/(norm(v(:,1)) * norm(m1)))^2)

min(acos((v(:,2)' * m2')/(norm(v(:,2)) * norm(m2)))^2,acos((v(:,2)' * -m2')/(norm(v(:,2)) * norm(m2)))^2)

