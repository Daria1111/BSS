clc;

m = pim(1:500);
m1 = u1(1:500);
K = 10;
L = 2000;
N = 500;
C = randn(L,N);
B = randn(L,K);
A00 = circulant(C(:,1)) * B;
for iter = 2 : N
    elem = circulant(C(:,iter)) * B;
    A00 = cat(2,A00,elem);
end
nosensors = K;
%nosources = size(s,1);
snrdef1 = [15 20];   % SNR for source 1
snrdef2 = [20 28];  % SNR for source 2
[X11,h1,snrgen1] = gen_model2(m,nosensors,snrdef1,nf(1:N));
[X12,h2,snrgen2] = gen_model2(m1,nosensors,snrdef2,nf(N:2*N));
X1 = cat(1,X11,X12);
[X21,h3,snrgen3] = gen_model2(m,nosensors,snrdef1,nf(1:N));
[X22,h4,snrgen4] = gen_model2(m1,nosensors,snrdef2,nf(N:2*N));
X2 = cat(1,X21,X22);
X = X1 + X2;
A0 = blkdiag(A00,A00);

y = A(X,A0);
b = y + 10;
delta = norm(y-b);
[nr,nc] = size(X);
m1 = 4000;
m2 = 0;
m3 = 0;
amap = @A;
atmap = @At;
[X_opt,iter,ttime,sd,runhist]=dualPPA(amap, atmap,A0,b,delta,nr,nc,m1,m2,m3);

[u,s,v] = svd(X_opt);
m1 = u1(1:500);
acos((v(:,1)' * -m1')/(norm(v(:,1)) * norm(m1)))^2
figure
plot(u1(1:500))
figure
plot(v(:,1))
figure
snr(v(:,1))

acos((v(:,2)' * -m')/(norm(v(:,2)) * norm(m)))^2
figure
plot(pim(1:500))
figure
plot(v(:,2))
figure
snr(v(:,2))

