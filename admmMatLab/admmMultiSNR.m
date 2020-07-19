clc;

K = 10;
L = 1000;
N = 500;
m = pim(1:N);
m1 = u1(1:N);
m2 = u2(1:N);
C = randn(L,N);
B = randn(L,K);
%s = pim(1:100);
%n = nf(1:100);
%N = size(s,2);
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
A = circulant(C(:,1)) * B;
for iter = 2 : N
    elem = circulant(C(:,iter)) * B;
    A = cat(2,A,elem);
end
A = blkdiag(A,A);
y1 = A * X1(:);
y2 = A * X2(:);

y = y1 + y2;
lambda = 0.001;
[Xk,E, errs, iter, times] = lrr2(y,A,lambda,K,N);
[u,s,v] = svds(Xk);


acos((v(:,1)' * -m1')/(norm(v(:,1)) * norm(m1)))^2
figure
plot(u1(1:N))
figure
plot(v(:,1))
figure
snr(v(:,1))
acos((v(:,2)' * m')/(norm(v(:,2)) * norm(m)))^2
figure
plot(pim(1:N))
figure
plot(v(:,2))
figure
snr(v(:,2))