clc;

[u,s,v] = svd(X1);
acos((v(:,2)' * -s2')/(norm(v(:,2)) * norm(s2)))^2
figure
plot(u1(1:500))
figure
plot(v(:,2))
figure
snr(v(:,2))