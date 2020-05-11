
figure
plot(Source(1,:))
figure
plot(Source(2,:))
figure
plot(signal1)
figure
plot(signal2)

acos((Source(2,:) * -signal2)/(norm(Source(2,:)) * norm(signal2)))^2