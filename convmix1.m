clc;
Q = 2;              % nb outputs
R = 2;              % nb sources
%N = 1e3;            % nb samples
filterlength = 10;  % Filter length for sources

signal1 = pim(1:900)';
signal2 = u1(1:900)';
noise1 = nf(1:900)';
noise2 = nf(901:1800)';
s = cat(2,signal1,signal2);
filt1 = randn(filterlength,Q);
x1 = conv(s(:,1),filt1(:,1),'same');
for q = 2:Q
    x1 = x1 + conv(s(:,1),filt1(:,q),'same');
end

filt2 = randn(filterlength,Q);
x2 = conv(s(:,2),filt2(:,1),'same');
for q = 2:Q
    x2 = x2 + conv(s(:,2),filt2(:,q),'same');
end
x1 = x1+noise1;
x2 = x2 +noise2;
x = cat(2,x1,x2);
tic
[Source,Contribution] = Deflation(x,ConvolutiveMixtureParameters);
toc