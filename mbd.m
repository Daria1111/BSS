%load test_data
%[u hr] = MCrestoration(G,[21 21]);

clc;

%signal = reshape(pim(1:900),[30,30]);
kernel1 = randn(10);
kernel2 = randn(10);

%blurred1 = reshape(csvread('y1.csv'),[30,30]);
%blurred2 = reshape(csvread('y2.csv'),[30,30]);
%I = checkerboard(100);
I = pim;
%I = reshape(I,[100,100]);
test = conv(I,kernel1(:),'same');
test0 = conv(I,kernel2(:),'same');
blurredcell = {reshape(test,[500,500]),reshape(test0,[500,500])};
[u hr] = MCrestoration(blurredcell, [10 10]);
