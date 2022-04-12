clear all; close all; clc;
maxdistance = 10;
a = 0:maxdistance;
% a=[4 4 5 5 5 5 ];
square = 1./(a.^2+1);
Q = zeros(1,length(a));
for i = 1: length(a)
    Q(i) = qfunc(a(i)/3.796);
end
% square./Q
scaledsquare = scalingfitness(square);
scaledQ = scalingfitness(Q);

% figure;
subplot(2,1,1);
plot(square,'r*');hold on;
plot(Q,'g*');
% plot(square./Q,'B^')

% figure;
subplot(2,1,2);
plot(scaledsquare,'ro');hold on;
plot(scaledQ,'go');

