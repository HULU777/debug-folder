x=0:18;
y1 = 1-1./(x.^2+1);
y2 = 1-qfunc(x/6.6);
plot(x,y1);hold on;
plot(x,y2)
y11 = y1/sum(y1);y22 = y2/sum(y2);
plot(x,y11);hold on;
plot(x,y22);hold on;