function [y] = Interpolar(x1,x2,x,y1,y2)

m=(y2-y1)/(x2-x1);
n = y2-m*x2;

y = m*x+n;

end