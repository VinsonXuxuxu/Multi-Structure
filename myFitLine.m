function xn = myFitLine(X)

Y = X(end, :)'; 
X = X(1:end-1, :)';
xn = X\Y; 