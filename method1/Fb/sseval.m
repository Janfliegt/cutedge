function sse = sseval(x, xdata, ydata)
% y = A * exp(-B*x + C)
% https://ch.mathworks.com/help/matlab/math/example-curve-fitting-via-optimization.html
A = x(1);
B = x(2);
C = x(3);
sse = sum((ydata - A * exp(-B*xdata + C)).^2);