function [x] = chord_x(y, c_0, lambda,b)
%x = c_0 - (1-lambda)*abs(2*c_0*y/b);
x = c_0*(1 - 2*(1-lambda)*abs(y/b));

end

