function [V] = horse(Gamma, xc, xa, xb)
a = xc - xa;
b = xc - xb;
eps = 2.2204e-16;

if abs((norm(a) * norm(b) + a*b')) < eps
    Cb = 0;
else
    Cb = cross(a,b)/(norm(a)* norm(b) + a*b') * (1/norm(a) + 1/norm(b));
end % if else

% Left Trailing Vortex
Ctl = +cross(a, [1 0 0]) / (norm(a) - a*[1 0 0]') * 1/norm(a);

% Right Trailing Vortex
Ctr = -cross(b, [1 0 0]) / (norm(b) - b*[1 0 0]') * 1/norm(b);

V = Gamma/(4*pi) * (Cb + Ctl + Ctr);
end
% [V] = horse(1,1, xc, xa, xb)
