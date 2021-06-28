function [x_m_w_sweep] = add_sweep(sweep,x_m, y_m)
if sweep == 0
    x_m_w_sweep = x_m;

elseif y_m == 0
    x_m_w_sweep = 0;
        
else
    
    theta = atan(x_m/abs(y_m));
    x_m_w_sweep = abs(y_m)* tan(theta + sweep);
end
end