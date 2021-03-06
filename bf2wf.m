function [bf_wf_transform_mat] = bf2wf(alfa,beta)
bf_wf_transform_mat = ...
    [cos(beta)*cos(alfa)  sin(beta)*cos(alfa) -sin(alfa); ...
    -sin(beta)            cos(beta)            0;          ...
     cos(beta)*sin(alfa)  sin(beta)*sin(alfa)  cos(alfa)];
    
end

