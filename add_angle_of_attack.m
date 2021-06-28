function [x_m_aoa, z_m_aoa] = add_angle_of_attack(aoa,x_m)

% add the angle of attack, note the z_component is at this point still zero
aoa = -aoa;  % since we are pointing the trailing edge down and leaving
% as the twisting point
x_m_aoa = x_m*cos(aoa);
z_m_aoa = x_m*sin(aoa);

end

