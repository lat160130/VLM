function [y_m__wdihedral,z_m_w_dihedral] = add_dihedral(y_m,z_m, dihedral)

if rad2deg(dihedral) >= 90
    error("dihedral < 90");
elseif dihedral == 0
    y_m__wdihedral = y_m;
    z_m_w_dihedral = z_m;
    
elseif y_m == 0
    y_m__wdihedral = 0;
    z_m_w_dihedral = 0;
else
    y_m = abs(y_m);
    y_m__wdihedral = y_m * cos(dihedral);
    z_m_w_dihedral = y_m * sin(dihedral);
end    
end
