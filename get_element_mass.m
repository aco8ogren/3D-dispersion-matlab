function m_ele = get_element_mass(rho,element_size)
    % Stiffness matrix for linear hexahedral element
    
    element_size = element_size(1);

    m = rho*element_size^3; % total mass of each element
    m_ele = (1/216)*m*...
            [
            8 0 0 4 0 0 4 0 0 2 0 0 4 0 0 2 0 0 2 0 0 1 0 0 ...
            0 8 0 0 4 0 0 4 0 0 2 0 0 4 0 0 2 0 0 2 0 0 1 0 ...
            0 0 8 0 0 4 0 0 4 0 0 2 0 0 4 0 0 2 0 0 2 0 0 1 ...
            4 0 0 8 0 0 2 0 0 4 0 0 2 0 0 4 0 0 1 0 0 2 0 0 ...
            0 4 0 0 8 0 0 2 0 0 4 0 0 2 0 0 4 0 0 1 0 0 2 0 ...
            0 0 4 0 0 8 0 0 2 0 0 4 0 0 2 0 0 4 0 0 1 0 0 2 ...
            4 0 0 2 0 0 8 0 0 4 0 0 2 0 0 1 0 0 4 0 0 2 0 0 ...
            0 4 0 0 2 0 0 8 0 0 4 0 0 2 0 0 1 0 0 4 0 0 2 0 ...
            0 0 4 0 0 2 0 0 8 0 0 4 0 0 2 0 0 1 0 0 4 0 0 2 ...
            2 0 0 4 0 0 4 0 0 8 0 0 1 0 0 2 0 0 2 0 0 4 0 0 ...
            0 2 0 0 4 0 0 4 0 0 8 0 0 1 0 0 2 0 0 2 0 0 4 0 ...
            0 0 2 0 0 4 0 0 4 0 0 8 0 0 1 0 0 2 0 0 2 0 0 4 ...
            4 0 0 2 0 0 2 0 0 1 0 0 8 0 0 4 0 0 4 0 0 2 0 0 ...
            0 4 0 0 2 0 0 2 0 0 1 0 0 8 0 0 4 0 0 4 0 0 2 0 ...
            0 0 4 0 0 2 0 0 2 0 0 1 0 0 8 0 0 4 0 0 4 0 0 2 ...
            2 0 0 4 0 0 1 0 0 2 0 0 4 0 0 8 0 0 2 0 0 4 0 0 ...
            0 2 0 0 4 0 0 1 0 0 2 0 0 4 0 0 8 0 0 2 0 0 4 0 ...
            0 0 2 0 0 4 0 0 1 0 0 2 0 0 4 0 0 8 0 0 2 0 0 4 ...
            2 0 0 1 0 0 4 0 0 2 0 0 4 0 0 2 0 0 8 0 0 4 0 0 ...
            0 2 0 0 1 0 0 4 0 0 2 0 0 4 0 0 2 0 0 8 0 0 4 0 ...
            0 0 2 0 0 1 0 0 4 0 0 2 0 0 4 0 0 2 0 0 8 0 0 4 ...
            1 0 0 2 0 0 2 0 0 4 0 0 2 0 0 4 0 0 4 0 0 8 0 0 ...
            0 1 0 0 2 0 0 2 0 0 4 0 0 2 0 0 4 0 0 4 0 0 8 0 ...
            0 0 1 0 0 2 0 0 2 0 0 4 0 0 2 0 0 4 0 0 4 0 0 8 ...
            ];
end