function wavevectors = get_IBZ_wavevectors(N_wv,unit_cell_size,symmetry_type,N_tesselations)
    a = unit_cell_size; % vector of length 3
    if ~exist('symmetry_type','var')
        symmetry_type = 'none';
    end
    if ~exist('num_tesselations','var')
        N_tesselations = 1;
    end
    if numel(N_wv) == 1
        N_wv = [N_wv N_wv];
    end
    
    switch symmetry_type
        case 'omit'
            [X,Y,Z] = ndgrid(linspace(-pi/a(1),pi/a(1),N_wv(1)),linspace(-pi/a(2),pi/a(2),N_wv(2)),linspace(-pi/a(3),pi/a(3),N_wv(3))); % a square centered at the origin of side length 2*pi/a
            gamma_x = X(true(size(X))); gamma_y = Y(true(size(Y))); gamma_z = Z(true(size(Z))); % rect
        case 'none'
            [X,Y,Z] = ndgrid(linspace(-pi/a(1),pi/a(1),N_wv(1)),linspace(-pi/a(2),pi/a(2),N_wv(2)),linspace(0,pi/a(3),N_wv(3))); % true asymmetric IBZ - note that this IBZ can be rotated arbitrarily!
            gamma_x = X(true(size(X))); gamma_y = Y(true(size(Y))); gamma_z = Z(true(size(Z)));% rect
        case 'p4mm'
            error('symmetry type not developed yet')
            [X,Y] = meshgrid(linspace(0,pi/a,N_wv(1)),linspace(0,pi/a,N_wv(2))); % these are more points than we need, so we chop it with triu
            gamma_x = X(triu(true(size(X)))); gamma_y = Y(triu(true(size(Y)))); % tri
        case 'c1m1'
            error('symmetry type not developed yet')
            if (floor(N_wv(2)/2) == N_wv(2)/2)
                error('For symmetry type c1m1, N_wv(2) must be an odd integer')
            end
            [X,Y] = meshgrid(linspace(0,pi/a,N_wv(1)),linspace(-pi/a,pi/a,N_wv(2))); % these are more points than we need, so we chop it with triu
            mask = triu(true([N_wv(1) (N_wv(2)+1)/2]));
            mask = [flipud(mask(2:end,:)); mask];
            gamma_x = X(mask); gamma_y = Y(mask);
        case 'p2mm'
            error('symmetry type not developed yet')
            [X,Y] = meshgrid(linspace(0,pi/a,N_wv(1)),linspace(0,pi/a,N_wv(2)));
            gamma_x = X(:); gamma_y = Y(:);
        otherwise
            error('symmetry_type not recognized')
    end
    
    wavevectors = N_tesselations.*cat(2,gamma_x,gamma_y,gamma_z);
end