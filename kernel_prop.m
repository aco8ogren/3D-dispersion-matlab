function prop = kernel_prop(kernel,N_pix,design_options)
    xx = linspace(0,1,N_pix(1)); yy = linspace(0,1,N_pix(2)); zz = linspace(0,1,N_pix(3));
    [X,Y,Z] = meshgrid(xx,yy,zz);
    points = [reshape(X,[],1),reshape(Y,[],1),reshape(Z,[],1)];
    %     scatter(points(:,1),points(:,2))

    switch kernel
        case 'matern52'
            error('not implemented')
            C = matern52(points,points,design_options.sigma_f,design_options.sigma_l);
        case 'periodic'
            C = periodic_kernel(points,points,design_options.sigma_f,design_options.sigma_l,design_options.period);
        otherwise
            error(['kernel name "' kernel '" not recognized'])
    end
    mu = 0.5*ones(1,size(points,1));
    prop = mvnrnd(mu, C);
    prop = reshape(prop,N_pix);
    prop(prop<0) = 0;
    prop(prop>1) = 1;
end