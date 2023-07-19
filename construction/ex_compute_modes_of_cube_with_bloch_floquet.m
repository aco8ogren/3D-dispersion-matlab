clear; close all;

addpath('../')

const.N_pix = [2 2 2];
const.N_ele = [1 1 1];
const.unit_cell_size = [2 2 2];

temp = ones(1,prod(const.N_pix));
% temp = 1:prod(const.N_pix); temp = temp./max(temp);
prop = reshape(temp,const.N_pix);
const.design = repmat(prop,1,1,1,3);

const.E_min = 200e6;
const.E_max = 200e6;
const.rho_min = 1e3;
const.rho_max = 1e3;
const.poisson_min = 0.3;
const.poisson_max = 0.3;
const.design_scale = 'linear';

const.N_eig = 20;
const.sigma_eig = 1;

wavevector = [-pi/const.unit_cell_size(1) ,0,0];

[K,M] = get_system_matrices(const);
[T,~] = get_transformation_matrix(wavevector,const);

Kr = T'*K*T;
Mr = T'*M*T;

[vecs,vals] = eigs(Kr,Mr,const.N_eig,const.sigma_eig);
% [vecs,vals] = eigs(K,M,const.N_eig,const.sigma_eig,'IsSymmetricDefinite',true);

eigenvalues = diag(vals);
eigenfrequencies = sqrt(eigenvalues)/(2*pi);

figure
imagesc(imag(vals))
title('imaginary component of eigenvals')
colorbar
daspect([1 1 1])

figure
imagesc(real(vals))
title('real component of eigenvals')
colorbar
daspect([1 1 1])

mfn = mfilename;
vars_to_save = {'eigenvalues','eigenfrequencies'}
save([mfn '_solution.mat'],vars_to_save{:})
