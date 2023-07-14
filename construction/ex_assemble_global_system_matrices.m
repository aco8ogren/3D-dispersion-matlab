clear; close all;

addpath('../')

const.N_pix = [2 2 2];
const.N_ele = [1 1 1];
const.unit_cell_size = [1 1 1];

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

[K,M] = get_system_matrices(const);

mfilename_var = mfilename;
save([mfilename_var '.mat'],'K','M')