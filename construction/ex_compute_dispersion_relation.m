clear; close all;

addpath('../') % Because this file lives in 3D-dispersion-matlab/construction

save_info.datetime_var = datetime;
save_info.script_start_time = replace(char(save_info.datetime_var),':','-');
save_info.mfilename_fullpath_var = mfilename('fullpath');
save_info.mfilename_var = mfilename;
save_info.output_folder = ['OUTPUT/output ' save_info.script_start_time];

t_total = tic;

isSaveOutput = false;

% c is the new const
c.N_pix = [10 10 10];
c.N_ele = [1 1 1];
c.N_eig = 20;
c.N_wv = [5 5 3];

c.E_min = 200e6; % 200e6
c.E_max = 200e6; % 200e9
c.rho_min = 1e3; % 8e2
c.rho_max = 1e3; % 8e3
c.nu_min = 0.3; % 0
c.nu_max = 0.3; % 0.49
c.unit_cell_size = [2 2 2]; % [m]
c.struct_idxs = [];

c.isUseGPU = false;
c.isSaveEigenvectors = false;
c.isUseParallel = true;
c.numParallelWorkers = 20;
c.sigma_eig = 1;

wavevectors = get_IBZ_wavevectors(c.N_wv,c.unit_cell_size,'none',1);
% wavevectors = [-pi/c.unit_cell_size(1) 0 0]; warning('''wavevectors'' is set to something unusual.')

% Generate design
design = repmat(ones(c.N_pix),1,1,1,3);

c.design = design;
c.design_scale = 'linear';

if ~isSaveOutput
    warning('isSaveOutput is set to false. Output will not be saved.')
end

[wv,fr,ev] = dispersion(c,wavevectors);

t_total = toc(t_total);
disp(['total time elapsed ' num2str(t_total) ' sec'])
