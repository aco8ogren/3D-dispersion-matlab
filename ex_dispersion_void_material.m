clear; close all;
% Compute dispersion relation for a single design

% Set the directory
cd('C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-matlab')

% Start the parallel pool
gcp;

isSaveOutput = false;
c.isSaveEigenvectors = false;

% Set dispersion problem parameters
% c is the new const
c.N_pix = [8 8 8];
c.N_ele = [1 1 1];
c.N_eig = 20;
c.N_struct = 1;
c.N_wv = [11 5 3];
c.rng_seed_offset = 0;

c.isAllowVoid = true;
c.E_min = 0; % 200e6
c.E_max = 200e9; % 200e9
c.rho_min = 0; % 8e2
c.rho_max = 1e4; % 8e3
c.nu_min = 0; % 0
c.nu_max = 0.45; % 0.49
c.unit_cell_size = [1 1 1]; % [m]
c.struct_idxs = [];
c.design_scale = 'linear';

c.isUseParallel = true;
c.isUseGPU = false;

c.sigma_eig = 1;

% wavevectors = get_IBZ_wavevectors(c.N_wv,c.unit_cell_size,'none',1);
% wavevectors = [-pi/c.unit_cell_size(1) 0 0]; warning('''wavevectors'' is set to something unusual.')
wavevectors = [linspace(0,pi/c.unit_cell_size(1),c.N_wv(1))' zeros(c.N_wv(1),1) zeros(c.N_wv(1),1)]; warning('''wavevectors'' is set to something unusual.');

% Generate a design
design = ones(c.N_pix);
idxs_1 = (floor(c.N_pix(1)/4) + 1):ceil(3*c.N_pix(1)/4);
idxs_2 = (floor(c.N_pix(1)/8) + 1):ceil(7*c.N_pix(1)/8);
idxs_3 = 1:c.N_pix(2);
[c1,c2,c3] = ndgrid(idxs_1,idxs_2,idxs_3);
linear_indices = sub2ind(size(design),c1(:),c2(:),c3(:));
design(linear_indices) = 0;

design = repmat(design,1,1,1,3);
c.design = design;

% Visualize the design
visualize_design(design,1,0.2)
title(['Unit cell design' newline 'Yellow indicates material, blue indicates void'])
xlabel('x')
ylabel('y')
zlabel('z')

% Give warning if isSaveOutput flag is false
if ~isSaveOutput
    warning('isSaveOutput is set to false. Output will not be saved.')
end

% Solve the dispersion problem
t_dispersion = tic; % Start the timer
[wv,fr,ev] = dispersion(c,wavevectors);
t_dispersion = toc(t_dispersion);
disp(['total time elapsed during dispersion computation: ' num2str(t_dispersion) ' sec'])

% Plot the resulting dispersion relation
fig = figure;
ax = axes(fig);
plot(wv(:,1),fr)
xlabel('wavevector [1/m]')
ylabel('frequency [1/s]')
title('dispersion relation')

ax.XLim = [min(wv(:,1)) max(wv(:,1))];

grid minor
