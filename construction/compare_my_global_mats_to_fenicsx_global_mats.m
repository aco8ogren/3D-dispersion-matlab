clear; close all;

% 1. Run setting-up-fenicsx/ex_get_global_stiffness_and_mass_matrices.py.
% Ensure it has the parameters set to properly define the mesh of interest.
% 
% 2. Open the vtu file output by the script in step 1. Its path should be
% setting-up-fenicsx/ex_get_global_stiffness_and_mass_matrices.vtu
%  
% 2. Copy and paste the coordinate in the vtu file into
% visualize_fenicsx_mesh_node_ordering.m and into the definition of
% fenicsx_nodal_coordinates below.
%
% 3. Run visualize_fenicsx_mesh_node_ordering.m to visualize the fenicsx node
% numbering scheme, and type their order into the definition of
% my_node_order below. Use the comment following that definition if
% confused on how to do this.
%
% 4. Run ex_assemble_global_system_matrices.m with parameters set to properly
% match the mesh of interest specified in the first step.

isSwapRowsAndCols = true; % For global dof ordering, I think this will always have to be true. I don't think fenicsx ever will have the same global dof ordering as me.

fenicsx_nodal_coordinates = [0 0 0 0.5 0 0 0 0.5 0 0.5 0.5 0 0 0 0.5 0.5 0 0.5 0 0.5 0.5 0.5 0.5 0.5 1 0 0 1 0.5 0 1 0 0.5 1 0.5 0.5 0 1 0 0.5 1 0 0 1 0.5 0.5 1 0.5 0 0 1 0.5 0 1 0 0.5 1 0.5 0.5 1 1 1 0 1 1 0.5 1 0 1 1 0.5 1 0 1 1 0.5 1 1 1 1 1];

fenicsx_nodal_coordinates = fenicsx_nodal_coordinates';
fenicsx_nodal_coordinates = reshape(fenicsx_nodal_coordinates,3,[])';

number_of_nodes = size(fenicsx_nodal_coordinates,1);

my_node_order = [1 2 9 3 4 10 13 14 21 5 6 11 7 8 12 15 16 22 17 18 23 19 20 24 25 26 27]; % my_node_order(i) gives the fenicsx node label that corresponds to the node labeled i in my ordering

% Load fenicsx system matrices
fenicsx_folder = '../../setting-up-fenicsx/';
fenicsx_fn = 'ex_get_global_stiffness_and_mass_matrices.mat';
fenicsx_data = load([fenicsx_folder fenicsx_fn]);

number_of_nodes = numel(my_node_order);
if isSwapRowsAndCols
    new_node_order = my_node_order;
    new_dof_order = repelem(new_node_order,1,3)*3-2;
    to_add = repmat([0 1 2],1,number_of_nodes);
    new_dof_order = new_dof_order + to_add;
    K_fenicsx = fenicsx_data.K(new_dof_order,new_dof_order);
    M_fenicsx = fenicsx_data.M(new_dof_order,new_dof_order);
end

% Load my system matrices
my_data = load('ex_assemble_global_system_matrices.mat');
K_my = my_data.K;
M_my = my_data.M;

% Plot for comparison
fig = figure;
tlo = tiledlayout(2,2);

nexttile
plot_matrix(K_my)
title('my global stiffness matrix')

nexttile
plot_matrix(K_fenicsx)
title('fenicsx global stiffness matrix')

nexttile
plot_matrix(M_my)
title('my global mass matrix')

nexttile
plot_matrix(M_fenicsx)
title('fenicsx global mass matrix')

function plot_matrix(A)
    imagesc(A)
    daspect([1 1 1])
    colorbar
end