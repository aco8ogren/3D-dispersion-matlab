clear; close all;

addpath('../')

isSwapRowsAndCols = false;

fenicsx_folder = '../../setting-up-fenicsx/';
fenicsx_filename = 'ex_element_stiffness_and_mass.mat';
fenicsx_data = load([fenicsx_folder fenicsx_filename]);

assert(~any(imag(fenicsx_data.K),'all')) % Ensure there are no nonzero imaginary components in the Fenicsx stiffness matrix
assert(~any(imag(fenicsx_data.M),'all')) % Ensure there are no nonzero imaginary components in the Fenicsx mass matrix

% Construct stiffness and mass matrices from maple computations
E = 200e6;
nu = 0.3;
rho = 1000;
element_size = 1;

K_maple = get_element_stiffness(E,nu,element_size);
M_maple = get_element_mass(rho,element_size);

if isSwapRowsAndCols
    new_node_order = [1 2 4 3 5 6 8 7];
    new_dof_order = repelem(new_node_order,1,3)*3-2;
    to_add = repmat([0 1 2],1,8);
    new_dof_order = new_dof_order + to_add;
    fenicsx_data.K = fenicsx_data.K(new_dof_order,new_dof_order);
end

fig = figure;
tlo = tiledlayout(2,2);

% Plot stiffness matrix from fenicsx
nexttile
imagesc(real(fenicsx_data.K))
colorbar
daspect([1 1 1])
title('stiffness matrix - fenicsx')

% Plot stiffness matrix from maple
nexttile
imagesc(K_maple)
colorbar
daspect([1 1 1])
title('stiffness matrix - maple')

% Plot mass matrix from fenicsx
nexttile
imagesc(real(fenicsx_data.M))
colorbar
daspect([1 1 1])
title('mass matrix - fenicsx')

% Plot stiffness matrix from maple
nexttile
imagesc(M_maple)
colorbar
daspect([1 1 1])
title('mass matrix - maple')

display_matrix_comparison_statistics(fenicsx_data.K,K_maple,'fenicsx','maple')
display_matrix_comparison_statistics(fenicsx_data.M,M_maple,'fenicsx','maple')

function display_matrix_comparison_statistics(A, B, name_of_A, name_of_B)
    disp(['max of ' name_of_A ' matrix is']);
    disp(num2str(max(real(A), [], 'all'), 16));

    disp(['max of ' name_of_B ' matrix is']);
    disp(num2str(max(real(B), [], 'all'), 16));

    disp(['min of ' name_of_A ' matrix is']);
    disp(num2str(min(real(A), [], 'all'), 16));

    disp(['min of ' name_of_B ' matrix is']);
    disp(num2str(min(real(B), [], 'all'), 16));

    disp('max abs val difference between the two matrices is');
    disp(num2str(max(abs(real(A) - real(B)), [], 'all')));
end


% function display_matrix_comparison_statistics(A,B,name_of_A,name_of_B)
%     disp('max of fenicsx matrix is')
%     disp(num2str(max(real(fenicsx_data.K),[],'all'),16))
% 
%     disp('max of maple matrix is')
%     disp(num2str(max(K_maple,[],'all'),16))
% 
%     disp('min of fenicsx matrix is')
%     disp(num2str(min(real(fenicsx_data.K),[],'all'),16))
% 
%     disp('min of maple matrix is')
%     disp(num2str(min(K_maple,[],'all'),16))
% 
%     disp('max abs val difference between the two matrices is')
%     disp(num2str(max(abs(real(fenicsx_data.K)-K_maple),[],'all')))
% end