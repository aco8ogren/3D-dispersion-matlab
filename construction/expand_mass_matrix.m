clear; close all;

mass_matrix = readmatrix('mass_matrix.txt');

disp('original mass matrix')
disp(mass_matrix)

N_dof_per_node = 3;
N_node = 8;
N_dof = N_dof_per_node*N_node; % Total number of degrees of freedom in element

expanded_mass_matrix = zeros(N_dof,N_dof);

for i = 1:8
    for j = 1:8
        row_inds = ((i-1)*3+1):((i-1)*3+3);
        disp(row_inds)
        col_inds = ((j-1)*3+1):((j-1)*3+3);
        disp(col_inds)
        expanded_mass_matrix(row_inds,col_inds) = mass_matrix(i,j)*eye(N_dof_per_node,N_dof_per_node);
    end
end

disp('expanded mass matrix:')
disp(expanded_mass_matrix)

writematrix(expanded_mass_matrix,'expanded_mass_matrix.txt','Delimiter',' ')
