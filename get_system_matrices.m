function [global_stiffness_matrix,global_mass_matrix] = get_system_matrices(const)
    N_node_per_element = 8;
    N_dof_per_node = 3;
    N_dof_per_element = N_node_per_element*N_dof_per_node; % 24

    N_ele_global = const.N_pix.*const.N_ele;
    number_of_elements = prod(N_ele_global);
    N_node_global = N_ele_global + 1;
    number_of_nodes = prod(N_node_global);

    element_size = const.unit_cell_size./N_ele_global;
    
    elementwise_design = repelem(const.design,const.N_ele(1),const.N_ele(2),const.N_ele(3),1);

    if strcmp(const.design_scale,'linear')
        E = (const.E_min + elementwise_design(:,:,:,1).*(const.E_max - const.E_min));
        nu = (const.poisson_min + elementwise_design(:,:,:,3).*(const.poisson_max - const.poisson_min));
        rho = (const.rho_min + elementwise_design(:,:,:,2).*(const.rho_max - const.rho_min));
    elseif strcmp(const.design_scale,'log')
        E = exp(elementwise_design(:,:,:,1));
        nu = (const.poisson_min + elementwise_design(:,:,:,3).*(const.poisson_max - const.poisson_min));
        rho = exp(elementwise_design(:,:,:,2));
    else
        error('const.design_scale not recognized as log or linear')
    end
    
    global_node_numbers_grid = reshape(1:number_of_nodes,N_node_global);
    global_element_numbers_grid = reshape(1:number_of_elements,N_ele_global);
    smallest_global_degree_of_freedom_in_each_element = N_dof_per_node*...
                                                        (reshape(global_node_numbers_grid(1:end-1,1:end-1,1:end-1),[],1) - 1)...
                                                        + 1; % This also corresponds to the x-dir displacement of the node with the smallest x,y,z coords
    global_degrees_of_freedom_in_each_element = smallest_global_degree_of_freedom_in_each_element + [...
                                                                                                    (0:5) ...
                                                                                                    (0:5)+N_node_global(1)*N_dof_per_node ...
                                                                                                    (0:5)+N_node_global(1)*N_node_global(2)*N_dof_per_node ...
                                                                                                    (0:5)+N_node_global(1)*N_node_global(2)*N_dof_per_node+N_node_global(1)*N_dof_per_node ...
                                                                                                    ];
    
    % ---------------------------------------------------------------------
    % This can be used if the global dof ordering is different from the
    % local ordering that the local stiffness/mass matrices use

    % local_to_global_order = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
    % global_degrees_of_freedom_in_each_element = global_degrees_of_freedom_in_each_element(:,local_to_global_order)
    % ---------------------------------------------------------------------
    
    % Because the system matrices are symmetric, the row indices and column
    % indices could be interchanged without consequence.
    %
    % N_dof_per_element^2*number_of_elements cannot be interchanged with
    % number_of_dof = number_of_nodes*N_dof_per_node because of the fact
    % that elements share nodes. The sparse construction allows 
    % row_idx,col_idx pairs to be repeated, and adds their corresponding
    % value contributions together.
    row_idxs = reshape(repelem(global_degrees_of_freedom_in_each_element,1,N_dof_per_element)',N_dof_per_element^2*number_of_elements,1);
    col_idxs = reshape(repelem(global_degrees_of_freedom_in_each_element,N_dof_per_element,1)',N_dof_per_element^2*number_of_elements,1);
    
    % temp1 = row_idxs;
    % temp2 = col_idxs;
    % row_idxs = temp2;
    % col_idxs = temp1;

    % We assume that the design is given in the column major order. That
    % means that entries in the design array that has a low first index
    % correspond to a pixel whose coordinates have a low x value, low
    % second index --> pixel coordinates have a low y value, and low third
    % index --> pixel coordinates have a low z value.
    %
    % This column major ordering can be visualized by visualize_designs,
    % using the helper function visualize_design.
    stiffness_values = get_element_stiffness(E(:),nu(:),element_size); % Each row gives the local stiffness matrix entries for each element
    mass_values = get_element_mass(rho(:),element_size); % Each row gives the local mass matrix entries for each element

    stiffness_values = stiffness_values'; % Transpose so that stiffness_values(:) unpacks along the right direction, since matlab is column major
    mass_values = mass_values'; % Transpose so that mass_values(:) unpacks along the right direction, since matlab is column major

    global_stiffness_matrix = sparse(row_idxs,col_idxs,stiffness_values);
    global_mass_matrix = sparse(row_idxs,col_idxs,mass_values);

    % fig = figure;
    % ax = axes(figure);
    % scatter3(0,0,0)
    % hold on
    % scatter3(4,4,4)
    % for i = 1:3
    %     for j = 1:3
    %         for k = 1:3
    %             disp('----')
    %             disp(i)
    %             disp(j)
    %             disp(k)
    %             text(i,j,k,num2str(global_node_numbers_grid(i,j,k)))
    %             hold on
    %         end
    %     end
    % end
end