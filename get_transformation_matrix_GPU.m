function [T,dTdwavevector] = get_transformation_matrix_GPU(wavevector,const)

    N_node = const.N_ele.*const.N_pix + 1;
    number_of_nodes_full = prod(N_node);
    N_node_reduced = N_node - 1;
    number_of_nodes_reduced = prod(N_node_reduced);

    N_dof_per_node = 3;

    % Master nodes across x = 0, agent nodes across x = a
    r = [const.unit_cell_size(1); 0; 0];
    xphase = exp(1i*dot(wavevector,r));
    if nargout == 2
        dxphasedwavevector = 1i*r*xphase;
    end
    
    % Master nodes across y = 0, agent nodes across y = a
    r = [0; const.unit_cell_size(2); 0];
    yphase = exp(1i*dot(wavevector,r));
    if nargout == 2
        dyphasedwavevector = 1i*r*yphase;
    end

    % Master nodes across z = 0, agent nodes across z = a
    r = [0; 0; const.unit_cell_size(3)];
    zphase = exp(1i*dot(wavevector,r));
    if nargout == 2
        dzphasedwavevector = 1i*r*zphase;
    end

    % Master nodes along (y,z) = (0,0), agent nodes along (y,z) = (a,a)
    r = [0; const.unit_cell_size(2); const.unit_cell_size(3)];
    yzphase = exp(1i*dot(wavevector,r));
    if nargout == 2
        dyzphasedwavevector = 1i*r*yzphase;
    end

    % Master nodes along (x,z) = (0,0), agent nodes along (x,z) = (a,a)
    r = [const.unit_cell_size(1); 0; const.unit_cell_size(3)];
    xzphase = exp(1i*dot(wavevector,r));
    if nargout == 2
        dxzphasedwavevector = 1i*r*xzphase;
    end

    % Master nodes along (x,y) = (0,0), agent nodes along (x,y) = (a,a)
    r = [const.unit_cell_size(1); const.unit_cell_size(2); 0];
    xyphase = exp(1i*dot(wavevector,r));
    if nargout == 2
        dxyphasedwavevector = 1i*r*xyphase;
    end
    
    % Master nodes at (x,y,z) = (0,0,0), agent nodes at (x,y,z) = (a,a,a)
    r = [const.unit_cell_size(1); const.unit_cell_size(2); const.unit_cell_size(3)];
    xyzphase = exp(1i*dot(wavevector,r));
    if nargout == 2
        dxyzphasedwavevector = 1i*r*xyzphase;
    end

    %% Define sets of nodes using the full (unreduced) global node numbering scheme
    % For this set, it is not okay to have preceeding set list a node that
    % is in one of the later sets, since each agent can only follow one
    % master.

    % Define the grid of unreduced node numbers
    unreduced_node_numbers_grid = reshape(gpuArray(1:number_of_nodes_full),N_node);
    
    % Define the nodes that are not agent nodes
    nonagent_unreduced_node_numbers_grid = unreduced_node_numbers_grid(1:(N_node(1)-1),1:(N_node(2)-1),1:(N_node(3)-1));

    % Define the nodes that are agents across x = a
    x_agent_node_numbers_grid = unreduced_node_numbers_grid(N_node(1),1:(N_node(2)-1),1:(N_node(3)-1));

    % Define the nodes that are agents across y = a
    y_agent_node_numbers_grid = unreduced_node_numbers_grid(1:(N_node(1)-1),N_node(2),1:(N_node(3)-1));

    % Define the nodes that are agents across z = a
    z_agent_node_numbers_grid = unreduced_node_numbers_grid(1:(N_node(1)-1),1:(N_node(2)-1),N_node(3));

    % Define the nodes that are agents along (y,z) = (a,a)
    yz_agent_node_numbers_grid = unreduced_node_numbers_grid(1:(N_node(1)-1),N_node(2),N_node(3));

    % Define the nodes that are agents along (x,z) = (a,a)
    xz_agent_node_numbers_grid = unreduced_node_numbers_grid(N_node(1),1:(N_node(2)-1),N_node(3));

    % Define the nodes that are agents along (x,y) = (a,a)
    xy_agent_node_numbers_grid = unreduced_node_numbers_grid(N_node(1),N_node(2),1:(N_node(3)-1));

    % Define the nodes that are agents at (x,y,z) = (a,a,a)
    xyz_agent_node_numbers_grid = unreduced_node_numbers_grid(N_node(1),N_node(2),N_node(3));
    
    % Concatenate into a list of all the agent nodes
    agent_node_numbers = [nonagent_unreduced_node_numbers_grid(:); x_agent_node_numbers_grid(:); y_agent_node_numbers_grid(:); z_agent_node_numbers_grid(:); yz_agent_node_numbers_grid(:); xz_agent_node_numbers_grid(:); xy_agent_node_numbers_grid(:); xyz_agent_node_numbers_grid(:)];

    %% Define sets of nodes using the reduced global node numbering scheme
    % For this set, it's okay to have a preceeding set list a node that is
    % in a following set, since each master can direct multiple agents
    
    % Define the grid of reduced node numbers
    reduced_node_numbers_grid = reshape(gpuArray(1:number_of_nodes_reduced),N_node_reduced);

    % Define the nodes that are not agent nodes (none of the nodes in the
    % reduced numbering scheme are agents, so we must label all of them)
    nonagent_reduced_node_numbers_grid = reduced_node_numbers_grid(1:N_node_reduced(1),1:N_node_reduced(2),1:N_node_reduced(3));
    
    % Define the nodes that are masters at x = 0
    x_master_node_numbers_grid = reduced_node_numbers_grid(1,1:N_node_reduced(2),1:N_node_reduced(3));

    % Define the nodes that are masters at y = 0
    y_master_node_numbers_grid = reduced_node_numbers_grid(1:N_node_reduced(1),1,1:N_node_reduced(3));

    % Define the nodes that are masters at z = 0
    z_master_node_numbers_grid = reduced_node_numbers_grid(1:N_node_reduced(1),1:N_node_reduced(2),1);
    
    % Define the nodes that are masters along (y,z) = (0,0)
    yz_master_node_numbers_grid = reduced_node_numbers_grid(1:N_node_reduced(1),1,1);

    % Define the nodes that are masters along (x,z) = (0,0)
    xz_master_node_numbers_grid = reduced_node_numbers_grid(1, 1:N_node_reduced(2), 1);

    % Define the nodes that are masters along (x,y) = (0,0)
    xy_master_node_numbers_grid = reduced_node_numbers_grid(1, 1, 1:N_node_reduced(3));

    % Define the master node that is at (x,y,z) = (0,0,0)
    xyz_master_node_numbers_grid = reduced_node_numbers_grid(1,1,1);
    
    % Concatenate into a list of all the master nodes
    master_node_numbers = [nonagent_reduced_node_numbers_grid(:); x_master_node_numbers_grid(:); y_master_node_numbers_grid(:); z_master_node_numbers_grid(:); yz_master_node_numbers_grid(:); xz_master_node_numbers_grid(:); xy_master_node_numbers_grid(:); xyz_master_node_numbers_grid];

    %% Map agents to masters
    % Map the agent nodes to their corresponding master nodes
    map_nodes_agent_to_master = [agent_node_numbers master_node_numbers];
    
    % Map the agent degrees of freedom to master degrees of freedom
    map_u_dof_agent_to_master = (map_nodes_agent_to_master - 1)*N_dof_per_node + 1; % This array only considers the u (x-dir) displacement
    map_dof_agent_to_master = repelem(map_u_dof_agent_to_master,3,1); % repelem the map_u_dof array in preparation for expansion from considering only u to considering u,v,w
    to_add = repmat([0; 1; 2],size(map_dof_agent_to_master,1)/3,2); % Create an array that can be added to the map_dof array to consider u,v,w
    map_dof_agent_to_master = map_dof_agent_to_master + to_add; % Add this to the degrees of freedom so that we consider u,v,w dofs, not just u

    % Define phase by which to map agent dofs to master dofs
    nonagent_phase_vector = gpuArray.ones(N_dof_per_node*numel(nonagent_unreduced_node_numbers_grid),1);
    xphase_vector = xphase*gpuArray.ones(N_dof_per_node*numel(x_agent_node_numbers_grid),1);
    yphase_vector = yphase*gpuArray.ones(N_dof_per_node*numel(y_agent_node_numbers_grid),1);
    zphase_vector = zphase*gpuArray.ones(N_dof_per_node*numel(z_agent_node_numbers_grid),1);
    yzphase_vector = yzphase*gpuArray.ones(N_dof_per_node*numel(yz_agent_node_numbers_grid),1);
    xzphase_vector = xzphase*gpuArray.ones(N_dof_per_node*numel(xz_agent_node_numbers_grid),1);
    xyphase_vector = xyphase*gpuArray.ones(N_dof_per_node*numel(xy_agent_node_numbers_grid),1);
    xyzphase_vector = xyzphase*gpuArray.ones(N_dof_per_node*numel(xyz_agent_node_numbers_grid),1);
    
    % Concatenate
    phase_vector = [nonagent_phase_vector; xphase_vector; yphase_vector; zphase_vector; yzphase_vector; xzphase_vector; xyphase_vector; xyzphase_vector];

    %% Construct the transformation matrix
    T = sparse(map_dof_agent_to_master(:,1),map_dof_agent_to_master(:,2),phase_vector);

    %% Construct the derivative of the transformation matrix with respect to wavevector if desired
    if nargout == 2
        dTdwavevector = cell(3,1); % Derivative with respect to each component of the wavevector
        % TODO complete this later - refer to get_transformation_matrix of
        % 2D-dispersion. The implementation shouldn't be that bad.
        % dphase_vector =  
    end
end