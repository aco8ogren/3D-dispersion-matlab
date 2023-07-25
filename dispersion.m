function [wv,fr,ev] = dispersion(const,wavevectors)

    N_dof_per_node = 3; % Linear hex element
    N_node = prod(const.N_ele.*const.N_pix); % Total number of nodes in the model *after* boundary conditions have been applied
    N_dof = N_node*N_dof_per_node; % Total number of degrees of freedom in the model *after* boundary conditions have been applied
    if const.N_eig>N_dof
        error('Number of requested eigenvalues is larger than the number of degrees of freedom in the finite element model.')
    end


    if ~const.isUseGPU
        % Perform all computations on the CPU
        disp('Using CPU.')

        fr = zeros(size(wavevectors,2),const.N_eig);
        if const.isSaveEigenvectors
            ev = zeros(N_dof,size(wavevectors,2),const.N_eig);
        else
            ev = [];
        end

        [K,M] = get_system_matrices(const);

        if const.isUseParallel
            if isfield(const,'numParallelWorkers')
                parforArg = const.numParallelWorkers;
                % disp(['Parallelizing with ' num2str(gcp().NumWorkers) ' workers.']);
            else
                parforArg = Inf;
            end
        else
            parforArg = 0;
        end
        % for k_idx = 1:size(wavevectors,1); if k_idx == 1; warning('parfor loop is commented out - performance will suffer!')'; end % USE THIS TO DEBUG
        parfor (k_idx = 1:size(wavevectors,1), parforArg); if k_idx == 1; disp(['Running in parallel.']);  end % USE THIS FOR PERFORMANCE
            wavevector = wavevectors(k_idx,:);
            T = get_transformation_matrix(wavevector,const);
            Kr = T'*K*T;
            Mr = T'*M*T;

            [eig_vecs,eig_vals] = eigs(Kr,Mr,const.N_eig,const.sigma_eig);
            [eig_vals,idxs] = sort(diag(eig_vals)); % Sort eigenvalues by magnitude
            eig_vecs = eig_vecs(:,idxs); % Sort eigenvectors in corresponding fashion


            if const.isSaveEigenvectors
                %             ev(:,k_idx,:) = (eig_vecs./(diag(eig_vecs'*Mr*eig_vecs)'))'; % normalize by mass matrix
                ev(:,k_idx,:) = (eig_vecs./vecnorm(eig_vecs,2,1)).*exp(-1i*angle(eig_vecs(1,:))); % normalize by p-norm, align complex angle
                %             ev(:,k_idx,:) = (eig_vecs./max(eig_vecs))'; % normalize by max
                %             ev(:,k_idx,:) = eig_vecs'; % don't normalize
            end


            fr(k_idx,:) = sqrt(real(eig_vals));
            fr(k_idx,:) = fr(k_idx,:)/(2*pi);
        end

    elseif const.isUseGPU
        % Perform eig computation on the GPU
        disp('Using GPU.')

        N_eig = N_dof; % Because eig spits out all eigenvalues of the matrix

        fr = gpuArray.zeros(size(wavevectors,2),N_eig);
        if const.isSaveEigenvectors
            ev = gpuArray.zeros(N_dof,size(wavevectors,2),N_eig);
        else
            ev = [];
        end

        [K,M] = get_system_matrices_GPU(const);

        if const.isUseParallel
            if isfield(const,'numParallelWorkers')
                parforArg = const.numParallelWorkers;
            else
                parforArg = Inf;
            end
        else
            parforArg = 0;
        end
        for k_idx = 1:size(wavevectors,1);
            wavevector = wavevectors(k_idx,:);
            T = get_transformation_matrix_GPU(wavevector,const);
            Kr = T'*K*T;
            Mr = T'*M*T;
            
            Kr = 0.5*(Kr + Kr');
            Mr = 0.5*(Mr + Mr');

            Kr = full(Kr);
            Mr = full(Mr);

            [eig_vecs,eig_vals] = eig(Kr,Mr);
            [eig_vals,idxs] = sort(diag(eig_vals)); % Sort eigenvalues by magnitude
            eig_vecs = eig_vecs(:,idxs); % Sort eigenvectors in corresponding fashion

            if const.isSaveEigenvectors
                %             ev(:,k_idx,:) = (eig_vecs./(diag(eig_vecs'*Mr*eig_vecs)'))'; % normalize by mass matrix
                ev(:,k_idx,:) = (eig_vecs./vecnorm(eig_vecs,2,1)).*exp(-1i*angle(eig_vecs(1,:))); % normalize by p-norm, align complex angle
                %             ev(:,k_idx,:) = (eig_vecs./max(eig_vecs))'; % normalize by max
                %             ev(:,k_idx,:) = eig_vecs'; % don't normalize
            end

            fr(k_idx,:) = sqrt(max(real(eig_vals),0));
            fr(k_idx,:) = fr(k_idx,:)/(2*pi);
        end

    end
    wv = wavevectors;
end


