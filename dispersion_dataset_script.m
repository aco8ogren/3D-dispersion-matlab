clear; close all;
save_info.datetime_var = datetime;
save_info.script_start_time = replace(char(save_info.datetime_var),':','-');
save_info.mfilename_fullpath_var = mfilename('fullpath');
save_info.mfilename_var = mfilename;
save_info.output_folder = ['OUTPUT/output ' save_info.script_start_time];

t_total = tic;

isSaveOutput = true;
c.isSaveEigenvectors = false;

checkpoint_chunk_size = 1;

% c is the new const
c.N_pix = [10 10 10];
c.N_ele = [1 1 1];
c.N_eig = 10;
c.N_struct = 10000;
c.N_wv = [11 11 6];
c.rng_seed_offset = 0;

c.E_min = 200e6; % 200e6
c.E_max = 200e9; % 200e9
c.rho_min = 1e3; % 8e2
c.rho_max = 1e4; % 8e3
c.nu_min = 0; % 0
c.nu_max = 0.45; % 0.49
c.unit_cell_size = [1 1 1]; % [m]
c.struct_idxs = [];
c.design_scale = 'linear';

c.isUseParallel = true;
c.isUseGPU = false;

c.sigma_eig = 1;

wavevectors = get_IBZ_wavevectors(c.N_wv,c.unit_cell_size,'none',1);
% wavevectors = [-pi/c.unit_cell_size(1) 0 0]; warning('''wavevectors'' is set to something unusual.')

EIGENVALUE_DATA = zeros(size(wavevectors,1),c.N_eig,c.N_struct);
WAVEVECTOR_DATA = wavevectors; % zeros(prod(N_wv,3));
CONSTITUTIVE_DATA = containers.Map; % TODO
DESIGN_DATA = zeros([c.N_pix 3 c.N_struct]);

% Generate design
design_params = design_parameters;
design_params.design_number = []; % leave empty
design_params.design_style = 'kernel';
design_params.design_options = struct('kernel','periodic','period',c.unit_cell_size,'sigma_f',1,'sigma_l',1,'symmetry_type','none','N_value',inf);
design_params.N_pix = c.N_pix;

if ~isSaveOutput
    warning('isSaveOutput is set to false. Output will not be saved.')
end

checkpoint_counter = 0;
wb_struct = waitbar(0,['Performing dispersion computation...' newline 'structure 0/' num2str(c.N_struct)],'Position',[585 504.7500 270 58.5000]); % Seems the position argument isn't working
wb_struct.Position = [585 504.7500 270 58.5000]; % So I guess I'll just include this redundant line
for struct_idx = 1:c.N_struct
    c.struct_idxs = [c.struct_idxs struct_idx];
    design_params.design_number = struct_idx + c.rng_seed_offset;
    design_params = design_params.prepare();
    design = get_design(design_params);
    c.design = design;

    DESIGN_DATA(:,:,:,:,struct_idx) = design;
    [wv,fr,ev] = dispersion(c,wavevectors);
    EIGENVALUE_DATA(:,:,struct_idx) = fr;
    waitbar(struct_idx/c.N_struct,wb_struct,['Performing dispersion computation...' newline 'structure ' num2str(struct_idx) '/' num2str(c.N_struct)])
    if isSaveOutput
        if struct_idx/checkpoint_chunk_size == floor(struct_idx/checkpoint_chunk_size) && struct_idx ~= c.N_struct
            % Set up save locations
            setup_output_folder(save_info)

            % Specify variables to save
            vars_to_save = {'WAVEVECTOR_DATA','EIGENVALUE_DATA','c','design_params'};

            % Specify filename
            save_info.checkpoint_output_file_path_curr = [save_info.output_folder '/checkpoint' num2str(checkpoint_counter) '.mat'];

            % Save the data
            save(save_info.checkpoint_output_file_path_curr,vars_to_save{:},'-v7.3');
            
            % Check for previous checkpoint, and delete it
            save_info.checkpoint_output_file_path_prev = [save_info.output_folder '/checkpoint' num2str(checkpoint_counter - 1) '.mat'];
            if isfile(save_info.checkpoint_output_file_path_prev)
                delete(save_info.checkpoint_output_file_path_prev)
            end

            checkpoint_counter = checkpoint_counter + 1;
        end
    end
end
close(wb_struct)

%% Save data
if isSaveOutput
    % Set up save locations
    setup_output_folder(save_info)
    
    % Specify variables to save
    vars_to_save = {'WAVEVECTOR_DATA','EIGENVALUE_DATA','c','design_params'};
    
    % Specify filename
    save_info.full_output_file_path = [save_info.output_folder '/DATA' ...
        ' N_pix' num2str(c.N_pix(1)) 'x' num2str(c.N_pix(2)) 'x' num2str(c.N_pix(3))...
        ' N_ele' num2str(c.N_ele(1)) 'x' num2str(c.N_ele(2)) 'x' num2str(c.N_ele(3))...
        ' N_wv' num2str(c.N_wv(1)) 'x' num2str(c.N_wv(2)) 'x' num2str(c.N_wv(3))...
        ' N_disp' num2str(c.N_struct)...
        ' N_eig' num2str(c.N_eig)...
        ' offset' num2str(c.rng_seed_offset) ' ' save_info.script_start_time '.mat'];
    
    % Save the data
    save(save_info.full_output_file_path,vars_to_save{:},'-v7.3');

    % Delete the existing checkpoint since the full data has been written
    delete(save_info.checkpoint_output_file_path_curr)
end

t_total = toc(t_total);
disp(['total time elapsed ' num2str(t_total) ' sec'])
