function design = get_design2(design_params)
    design = zeros(design_params.N_pix);
    for prop_idx = 1:3
        design(:,:,:,prop_idx) = get_prop(design_params,prop_idx);
    end
end
