function setup_output_folder(save_info)
    if ~isfolder(save_info.output_folder)
        mkdir(save_info.output_folder);
        copyfile([save_info.mfilename_fullpath_var '.m'],[save_info.output_folder '/' save_info.mfilename_var '.m']);
    end
end