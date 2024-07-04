function [fig,ax,scat] = visualize_dispersion_band(WAVEVECTOR_DATA,EIGENVALUE_DATA,eig_idx,struct_idx,ax_in)

    if ~exist('ax_in','var')
        fig_temp = figure;
        ax_temp = axes(fig_temp);
    else
        ax_temp = ax_in;
    end

    scat_temp = scatter3(WAVEVECTOR_DATA(:,1),WAVEVECTOR_DATA(:,2),WAVEVECTOR_DATA(:,3),[],EIGENVALUE_DATA(:,eig_idx,struct_idx),'filled');
    scat_temp.MarkerFaceAlpha = 0.3;
    daspect([1 1 1])

    if nargout > 0
        fig = fig_temp;
        if nargout > 1
            ax = ax_temp;
            if nargout > 2
                scat = scat_temp;
            end
        end
    end
end