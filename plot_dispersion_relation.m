function [fig,ax] = plot_dispersion_relation(wv,fr,eig_idx)

    if nargout > 0
        fig = figure;
        ax = axes(fig);
    end

    s = scatter3(wv(:,1),wv(:,2),wv(:,3),[],fr(:,eig_idx),'filled');
    s.MarkerFaceAlpha = 0.3;

    daspect([1 1 1])

    xlabel('wavevector x')
    ylabel('wavevector y')
    zlabel('wavevector z')
    colorbar
end