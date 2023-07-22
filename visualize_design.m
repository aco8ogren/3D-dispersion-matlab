function [fig,ax] = visualize_design(design, property_idx, layer_spacing, ax)
    % Function to visualize a design using cubes in a 3D plot

    % Set default values for property_idx, layer_spacing, and target axes
    if nargin < 2
        property_idx = 1;
    end
    if nargin < 3
        layer_spacing = 0;
    end
    if nargin < 4
        % Create a new figure and axes
        fig_handle = figure;
        ax_handle = axes(fig_handle);
    end

    % Get the size of the design
    sz = size(design);

    % Calculate the dimensions of each cube
    d = 1./sz(1:3);
    edges = d;

    % Alpha represents the transparency value of the cubes
    alpha = 1;

    % Flag to indicate whether to plot only the interior cubes
    isPlotInteriorOnly = true;

    % Iterate over the design elements
    for i1 = 1:sz(1)
        for i2 = 1:sz(2)
            for i3 = 1:sz(3)
                % Check if the cube should be plotted based on the boundary condition
                if isPlotInteriorOnly
                    if ~isOnBoundary(i1,i2,i3,sz)
                        continue
                    end
                end

                % Calculate the origin of the cube
                origin = [(i1 - 1)*d(1) (i2 - 1)*d(2) (i3 - 1)*(d(3) + layer_spacing)];

                % Get the color of the cube from the design array
                color = design(i1,i2,i3,property_idx);

                % Plot the cube
                plotcube(edges,origin,alpha,color)
            end
        end
    end

    % Set the aspect ratio of the plot to be equal
    daspect([1 1 1])

    % Check the number of output arguments and assign the figure and axes handles accordingly
    if nargout >= 1
        fig = fig_handle;
        ax = ax_handle;
    end

    % Nested function to check if a cube is on the boundary
    function out = isOnBoundary(i1,i2,i3,sz)
        indices = [i1 i2 i3];
        out = any(indices == 1) | any(indices == sz(1:3));
    end
end
