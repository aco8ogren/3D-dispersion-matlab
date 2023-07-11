clear; close all;

my_data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-matlab\construction\nodal_coordinates.mat");

my_nodal_coordinates = my_data.mat;

fenicsx_nodal_coordinates = [
    0 0 0;
    1 0 0;
    0 1 0;
    1 1 0;
    0 0 1;
    1 0 1;
    0 1 1;
    1 1 1];


fig = figure;
tlo = tiledlayout(1,2);

node_labels = cellfun(@(x) num2str(x),num2cell(1:8),'uniformoutput',false);
spacing = 0.1;

coords = {my_nodal_coordinates,fenicsx_nodal_coordinates};
titles = {'my node ordering','fenicsx node ordering'};

for i = 1:2
    ax = nexttile;
    scatter3(coords{i}(:,1),coords{i}(:,2),coords{i}(:,3),'k','filled')
    text(coords{i}(:,1)+spacing,coords{i}(:,2)+spacing,coords{i}(:,3)+spacing,node_labels)
    daspect([1 1 1])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(titles{i})
end
