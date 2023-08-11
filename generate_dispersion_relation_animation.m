clear; close all;

mfn = mfilename;

isWriteGif = false;
isSetBackgroundColor = true;
% BackgroundColor = [218,227,243]; % mylightblue
BackgroundColor = [0 0 0]; % Black

% data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\test set output 06-Jun-2023 19-30-12\checkpoint37.mat");
data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\IBZ output 18-May-2023 18-16-56\checkpoint233_ndgrid.mat");
% data = load("C:\Users\alex\OneDrive - California Institute of Technology\Documents\Graduate\Research\3D-dispersion-comsol\OUTPUT\RECYCLE\output 29-May-2023 17-01-51\checkpoint15 - Copy.mat");

fig = figure('color',BackgroundColor/255);
ax = axes(fig);

struct_idx = 1;
eig_idx = 10;
[fig,ax,scat] = visualize_dispersion_band(data.WAVEVECTOR_DATA,data.EIGENVALUE_DATA,eig_idx,struct_idx);
daspect([1 1 1])
set(ax,'Visible','off')
set(ax,'CameraViewAngleMode','Manual')
axis equal
zoom(0.5)

N_t = 60;
N_cycle = 1;
cycle_period = 6; % [seconds]
t = linspace(0,1*N_cycle,N_t*N_cycle);

if cycle_period/N_t<0.02
    warning('Matlab website warns that DelayTime lower than 0.02 may cause GIFs to play more slowly in many players and image viewers')
end

azimuth_center = 15;
azimuth_oscillation_amplitude = 25;
elevation_center = 25;
elevation_oscillation_amplitude = 15;

azimuths = azimuth_center + azimuth_oscillation_amplitude*sin(2*pi*t);

elevations = elevation_center + elevation_oscillation_amplitude*cos(2*pi*t);

for i = 1:length(t)
    az = azimuths(i);
    el = elevations(i);
    view(az,el)
    % pause(0.1)
    % F = getframe;
    % myframes{i} = F;
    % imwrite(F.cdata,'gifoutput.gif','DelayTime',1/N_t)
    if isWriteGif
        if i == 1
            gif([mfn '_struct' num2str(struct_idx) '_eig' num2str(eig_idx) '.gif'],'DelayTime',cycle_period/N_t,'frame',fig,'Overwrite',true)
        else
            gif
        end
    else
        pause(0.1)
    end
end