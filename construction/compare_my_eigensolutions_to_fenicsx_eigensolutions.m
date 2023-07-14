clear; close all;

my_data = load('ex_compute_modes_of_cube_solution.mat');

fenicsx_folder = '../../setting-up-fenicsx/';
fenicsx_fn = 'ex_eigenvalues_of_K_and_M_solution';
fenicsx_data = load([fenicsx_folder fenicsx_fn]);

disp('my eigenfrequencies')
disp(my_data.eigenfrequencies)

disp('fenicsx eigenfrequencies')
disp(fenicsx_data.eigenfrequencies)

figure
p(1) = plot(real(my_data.eigenfrequencies),'bo');
p(1).DisplayName = 'my eigenfrequencies';

hold on
p(2) = plot(fenicsx_data.eigenfrequencies,'k.');
p(2).DisplayName = 'fenicsx eigenfrequencies';

legend
title('my computed eigenfrequencies compared to fenicsx computed eigenfrequencies')

xlabel('eigenfrequency index')
ylabel('eigenfrequency value [Hz]')



