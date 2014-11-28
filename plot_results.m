%% the timesteps plot
clear all
clc

data = importdata('energy_stort.data');

time1 = data(:,1);
energy1 = data(:,2);
potentialEnergy1 = data(:,3);
kineticEnergy1 = data(:,4);

data = importdata('energy_mellan.data');
time2 = data(:,1);
energy2 = data(:,2);
potentialEnergy2 = data(:,3);
kineticEnergy2 = data(:,4);

data = importdata('energy_litet.data');
time3 = data(:,1);
energy3 = data(:,2);
potentialEnergy3 = data(:,3);
kineticEnergy3 = data(:,4);


%total energy
figure
textStorlek = 14;
legendStorlek = 11;
hold on
plot(time1, energy1, 'b -.');
plot(time2, energy2, 'g');
plot(time3, energy3, 'r --');

text =legend('timestep = 0.01', 'timestep = 0.001' ...
    , 'timestep = 0.0001');
set(text, 'FontSize', legendStorlek);
title('Total energy', 'FontSize',textStorlek)
ylabel('energy [eV]', 'FontSize', textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);


% potential
figure
hold on
plot(time1, potentialEnergy1, 'b -.');
plot(time2, potentialEnergy2, 'g');
plot(time3, potentialEnergy3, 'r --');
text = legend('timestep = 0.01', 'timestep = 0.001', ...
    'timestep = 0.0001');
set(text, 'FontSize', legendStorlek);
title('Potential energy', 'FontSize',textStorlek);
ylabel('energy [eV]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);


% Kinetik
figure
hold on
plot(time1, kineticEnergy1, 'b -.');
plot(time2, kineticEnergy2, 'g');
plot(time3, kineticEnergy3, 'r--');
text = legend('timestep = 0.01', 'timestep = 0.001', ...
    'timestep = 0.0001');
set(text, 'FontSize', legendStorlek);
title('Kinetic energy', 'FontSize',textStorlek);
ylabel('energy [eV]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);


%set(a, 'markers', 12)
%% Equilibration - T & P
clc
scaleLatPar = 0.1;

data = importdata('pt500_eq.data');
pressure500 = data(:,2);
temp500 = data(:,3);
latpar500 = data(:,4);

time = data(:,1);


%%
data = importdata('pt700.data');
pressure700 = data(:,2);
temp700 = data(:,3);
latpar700 = data(:,4);

textStorlek = 14;
legendStorlek = 11;

subplot(2,1,1);
hold on
%plot(data(:,1), pressure500, 'r--');
plot(time, press500_1, 'r--');
plot(time, press500_2, 'g--');
%plot(data(:,1), (latpar500_1-4.05)*scaleLatPar,'b');
plot(time, (latpar500_1-4.05)*scaleLatPar,'r');
plot(time, (latpar500_2-4.05)*scaleLatPar,'g');
text  = legend('Pressure','Lattice parameter');
set(text, 'FontSize', legendStorlek);
ylabel('Pressure [eV/Å^{3}]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

hold on

subplot(2,1,2);
plot(data(:,1), temp500, 'b:');
text = legend('Temperature');
set(text, 'FontSize', legendStorlek);
ylabel('temperature [K]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

figure
subplot(2,1,1);
hold on
plot(data(:,1), pressure700, 'g');
plot(data(:,1), (latpar700-4.05)*scaleLatPar,'b');
text  = legend('Pressure');
set(text, 'FontSize', legendStorlek);
ylabel('Pressure [eV/Å^{3}]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

hold on

subplot(2,1,2);
plot(data(:,1), temp700, 'b');
text = legend('Temperature');
set(text, 'FontSize', legendStorlek);
ylabel('temperature [K]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);


%% Velocity correlation function - production

clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;

data_500 = importdata('velcor500_eq1.data');
data_700 = importdata('velcor700.data');

hold on
plot(data_500(:,1),data_500(:,2),'b')
plot(data_700(:,1),data_700(:,2),'r -.')

ylabel('Velocity correlation [Å^2/ps^2]','FontSize',textStorlek)
xlabel('\Delta t [ps]','FontSize',textStorlek)

text = legend('T=500 ^\circC', 'T=700 ^\circC');
set(text, 'FontSize', textStorlek);

saveas(gcf,'velcor.png','png')

%% Spectrum analysis
clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;

data_500 = importdata('spectrum500.data');
data_700 = importdata('spectrum700.data');

hold on
plot(data_500(:,1),data_500(:,2),'b')
plot(data_700(:,1),data_700(:,2),'r.-')

ylabel('Spectrum of velocity correlation function [Å^2/ps]', ...
    'FontSize',textStorlek)
xlabel('\omega [rad/ps]','FontSize',textStorlek)

text = legend('T=500 ^\circC', 'T=700 ^\circC');
set(text, 'FontSize', legendStorlek);
%axis([0 0.5 -5 25])

saveas(gcf,'spectrum.png','png')

%% MSD function - production

clear all
clc
clf

data = importdata('msd500_eq1.data');
data2 = importdata('msd700.data');

textStorlek = 14;
legendStorlek = 11;

% to show MSD converged when averaging
ds_line = 6*data2(:,1)*0.5422;    

hold on
plot(data(:,1),data(:,2), 'b');
plot(data2(:,1),data2(:,2),'r -.');
plot(data2(:,1),ds_line,'k--');

ylabel('MSD [Å^2]', 'FontSize',textStorlek)
xlabel('\Delta t [ps]', 'FontSize',textStorlek)
text = legend('T=500 ^\circC', 'T=700 ^\circC', ...
'6 D_s \Delta t');
set(text, 'FontSize', legendStorlek);

saveas(gcf,'MSD.png','png')

%% Plot positions (3D-plot) of one particle
clear all
clc
clf

pos1 = importdata('position500_eq1.data');
pos2 = importdata('position700.data');


textStorlek = 14;
legendStorlek = 11;

hold on
plot3(pos1(:,2), pos1(:,3),pos1(:,4),'b');

plot3(pos2(:,2), pos2(:,3),pos2(:,4),'r -.');
text = legend('T=500 ^\circC', 'T=700 ^\circC');

set(text, 'FontSize', legendStorlek);

xlabel('[Å]', 'FontSize',textStorlek);
ylabel('[Å]', 'FontSize',textStorlek);
zlabel('[Å]', 'FontSize',textStorlek);

