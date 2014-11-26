% in some places, data for several temperatures is required
% so you need to run the main program several times with
% different temperature settings and manually rename data files

clear all
clc

data = importdata('energyT1.data');

time1 = data(:,1);
energy1 = data(:,2);
potentialEnergy1 = data(:,3);
kineticEnergy1 = data(:,4);

data = importdata('energyT2.data');
time2 = data(:,1);
energy2 = data(:,2);
potentialEnergy2 = data(:,3);
kineticEnergy2 = data(:,4);

data = importdata('energyT3.data');
time3 = data(:,1);
energy3 = data(:,2);
potentialEnergy3 = data(:,3);
kineticEnergy3 = data(:,4);

%%
figure

textStorlek = 14;
legendStorlek = 11;
subplot(2,1,1)
hold on
plot(time1, energy1, 'b.-');
plot(time2, energy2, 'g');
plot(time3, energy3, 'r--');
text =legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Total energy', 'FontSize',textStorlek)
ylabel('energy [eV]', 'FontSize', textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

subplot(2,1,2)
hold on
plot(time2, energy2, 'g');
plot(time3, energy3, 'r--');
text =legend('timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Total energy', 'FontSize',textStorlek)
ylabel('energy [eV]', 'FontSize', textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

% potential
figure
subplot(2,1,1)
hold on
plot(time1, potentialEnergy1, 'b.-');
b = plot(time2, potentialEnergy2, 'g');
plot(time3, potentialEnergy3, 'r--');
text = legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Potential energy', 'FontSize',textStorlek);
ylabel('energy [eV]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

subplot(2,1,2)
hold on
plot(time2, potentialEnergy2, 'g');
plot(time3, potentialEnergy3, 'r--');
text =legend('timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Potential energy', 'FontSize',textStorlek)
ylabel('energy [eV]', 'FontSize', textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

% Kinetik
figure

subplot(2,1,1)
hold on
plot(time1, kineticEnergy1, 'b.-');
c = plot(time2, kineticEnergy2, 'g');
plot(time3, kineticEnergy3, 'r--');
text = legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Kinetic energy', 'FontSize',textStorlek);
ylabel('energy [eV]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

subplot(2,1,2)
hold on
plot(time2, kineticEnergy2, 'g');
plot(time3, kineticEnergy3, 'r--');
text =legend('timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Kinetic energy', 'FontSize',textStorlek)
ylabel('energy [eV]', 'FontSize', textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

%set(a, 'markers', 12)
%% Equilibration - T & P
data = importdata('pt500.data');
temp500 = data(:,2);
pressure500 = data(:,3);

data = importdata('pt700.data');
temp700 = data(:,2);
pressure700 = data(:,3);

textStorlek = 14;
legendStorlek = 11;


hold on
subplot(2,1,1);
plot(data(:,1), pressure500, 'g');
text  = legend('Pressure');
set(text, 'FontSize', legendStorlek);
ylabel('Pressure [eV/Å^{3}]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

hold on

subplot(2,1,2);
plot(data(:,1), temp500, 'b');
text = legend('Temperature');
set(text, 'FontSize', legendStorlek);
ylabel('temperature [K]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

figure
hold on
subplot(2,1,1);
plot(data(:,1), pressure700, 'g');
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

useFrom = 0.6;
start = fix(length(energy)*useFrom);

meanTemp = mean(temp(start:end))
meanPress = mean(pressure(start:end))

%% Velocity correlation function

clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;

data_500 = importdata('velcor_500.data');
data_700 = importdata('velcor_700.data');
data_900 = importdata('velcor_900.data');

hold on
plot(data_500(:,1),data_500(:,2),'b')
plot(data_700(:,1),data_700(:,2),'g.')
plot(data_900(:,1),data_900(:,2),'r--')

ylabel('Velocity correlation [(Å/ps^2)^2]','FontSize',textStorlek)
xlabel('\Delta t [ps]','FontSize',textStorlek)

text = legend('T = 773 K', 'T = 973 K', 'T = 1173 K');
set(text, 'FontSize', legendStorlek);

saveas(gcf,'velcor.png','png')

%% Spectrum analysis

clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;

data_500 = importdata('spectrum_500.data');
data_700 = importdata('spectrum_700.data');
data_900 = importdata('spectrum_900.data');

hold on
plot(data_500(:,1),data_500(:,2),'b')
plot(data_700(:,1),data_700(:,2),'g.')
plot(data_900(:,1),data_900(:,2),'r--')

ylabel('Spectrum of velocity correlation function [(Å/ps^2)^2]', ...
    'FontSize',textStorlek)
xlabel('\omega [rad/ps]','FontSize',textStorlek)

text = legend('T = 773 K', 'T = 973 K', 'T = 1173 K');
set(text, 'FontSize', legendStorlek);

saveas(gcf,'spectrum.png','png')

%% MSD function

clear all
clc
clf

%data = importdata('msd.data');
%data2 = importdata('msd2.data');
data3 = importdata('msd3.data');

textStorlek = 14;
legendStorlek = 11;

hold on
%plot(data(:,1),data(:,2), 'b');
%plot(data2(:,1),data2(:,2),'r--');
plot(data3(:,1),data3(:,2),'g.');

ylabel('MSD ', 'FontSize',textStorlek)
xlabel('\Delta t [ps]', 'FontSize',textStorlek)
text = legend('T=773K', 'T=973', 'T=1173K');
set(text, 'FontSize', legendStorlek);

saveas(gcf,'MSD.png','png')

%% Plot positions (3D-plot) of one particle
clear all
clc
pos1 = importdata('position.data');
%pos2 = importdata('position2.data');
%pos3 = importdata('position3.data');

textStorlek = 14;
legendStorlek = 11;

hold on
plot3(pos1(:,2), pos1(:,3),pos1(:,4),'b');
%plot3(pos2(:,2), pos2(:,3),pos2(:,4),'g.');
%plot3(pos3(:,2), pos3(:,3),pos3(:,4),'r.-');
text = legend('T=773K', 'T=973K', 'T=1173K');
set(text, 'FontSize', legendStorlek);

xlabel('[Å]', 'FontSize',textStorlek);
ylabel('[Å]', 'FontSize',textStorlek);
zlabel('[Å]', 'FontSize',textStorlek);

