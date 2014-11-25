%% Equilibration - T & P

clear all
clc
clf


data = importdata('pt.data');
pressure = data(:,2);
temp = data(:,3);
latpar = data(:,4);

textStorlek = 14;
legendStorlek = 11;


subplot(2,1,1)
hold on
plot(data(:,1), temp, 'b');
plot(data(:,1), (latpar-4.05)*10000, 'r');
xlabel('time [ps]', 'FontSize',textStorlek);


subplot(2,1,2)
hold on
plot(data(:,1), pressure, 'g');
xlabel('time [ps]', 'FontSize',textStorlek);


useFrom = 0.6;
start = fix(length(pressure)*useFrom);

meanTemp = mean(temp(start:end))
meanPress = mean(pressure(start:end))

%% MSD function

clear all
clc
clf

data = importdata('msd.data');

textStorlek = 14;
legendStorlek = 11;

hold on
plot(data(:,1),data(:,2), 'b');

ylabel('MSD ', 'FontSize',textStorlek)
xlabel('\Delta t [ps]', 'FontSize',textStorlek)

%% Plot positions (3D-plot) of one particle
clear all
clf
clc
pos1 = importdata('position.data');

textStorlek = 14;
legendStorlek = 11;

hold on
plot3(pos1(:,2), pos1(:,3),pos1(:,4),'b');

xlabel('[Å]', 'FontSize',textStorlek);
ylabel('[Å]', 'FontSize',textStorlek);
zlabel('[Å]', 'FontSize',textStorlek);


%% Velocity correlation function

clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;

data = importdata('velcor.data');

hold on
plot(data(:,1),data(:,2),'b.')

ylabel('Velocity correlation [(Å/ps^2)^2]','FontSize',textStorlek)
xlabel('\Delta t [ps]','FontSize',textStorlek)


%% bla blah


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
textStorlek = 14;
legendStorlek = 11;
subplot(3,1,1)
hold on
plot(time1, energy1, 'b');
plot(time2, energy2, 'g.');
plot(time3, energy3, 'r--');
text =legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Total energy', 'FontSize',textStorlek)
ylabel('energy [eV]', 'FontSize', textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);
subplot(3,1,2)
hold on
plot(time1, potentialEnergy1, 'b');
plot(time2, potentialEnergy2, 'g.');
plot(time3, potentialEnergy3, 'r--');
text = legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Potential energy', 'FontSize',textStorlek);
ylabel('energy [eV]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);
subplot(3,1,3)
hold on
plot(time1, kineticEnergy1, 'b');
plot(time2, kineticEnergy2, 'g.');
plot(time3, kineticEnergy3, 'r--');
text = legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
set(text, 'FontSize', legendStorlek);
title('Kinetic energy', 'FontSize',textStorlek);
ylabel('energy [eV]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);

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

