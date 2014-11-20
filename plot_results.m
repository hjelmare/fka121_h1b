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

data = importdata('pt.data');

temp = data(:,2);
pressure = data(:,3);

%%
hold on
subplot(3,1,1)
plot(time1, energy1, 'b');
plot(time2, energy2, 'g');
plot(time3, energy3, 'r');
legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
title('Total energy')
ylabel('energy [eV]');
xlabel('time [ps]');
subplot(3,1,2)
plot(time1, potentialEnergy1, 'b');
plot(time2, potentialEnergy2, 'g');
plot(time3, potentialEnergy3, 'r');
legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
title('Potential energy');
ylabel('energy [eV]');
xlabel('time [ps]');
subplot(3,1,3)
plot(time1, kineticEnergy1, 'b');
plot(time2, kineticEnergy2, 'g');
plot(time3, kineticEnergy3, 'r');
legend('timestep = 0.1', 'timestep = 0.01', 'timestep = 0.001');
title('Kinetic energy');
ylabel('energy [eV]');
xlabel('time [ps]');
%%
hold on
subplot(2,1,1);
plot(data(:,1), pressure, 'g');
legend('Pressure');
ylabel('Pressure [eV/Å^{3}]');
xlabel('time [ps]');

hold on

subplot(2,1,2);
plot(data(:,1), temp, 'b');
legend('Temperature');
ylabel('energy [K]');
xlabel('time [ps]');

useFrom = 0.6;
start = fix(length(energy)*useFrom);

meanTemp = mean(temp(start:end))
meanPress = mean(pressure(start:end))

%%

clear all
clc
clf

data = importdata('velcor.data');

plot(data(:,1),data(:,2))

%%

clear all
clc

data = importdata('spectrum.data');

hold on
plot(data(:,1),data(:,2),'b')

%%

clear all
clc
clf

data = importdata('pt2.data');

hold on
plot(data(:,1),data(:,2),'r')
plot(data(:,1),data(:,3),'b')
hold off

%% Velocity correlation function

clear all
clc
clf

data = importdata('velcor.data');
plot(data(:,1),data(:,2));

ylabel('Velocity correlation [(Å/ps^2)^2]')
xlabel('\Delta t [ps]')

saveas(gcf,'velcor.png','png')

%% Spectrum analysis

clear all
clc
clf

data = importdata('spectrum.data');
plot(data(:,1),data(:,2));

ylabel('Spectrum of velocity correlation function [(Å/ps^2)^2]')
xlabel('\omega [rad]')

saveas(gcf,'spectrum.png','png')

%% nothing here to see



