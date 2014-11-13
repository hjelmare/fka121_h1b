clear all
clc

data = importdata('energy.data');

time = data(:,1);
energy = data(:,2);
potentialEnergy = data(:,3);
kineticEnergy = data(:,4);

data = importdata('pt.data');

pressure = data(:,2);
temp = data(:,3);

clf

hold on
plot(time, energy, 'b');
plot(time, potentialEnergy, 'g');
plot(time, kineticEnergy*1e1, 'r');
plot(time(2:end), pressure*1e4, 'c');
plot(time(2:end), temp, 'k');

legend('total energy', 'potential energy', 'kinetic energy', 'pressure', 'temp');
ylabel('energy [eV]');
xlabel('time (ps)');

startTime = 40;

meanTemp = mean(temp(startTime:end))
meanPress = mean(pressure(startTime:end))