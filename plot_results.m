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


%%
hold on
plot(time, energy, 'b');
plot(time, potentialEnergy, 'g');
plot(time, kineticEnergy, 'r');
legend('total energy', 'potential energy', 'kinetic energy');
ylabel('energy [eV]');
xlabel('time (ps)');
%%
plot(data(:,1), pressure, 'c');
plot(data(:,1), temp, 'k');

legend('total energy', 'potential energy', 'kinetic energy', 'pressure', 'temp');
ylabel('energy [eV]');
xlabel('time (ps)');

%useFrom = 0.4;
%start = fix(length(energy)*useFrom);

<<<<<<< HEAD
%meanTemp = mean(temp(start:end))
%meanPress = mean(pressure(start:end))
=======
meanTemp = mean(temp(start:end))
meanPress = mean(pressure(start:end))

%%

clear all
clc
clf

data = importdata('velcor.data');

plot(data(:,1),data(:,2))
>>>>>>> 2ea529a6a51b442c974fbdd948fd05dcd4568c6f
