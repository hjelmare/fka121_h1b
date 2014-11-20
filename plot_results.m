clear all
clc

data = importdata('energy.data');

time = data(:,1);
energy = data(:,2);
potentialEnergy = data(:,3);
kineticEnergy = data(:,4);

data = importdata('pt.data');

temp = data(:,2);
pressure = data(:,3);

%%
hold on
plot(time, energy, 'b');
plot(time, potentialEnergy, 'g');
plot(time, kineticEnergy, 'r');
legend('total energy', 'potential energy', 'kinetic energy');
ylabel('energy [eV]');
xlabel('time (ps)');
%%
hold on
subplot(2,1,1);
plot(data(:,1), pressure, 'g');
legend('Pressure');
ylabel('Pressure [ENHET?]');
xlabel('time [ps]');

hold on

subplot(2,1,2);
plot(data(:,1), temp, 'm');
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
