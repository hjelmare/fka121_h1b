clear all
clc
energy = importdata('totEnergy.data');
potentialEnergy = importdata('potentialEnergy.data');
kineticEnergy = importdata('kineticEnergy.data');
%%
plot(potentialEnergy(:,1), potentialEnergy(:,2), 'b');
hold on
plot(kineticEnergy(:,1), kineticEnergy(:,2), 'r');
plot(energy(:,1), energy(:,2), 'g');

legend('potential energy', 'kinetic energy', 'total energy');
ylabel('energy [eV]');
xlabel('time (ps)');