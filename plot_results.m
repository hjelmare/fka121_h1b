%% the timesteps plot
clear all
clc

data = importdata('stort/energy_stort.data');

time1 = data(:,1);
energy1 = data(:,2);
potentialEnergy1 = data(:,3);
kineticEnergy1 = data(:,4);

data = importdata('mellan/energy_mellan.data');
time2 = data(:,1);
energy2 = data(:,2);
potentialEnergy2 = data(:,3);
kineticEnergy2 = data(:,4);

data = importdata('litet/energy_litet.data');
time3 = data(:,1);
energy3 = data(:,2);
potentialEnergy3 = data(:,3);
kineticEnergy3 = data(:,4);

%%

a=round(length(time1)/2); % divides the equilibrationsteps 
t1_1 = time1(1: a);
t1_2 = time1(a+1:end);

b=round(length(time2)/2); % divides the equilibrationsteps 
t2_1 = time2(1:b);
t2_2 = time2(b+1:end);


c=round(length(time3)/2); % divides the equilibrationsteps 
t3_1 = time3(1: c);
t3_2 = time3(c+1:end);
%%
%total energy
figure
textStorlek = 14;
legendStorlek = 11;
hold on
plot(t1_1, energy1(1:a), 'b.-');
plot(t2_1, energy2(1:b), 'g');
plot(t3_1, energy3(1:c), 'r--');

plot(t1_2, energy1(a+1:end), 'b.-')

plot(t2_2, energy2(b+1:end), 'g');

plot(t3_2, energy3(c+1:end), 'r--');
text =legend('timestep = 0.01', 'timestep = 0.001', 'timestep = 0.0001');
set(text, 'FontSize', legendStorlek);
title('Total energy', 'FontSize',textStorlek)
ylabel('energy [eV]', 'FontSize', textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);


% potential
figure
hold on
plot(t1_1, potentialEnergy1(1:a), 'b.-');
plot(t2_1, potentialEnergy2(1:b), 'g');
plot(t3_1, potentialEnergy3(1:c), 'r--');

plot(t1_2, potentialEnergy1(a+1:end), 'b.-');
plot(t2_2, potentialEnergy2(b+1:end), 'g');
plot(t3_2, potentialEnergy3(c+1:end), 'r--');
text = legend('timestep = 0.01', 'timestep = 0.001', 'timestep = 0.0001');
set(text, 'FontSize', legendStorlek);
title('Potential energy', 'FontSize',textStorlek);
ylabel('energy [eV]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);


% Kinetik
figure
hold on
plot(t1_1, kineticEnergy1(1:a), 'b.-');
plot(t2_1, kineticEnergy2(1:b), 'g');
plot(t3_1, kineticEnergy3(1:c), 'r--');

plot(t1_2, kineticEnergy1(a+1:end), 'b.-');
plot(t2_2, kineticEnergy2(b+1:end), 'g');
plot(t3_2, kineticEnergy3(c+1:end), 'r--');
text = legend('timestep = 0.01', 'timestep = 0.001', 'timestep = 0.0001');
set(text, 'FontSize', legendStorlek);
title('Kinetic energy', 'FontSize',textStorlek);
ylabel('energy [eV]', 'FontSize',textStorlek);
xlabel('time [ps]', 'FontSize',textStorlek);


%set(a, 'markers', 12)
%% Equilibration - T & P
clc
scaleLatPar = 0.1;

data = importdata('pt500.data');
pressure500 = data(:,2);
temp500 = data(:,3);
latpar500 = data(:,4);

time_1 = data(1:end/2,1);
time_2 = data(end/2:end,1);

latpar500_1 = latpar500(1:end/2);
latpar500_2 = latpar500(end/2:end);

press500_1 = pressure500(1:end/2);
press500_2 = pressure500(end/2:end);

data = importdata('pt700.data');
pressure700 = data(:,2);
temp700 = data(:,3);
latpar700 = data(:,4);

textStorlek = 14;
legendStorlek = 11;

subplot(2,1,1);
hold on
%plot(data(:,1), pressure500, 'r--');
plot(time_1, press500_1, 'r--');
plot(time_2, press500_2, 'g--');
%plot(data(:,1), (latpar500_1-4.05)*scaleLatPar,'b');
plot(time_1, (latpar500_1-4.05)*scaleLatPar,'r');
plot(time_2, (latpar500_2-4.05)*scaleLatPar,'g');
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

useFrom = 0.6;
start = fix(length(energy)*useFrom);

meanTemp = mean(temp(start:end))
meanPress = mean(pressure(start:end))

%% Velocity correlation function - production

clear all
clc
clf

textStorlek = 14;
legendStorlek = 11;

data_500 = importdata('velcor500.data');
data_700 = importdata('velcor700.data');
%data_900 = importdata('velcor900.data');

hold on
plot(data_500(:,1),data_500(:,2),'b')
plot(data_700(:,1),data_700(:,2),'r.-')
%plot(data_900(:,1),data_900(:,2),'r--')

ylabel('Velocity correlation [Å^2/ps^2]','FontSize',textStorlek)
xlabel('\Delta t [ps]','FontSize',textStorlek)

text = legend('T=500 ^\circC', 'T=700 ^\circC');
set(text, 'FontSize', legendStorlek);

saveas(gcf,'velcor.png','png')

%% Spectrum analysis - check axes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% and probably incorrect data as well.... fuck.

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
axis([0 0.5 -5 25])

saveas(gcf,'spectrum.png','png')

%% MSD function - production

clear all
clc
clf

data = importdata('msd500.data');
data2 = importdata('msd700.data');

textStorlek = 14;
legendStorlek = 11;

ds_line = 6*data2(:,1)*0.5485;    % to show MSD converged when averaging

hold on
plot(data(:,1),data(:,2), 'b');
plot(data2(:,1),data2(:,2),'r.-');
plot(data2(:,1),ds_line,'k--');

ylabel('MSD [Å^2]', 'FontSize',textStorlek)
xlabel('\Delta t [ps]', 'FontSize',textStorlek)
text = legend('T=500 ^\circC', 'T=700 ^\circC', '6 D_s \Delta t');
set(text, 'FontSize', legendStorlek);

saveas(gcf,'MSD.png','png')

%% Plot positions (3D-plot) of one particle - manual save - production
clear all
clc
clf

pos1 = importdata('position500.data');
pos2 = importdata('position700.data');


textStorlek = 14;
legendStorlek = 11;

hold on
plot3(pos1(:,2), pos1(:,3),pos1(:,4),'b');
<<<<<<< HEAD
plot3(pos2(:,2), pos2(:,3),pos2(:,4),'r--');
%plot3(pos3(:,2), pos3(:,3),pos3(:,4),'r.-');
text = legend('T=500 ^\circ C', 'T=973 ^\circ C', 'T=1173K');
=======
plot3(pos2(:,2), pos2(:,3),pos2(:,4),'r.-');
text = legend('T=500 ^\circC', 'T=700 ^\circC');
>>>>>>> 4f4830c2e6ed7fb90555c7c3f8db220074d9c24f
set(text, 'FontSize', legendStorlek);

xlabel('[Å]', 'FontSize',textStorlek);
ylabel('[Å]', 'FontSize',textStorlek);
zlabel('[Å]', 'FontSize',textStorlek);

