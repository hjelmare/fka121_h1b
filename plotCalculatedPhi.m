clear all
clc

phiP = importdata('phiPressure.data');
phiT = importdata('phiTemp.data');
%%

plot(phiT(:,1), phiT(:,2), 'b');
hold on
plot(phiP(:,1), phiP(:,2), 'g');
ylim([-1, 1])
xlabel('time [ps]')
ylabel('\phi (k)')
plot(phiT(:,1), exp(-2),'r');
legend('correlationfunction (\phi) - temperature', 'correlationfunction (\phi) - pressure', 'e^{-2}');

%% Plot positions (3D-plot) of one particle
clear all
clc
pos = importdata('position.data');

hold on
plot3(pos(:,2), pos(:,3),pos(:,4),'g');
xlabel('[Å]');
ylabel('[Å]');
zlabel('[Å]');
