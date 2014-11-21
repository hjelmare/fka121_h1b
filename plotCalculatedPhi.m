%% Plot phi
clear all
clc

phiP = importdata('phiPressure.data');
phiT = importdata('phiTemp.data');


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
pos1 = importdata('position1.data');
pos2 = importdata('position2.data');
pos3 = importdata('position3.data');
<<<<<<< HEAD
=======

>>>>>>> fc889f6b05cf692861d44aa2b72b54e584d117e4

textStorlek = 14;
legendStorlek = 11;

hold on
plot3(pos1(:,2), pos1(:,3),pos1(:,4),'b');
plot3(pos2(:,2), pos2(:,3),pos2(:,4),'g.');
plot3(pos3(:,2), pos3(:,3),pos3(:,4),'r.-');
text = legend('T=773K', 'T=973K', 'T=1173K');
set(text, 'FontSize', legendStorlek);

xlabel('[Å]', 'FontSize',textStorlek);
ylabel('[Å]', 'FontSize',textStorlek);
zlabel('[Å]', 'FontSize',textStorlek);
