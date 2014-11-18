clear all
clc

phiP = importdata('phiPressure.data');
phiT = importdata('phiTemp.data');
%%

plot(phiT(:,1), phiT(:,2), 'b');
hold on
plot(phiP(:,1), phiP(:,2), 'g');
legend('phi för temp', 'phi för pressure');
ylim([-1, 1])

plot(phiT(:,1), exp(-2),'r');