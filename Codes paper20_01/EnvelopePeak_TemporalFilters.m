% 2020-01-23

% Synaptic depression and facilitation models 
% (Dayan & Abbott, see Ermentrout & Terman book, Chapter 7) 

% Inherited from "SynapticShortTermPlasticity_ToyModels.m"

% Envelope-peak curves (F, G and H) devoif from any biophsyical meaning
% with the specific goal of understanding how the productx of decreasing 
% (depression) and increasing (facilitation) curves interact to produce
% temporal filters.

% F = A + (1-A)*exp(-t/sgmad)
% G = B (1 - C*exp(-t/sgmaf)
% H = F*G

% B = 1 (scales the hight of the temporal filter H)
% Rescaling: \hat{t} = t/sgmaf 
% eta = sgmad/sgmaf
% hat{H} --> H

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightgray = [.7 .7 .7];

% Parameters

B = 1;
A = 0.2;
C = 0.1;
eta = 1;

Tmax = 25;
dt = 0.01;
t = 0:dt:Tmax;

F = A+(1-A)*exp(-t/eta);
G = B*(1-C*exp(-t));
H = F.*G;

figure
hold on
plot(t,F,'-r','linewidth',2);
plot(t,G,'-g','linewidth',2);
plot(t,H,'-b','linewidth',2);
%plot([0 Tmax],[A A],'--','Color',lightgray,'linewidth',1);
set(gca,'fontsize',24);
axis([0 Tmax 0 1.1]);
xlabel('t');
legend('F','G','H = F G','Location','SouthEast');





