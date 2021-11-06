% 2021-04-05

% DA model: evolution governed by the solution to the difference equations
% given by
% 
% Xn = Xss + Q(adep,taudep)^(n-1) (X1-Xss)
%
% Zn = Zss = Q(afac,taufac)^(n-1) (Z1-Zss)
%
% Q(astd,taustd) = (1-astd)*exp(-Deltaspk/taustd)
%
% Analysis of the dynamics of Q(astd,taustd)
%
% DA: Dayan & Abbott (2001)
%

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
mediumacquamarine = [0.4 0.8 0.6];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];


% Functions

Q=@(a,tau,Delta) (1-a).*exp(-Delta/tau);
Xbar=@(a,tau,Delta,Xinf) (1-exp(-Delta/tau))*Xinf./(1-(1-a).*exp(-Delta/tau));
Zbar=@(a,tau,Delta,Zinf) ((1-exp(-Delta/tau))*(1-a)*Zinf+a)./(1-(1-a).*exp(-Delta/tau));

% Parameters

% Reference parameter values:
%
% adep = 0.1;
% afac = 0.1;
Xinf = 1;
Zinf = 0;

astp = 0.1;
taustp = 100;

SpkFreqin = 1000;
SpkPerin = 1000/SpkFreqin;
Deltaspk = SpkPerin;


Nspk = 10000;
nvec = 1:1:Nspk;
Qn = Q(astp,taustp,Deltaspk).^(nvec-1);

% figure
% hold on
% plot(nvec,Qn,'o','Color',lightgray,'linewidth',2);
% axis([1 Nspk -0.1 1.1]);
% set(gca,'fontsize',24);
% xlabel('Spike #');

Xn = Xbar(astp,taustp,Deltaspk,Xinf)*(1-Qn)+Qn*1; 
Zn = Zbar(astp,taustp,Deltaspk,Zinf)*(1-Qn)+Qn*astp; 

figure
hold on
plot(-100,-100,'or','linewidth',2);
plot(-100,-100,'og','linewidth',2);
plot(-100,-100,'ob','linewidth',2);
plot(-100,-100,'o','Color',lightgray','linewidth',2);
plot([0 Nspk],[0 0],'--','Color',lightgray','linewidth',1);
plot([0 Nspk],[1 1],'--','Color',lightgray','linewidth',1);
plot(nvec,Qn,'o','Color',lightgray,'linewidth',1);
plot(nvec,Xn,'or','linewidth',2);
plot(nvec,Zn,'og','linewidth',2);
plot(nvec,Xn.*Zn,'ob','linewidth',2);
axis([1 80 -0.1 1.1]);
set(gca,'fontsize',24);
xlabel('Spike #');
legend('X_n','Z_n','X_n Z_n','Q_n');

tvec = nvec*Deltaspk;

figure
hold on
plot(-100,-100,'or','linewidth',2);
plot(-100,-100,'og','linewidth',2);
plot(-100,-100,'ob','linewidth',2);
plot(-100,-100,'o','Color',lightgray','linewidth',2);
plot([0 Nspk],[0 0],'--','Color',lightgray','linewidth',1);
plot([0 Nspk],[1 1],'--','Color',lightgray','linewidth',1);
plot(tvec,Qn,'o','Color',lightgray,'linewidth',1);
plot(tvec,Xn,'or','linewidth',2);
plot(tvec,Zn,'og','linewidth',2);
plot(tvec,Xn.*Zn,'ob','linewidth',2);
axis([0 Nspk*Deltaspk -0.1 1.1]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
legend('X_n','Z_n','X_n Z_n','Q_n');

figure
hold on
plot(-100,-100,'-r','linewidth',2);
plot(-100,-100,'-g','linewidth',2);
plot(-100,-100,'-b','linewidth',2);
plot(-100,-100,'o','Color',lightgray','linewidth',2);
plot([0 1000],[0 0],'--','Color',lightgray','linewidth',1);
plot([0 1000],[1 1],'--','Color',lightgray','linewidth',1);
plot(tvec,Xn,'-r','linewidth',2);
plot(tvec,Zn,'-g','linewidth',2);
plot(tvec,Xn.*Zn,'-b','linewidth',2);
axis([0 250 -0.1 1.2]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
legend('X_n','Z_n','\Delta S_n = X_n Z_n');




% Time scale for Q (time it takes for Q to decrease from Q=1 to Q = 0.37)
% sgma

alpha = 0.37;
sgma = log(alpha)/log(Q(astp,taustp,Deltaspk));


% Time scales for Xn and Zn (time it takes for Xn and  Zn to decrease from
% X1 and Z1 to X1*0.37 and Z1*0.37, respectively
% sgmad and sgmaf

adep = astp;
afac = astp;
taudep = taustp;
taufac = taustp;

sgmad = log(alpha)/log(Q(adep,taudep,Deltaspk));

% log((alpha*Xn(1)-Xbar(adep,taudep,Deltaspk,Xinf))/(Xn(1)-Xbar(adep,taudep,Deltaspk,Xinf)))/log(Q(adep,taudep,Deltaspk));
% sgmadf= log((alpha*Zn(1)-Zbar(afac,taufac,Deltaspk,Zinf))/(Zn(1)-Zbar(afac,taufac,Deltaspk,Zinf)))/log(Q(afac,taufac,Deltaspk));


% % Continuous version
% 
% x = 1:0.01:Nspk;
% Qx = Q(astp,taustp,Deltaspk).^(x-1);
% Qx1 = Qx.*log(Q(astp,taustp,Deltaspk));
% 
% plot(x,Qx,'-r','linewidth',2);
% plot(x,Qx1,'-g','linewidth',2);









