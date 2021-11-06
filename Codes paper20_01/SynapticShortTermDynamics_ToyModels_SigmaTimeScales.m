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
% Analysis of the the time constants sgma_{d}, sgma_{f} and sigma_{d+f}
%
% DA: Dayan & Abbott (2001)
%
% tau_stp, a_stp: generic notation for tau_{dep/fac}, a_{d/f}

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

SpkFreqin = 0:1:200;
SpkPerin = 1000./SpkFreqin;
Deltaspk = SpkPerin;


astp = 0.1;
taustp = 100;
Sgma_stp1 = Deltaspk./(Deltaspk./taustp-log(1-astp));
taustp = 200;
Sgma_stp2 = Deltaspk./(Deltaspk./taustp-log(1-astp));
taustp = 300;
Sgma_stp3 = Deltaspk./(Deltaspk./taustp-log(1-astp));

figure
hold on
plot(SpkFreqin,Sgma_stp1,'-b','linewidth',2);
plot(SpkFreqin,Sgma_stp2,'-r','linewidth',2);
plot(SpkFreqin,Sgma_stp3,'-g','linewidth',2);
set(gca,'fontsize',24);
xlabel('f_{spk}  [Hz]');
ylabel('\sigma_{dep/fac}  [ms]');
legend('\tau_{dep/fac}=100','\tau_{dep/fac}=200','\tau_{dep/fac}=300');




astp = 0.2;
taustp = 100;
Sgma_stp1 = Deltaspk./(Deltaspk./taustp-log(1-astp));
taustp = 200;
Sgma_stp2 = Deltaspk./(Deltaspk./taustp-log(1-astp));
taustp = 300;
Sgma_stp3 = Deltaspk./(Deltaspk./taustp-log(1-astp));

figure
hold on
plot(SpkFreqin,Sgma_stp1,'-b','linewidth',2);
plot(SpkFreqin,Sgma_stp2,'-r','linewidth',2);
plot(SpkFreqin,Sgma_stp3,'-g','linewidth',2);
set(gca,'fontsize',24);
xlabel('f_{spk}  [Hz]');
ylabel('\sigma_{dep/fac}  [ms]');
legend('\tau_{dep/fac}=100','\tau_{dep/fac}=200','\tau_{dep/fac}=300');


adep = 0.1;
afac = 0.1;
taudep = 100;
taufac = 100;
Sgma_depfac1 = Deltaspk./(Deltaspk*(1/taudep+1/taufac)-log(1-adep)-log(1-afac));
taudep = 200;
taufac = 200;
Sgma_depfac2 = Deltaspk./(Deltaspk*(1/taudep+1/taufac)-log(1-adep)-log(1-afac));
taudep = 300;
taufac = 300;
Sgma_depfac3 = Deltaspk./(Deltaspk*(1/taudep+1/taufac)-log(1-adep)-log(1-afac));


figure
hold on
plot(SpkFreqin,Sgma_depfac1,'-b','linewidth',2);
plot(SpkFreqin,Sgma_depfac2,'-r','linewidth',2);
plot(SpkFreqin,Sgma_depfac3,'-g','linewidth',2);
set(gca,'fontsize',24);
xlabel('f_{spk}  [Hz]');
ylabel('\sigma_{dep+fac}  [ms]');
legend('\tau_{dep}=\tau_{fac}=100','\tau_{dep}=\tau_{fac}=200','\tau_{dep}=\tau_{fac}=300');

adep = 0.1;
afac = 0.1;
taudep = 100;
taufac = 100;
Sgma_depfac1 = Deltaspk./(Deltaspk*(1/taudep+1/taufac)-log(1-adep)-log(1-afac));
taudep = 100;
taufac = 200;
Sgma_depfac2 = Deltaspk./(Deltaspk*(1/taudep+1/taufac)-log(1-adep)-log(1-afac));
taudep = 100;
taufac = 300;
Sgma_depfac3 = Deltaspk./(Deltaspk*(1/taudep+1/taufac)-log(1-adep)-log(1-afac));


figure
hold on
plot(SpkFreqin,Sgma_depfac1,'-b','linewidth',2);
plot(SpkFreqin,Sgma_depfac2,'-r','linewidth',2);
plot(SpkFreqin,Sgma_depfac3,'-g','linewidth',2);
set(gca,'fontsize',24);
xlabel('f_{spk}  [Hz]');
ylabel('\sigma_{dep+fac}  [ms]');
legend('\tau_{dep}=\tau_{fac}=100','\tau_{dep}=\tau_{fac}=200','\tau_{dep}=\tau_{fac}=300');











