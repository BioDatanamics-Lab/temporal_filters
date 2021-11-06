% 2020-01-08

% Synaptic depression and facilitation models 
% (Dayan & Abbott, see Ermentrout & Terman book, Chapter 7) 

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightgray = [.7 .7 .7];

% Parameters

taudec = 5;
taudep = 400;
taufac = 50;
adep = 0.1;
afac = 0.2;
Xinf = 1;
Zinf = 0;

% Time definitions

Tmax = 5000;
dt = 0.01;
t = 0:dt:Tmax;

% Input spike definitions

SpkFreqin = 80;
SpkPerin = 1000/SpkFreqin;
Nspk = floor(Tmax/SpkPerin);
tspk = (1:Nspk)*SpkPerin;

% Synaptic dynamics (no short-term plasticity)

S = zeros(1,length(t));

jspk = 1;
for j=1:length(t)-1
    k1s = -S(j)/taudec;
    as = S(j)+k1s*dt;
    k2s = -as/taudec;
    S(j+1) = S(j)+(k1s+k2s)*dt/2;
    if t(j)> tspk(jspk)
        S(j+1) = S(j+1)+1;
        jspk = jspk+1;
    end
end

figure
hold on
plot(t,S,'-b','linewidth',2);
set(gca,'fontsize',24);
xlabel('t  [ms]');

% Synaptic short-term plasticity

x = zeros(1,length(t));
z = zeros(1,length(t));

x(1) = 1;
z(1) = 0;
jspk = 1;
for j=1:length(t)-1
    k1x = (Xinf-x(j))/taudep;
    k1z = (Zinf-z(j))/taufac;
    ax = x(j)+k1x*dt;
    az = z(j)+k1z*dt;
    k2x = (Xinf-ax)/taudep;
    k2z = (Zinf-az)/taufac;
    x(j+1) = x(j)+(k1x+k2x)*dt/2;
    z(j+1) = z(j)+(k1z+k2z)*dt/2;
    if t(j)>= tspk(jspk)
        x(j+1) = x(j+1)-adep*x(j);
        z(j+1) = z(j+1)+afac*(1-z(j+1));
        jspk = jspk+1;
    end
end

m(1) = 0;
for j=2:length(t)
    m(j) = x(j-1)*z(j);
end

% x-, z- and m-traces (m=DeltaS)

figure
hold on
plot(t,x,'-r','linewidth',2);
plot(t,z,'-g','linewidth',2);
plot(t,x.*z,'-b','linewidth',2);
axis([0 Tmax 0 1.2]);
set(gca,'fontsize',24);
xlabel('t  [ms]');

% x-, z- and m-traces (m=DeltaS) + peaks

figure
hold on
plot(t,x,'-r','linewidth',2);
plot(t,z,'-g','linewidth',2);
plot(t,x.*z,'-b','linewidth',2);
plot(t,m,'-','Color',lightblueish,'linewidth',1);
%plot(t,m,'-b','linewidth',2);
axis([0 Tmax 0 1.2]);
set(gca,'fontsize',24);
xlabel('t  [ms]');


% Sequence of peaks computed from the analytical solutions to the
% differential equations using the previous sequence values to generate the
% current ones (instead of the previous values in the differential
% equation discretization). 

Dtatspkin = SpkPerin;
Tspkmax = floor(Tmax/SpkPerin);
tspk = [SpkPerin:SpkPerin:Tspkmax*SpkPerin];


    
Pdep =  exp(-Dtatspkin/taudep);
Pfac =  exp(-Dtatspkin/taufac);
X = zeros(1,Tspkmax);
Z = zeros(1,Tspkmax);
Zaux = zeros(1,Tspkmax);
X(1) = 1;
Z(1) = afac;
for j=1:Tspkmax-1
    X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*Pdep;       
    Z(j+1) = (1-afac)*(Zinf+(Z(j)-Zinf)*Pfac)+afac;
end
M = X.*Z;
jspk = floor(tspk/dt);
Xfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
Zfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);

 
% Original version for Z-sequence
% X(1) = 1;
% Z(1) = 0;
% Zaux(1) = afac;
% for j=1:Tspkmax-1
%     X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*Pdep;   
%     Z(j+1) = Zinf+((1-afac)*Z(j)+afac-Zinf)*Pfac;   
%     Zaux(j+1) = Z(j+1)+afac*(1-Z(j+1));   
% end
% Zfp = (Zinf+(afac-Zinf)*Pfac)/(1-(1-afac)*Pfac);
% Maux = X.*Zaux;




plot(tspk,X,'o','Color',lightblueish,'linewidth',1);
plot(tspk,Z,'o','Color',lightblueish,'linewidth',1);
plot(tspk,M,'o','Color',lightblueish,'linewidth',1);



% Approximation of the Peak sequences by the toy model functions
% F = A+(1-A)*exp(-(t-tspk(1))/sigmadep)
% G = B*(1-C*exp(-(t-tspk(1))/sigmadep)

Fcut = 0.5*(X(1)-Xfp)+Xfp;
for k=2:length(tspk)
    if X(k)<Fcut
        kaux = k;    
        break
    end
end
sgmadep = (tspk(kaux)-tspk(1))/log((1-Xfp)/(X(kaux)-Xfp));
F = Xfp+(1-Xfp)*exp(-(tspk-tspk(1))/sgmadep);
E = sum((X-F).^2);

Gcut = 0.5*(Zfp-Z(1))+Z(1);
for k=2:length(tspk)
    if Z(k)>Gcut
        kaux = k;    
        break
    end
end

C = 1-afac/Zfp;
sgmafac = (tspk(kaux)-tspk(1))/log(Zfp*C/(Zfp-Z(kaux)));
G = Zfp*(1-C*exp(-(tspk-tspk(1))/sgmafac));
E = sum((Z-G).^2);
H = F.*G;
plot(tspk,F,'-','Color',lightcoral,'linewidth',2);
plot(tspk,G,'-','Color',lightblueish,'linewidth',2);
plot(tspk,H,'-','Color',lightblueish,'linewidth',2);
Hexpanded = Zfp*(Xfp-Xfp*C*exp(-(tspk-tspk(1))/sgmafac)+(1-Xfp)*exp(-(tspk-tspk(1))/sgmadep)-(1-Xfp)*C*exp(-(tspk-tspk(1))*(1/sgmadep+1/sgmafac)));
Hcut = Zfp*(Xfp-Xfp*C*exp(-(tspk-tspk(1))/sgmafac)+(1-Xfp)*exp(-(tspk-tspk(1))/sgmadep));
%plot(tspk,Hexpanded,'-k','linewidth',1);
%plot(tspk,Hcut,'-','Color',lightgray,'linewidth',1);



figure
hold on
plot(tspk,F,'-','Color',lightcoral,'linewidth',2);
plot(tspk,G,'-','Color',lightblueish,'linewidth',2);
plot(tspk,H,'-','Color',lightblueish,'linewidth',2);
plot(tspk,X,'or','linewidth',1);
plot(tspk,Z,'og','linewidth',1);
plot(tspk,M,'ob','linewidth',1);
axis([0 Tmax 0  1.1]);
set(gca,'fontsize',24);
xlabel('t_{spk}  [ms]');

figure
hold on
plot(tspk,F,'-r','linewidth',2);
plot(tspk,G,'-g','linewidth',2);
plot(tspk,H,'-b','linewidth',2);
plot(tspk,X,'or','linewidth',1);
plot(tspk,Z,'og','linewidth',1);
plot(tspk,M,'ob','linewidth',1);
axis([0 Tmax 0  1.1]);
set(gca,'fontsize',24);
xlabel('t_{spk}  [ms]');



% figure
% hold on
% % plot([1:Tspkmax],X,'or','linewidth',2);
% % plot([1:Tspkmax],Z,'og','linewidth',2);
% % plot([1:Tspkmax],M,'ob','linewidth',2);
% plot(tspk,X,'-or','linewidth',1);
% plot(tspk,Z,'-og','linewidth',1);
% plot(tspk,M,'-ob','linewidth',1);
% axis([0 Tmax 0 1.2]);
% set(gca,'fontsize',24);
% xlabel('t_{spk}  [ms]');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TOY = 0;
if TOY == 1

    % Toy-function model

    taudep=10;
    taufac = 1;
    B = 1;
    A = 0.2;
    C = 0.8;
    alpha=10;
    Ttmax = 10;
    tt = 0:0.01:Ttmax;
    f = A+(1-A)*exp(-tt/taudep);
    g = B*(1-C*exp(-tt/taufac));
    H = f.*g;


    % Rescaled toy-function model

    A = 0.2;
    C = 0.8;
    tau = 1;
    f = A+(1-A)*exp(-tt/tau);
    g = 1-C*exp(-tt);
    H = f.*g;

    figure
    hold on
    plot(tt,f,'-r','linewidth',2);
    plot(tt,g,'-g','linewidth',2);
    plot(tt,H,'-b','Color',lightblueish,'linewidth',2)
    plot(tt,H,'-b','linewidth',2)
    axis([0 Ttmax 0 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');

end
