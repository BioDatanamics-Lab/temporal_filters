% 2021-06-22

% % Synaptic depression and facilitation models 
% (Dayan & Abbott, see Ermentrout & Terman book, Chapter 7) 

% Extension of the standard synaptic depression-facilitation models to
% include two depression-facilitation processes interacting at the filter
% level. 

% Based on SynapticShortTermPlasticity_ToyModels.m (for the single
% depression-facilitation process)

% Cross-interactions model
%
% Delta S_n = ((1-etadep)*X_{1,n}+etadep*X_{2,n})*((1-etafac)*Z_{1,n}+etafac*Z_{2,n})

clearvars;
close all;

KSE = 2;
        % 1: Single set of depression-facilitation processes 
        %    It always shows this case
        % 2: Dependence of the filter time scales with the single-event
        % time scale tau_{dep,2}

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

taudec = 5;
adep = 0.1;
afac = 0.2;
Xinf = 1;
Zinf = 0;

taudep1 = 10;
taufac1 = 10;

taudep2 = 1000;
taufac2 = 1000;

etadep = 0.5;
etafac = 0.5;

% Time definitions

Tmax = 1000;
dt = 0.01;
t = 0:dt:Tmax;

% Input spike definitions

SpkFreqin = 100;
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
xlabel('S  [ms]');

% Synaptic short-term plasticity

x1 = zeros(1,length(t));
z1 = zeros(1,length(t));
x2 = zeros(1,length(t));
z2 = zeros(1,length(t));

x1(1) = 1;
z1(1) = 0;
x2(1) = 1;
z2(1) = 0;
jspk = 1;
for j=1:length(t)-1
    k1x1 = (Xinf-x1(j))/taudep1;
    k1z1 = (Zinf-z1(j))/taufac1;
    k1x2 = (Xinf-x2(j))/taudep2;
    k1z2 = (Zinf-z2(j))/taufac2;
    ax1 = x1(j)+k1x1*dt;
    az1 = z1(j)+k1z1*dt;
    ax2 = x2(j)+k1x2*dt;
    az2 = z2(j)+k1z2*dt;
    k2x1 = (Xinf-ax1)/taudep1;
    k2z1 = (Zinf-az1)/taufac1;
    k2x2 = (Xinf-ax2)/taudep2;
    k2z2 = (Zinf-az2)/taufac2;
    x1(j+1) = x1(j)+(k1x1+k2x1)*dt/2;
    z1(j+1) = z1(j)+(k1z1+k2z1)*dt/2;
    x2(j+1) = x2(j)+(k1x2+k2x2)*dt/2;
    z2(j+1) = z2(j)+(k1z2+k2z2)*dt/2;
    if t(j)>= tspk(jspk)
        x1(j+1) = x1(j+1)-adep*x1(j);
        z1(j+1) = z1(j+1)+afac*(1-z1(j+1));
        x2(j+1) = x2(j+1)-adep*x2(j);
        z2(j+1) = z2(j+1)+afac*(1-z2(j+1));
        jspk = jspk+1;
    end
end

m1 = zeros(1,length(t));
m2 = zeros(1,length(t));
mx = zeros(1,length(t));
for j=2:length(t)
    mx(j) = ((1-etadep)*x1(j-1)+etadep*x2(j-1))*((1-etafac)*z1(j)+etafac*z2(j));
    m1(j) = x1(j-1)*z1(j);
    m2(j) = x2(j-1)*z2(j);
end

figure
hold on
plot(t,x1,'-r','linewidth',2);
plot(t,z1,'-g','linewidth',2);
plot(t,x1.*z1,'-b','linewidth',2);
plot(t,m1,'-','Color',lightblueish,'linewidth',1);
plot(t,x2,'-r','linewidth',2);
plot(t,z2,'-g','linewidth',2);
plot(t,x2.*z2,'-b','linewidth',2);
plot(t,m2,'-','Color',lightblueish,'linewidth',1);
axis([0 Tmax 0 2.4]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
legend('x','z','x z','M');



% Sequence of peaks computed from the analytical solutions to the
% % differential equations using the previous sequence values to generate the
% current ones (instead of the previous values in the differential
% equation discretization). 


Dtatspkin = SpkPerin;
Tspkmax = floor(Tmax/SpkPerin);
tspk = SpkPerin:SpkPerin:Tspkmax*SpkPerin;
jspk = floor(tspk/dt);

Pdep1 =  exp(-Dtatspkin/taudep1);
Pfac1 =  exp(-Dtatspkin/taufac1);
Pdep2 =  exp(-Dtatspkin/taudep2);
Pfac2 =  exp(-Dtatspkin/taufac2);
X1 = zeros(1,Tspkmax);
Z1 = zeros(1,Tspkmax);
Zaux1 = zeros(1,Tspkmax);
X2 = zeros(1,Tspkmax);
Z2 = zeros(1,Tspkmax);
Zaux2 = zeros(1,Tspkmax);
X1(1) = 1;
Z1(1) = afac;
X2(1) = 1;
Z2(1) = afac;
for j=1:Tspkmax-1
    X1(j+1) = Xinf+((1-adep)*X1(j)-Xinf)*Pdep1;       
    Z1(j+1) = (1-afac)*(Zinf+(Z1(j)-Zinf)*Pfac1)+afac;
    X2(j+1) = Xinf+((1-adep)*X2(j)-Xinf)*Pdep2;       
    Z2(j+1) = (1-afac)*(Zinf+(Z2(j)-Zinf)*Pfac2)+afac;
end
X1fp = (1-Pdep1)*Xinf/(1-(1-adep)*Pdep1);
Z1fp = ((1-Pfac1)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac1);
X2fp = (1-Pdep2)*Xinf/(1-(1-adep)*Pdep2);
Z2fp = ((1-Pfac2)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac2);


M1 = X1.*Z1;
M2 = X2.*Z2;
Mcross_dep = ((1-etadep)*X1+etadep*X2);
Mcross_fac = ((1-etafac)*Z1+etafac*Z2);
Mcross = Mcross_dep.*Mcross_fac;

plot(tspk,X1,'o','Color',lightblueish,'linewidth',1);
plot(tspk,Z1,'o','Color',lightblueish,'linewidth',1);
plot(tspk,M1,'o','Color',lightblueish,'linewidth',1);
plot(tspk,X2,'s','Color',lightblueish,'linewidth',1);
plot(tspk,Z2,'s','Color',lightblueish,'linewidth',1);
plot(tspk,M2,'s','Color',lightblueish,'linewidth',1);
legend('x','z','x z','M');

figure
hold on
plot(-100,-100,'or','linewidth',2);
plot(-100,-100,'og','linewidth',2);
plot(-100,-100,'ob','linewidth',2);
plot(-100,-100,'o','Color',lightcoral,'linewidth',2);
plot(-100,-100,'o','Color',mediumacquamarine,'linewidth',2);
plot(-100,-100,'o','Color',lightblueish,'linewidth',2);
plot(-100,-100,'o','Color',darkgray,'linewidth',2);
plot(tspk,X1,'or','linewidth',1);
plot(tspk,Z1,'og','linewidth',1);
plot(tspk,M1,'ob','linewidth',1);
plot(tspk,X2,'o','Color',lightcoral,'linewidth',1);
plot(tspk,Z2,'o','Color',mediumacquamarine,'linewidth',1);
plot(tspk,M2,'o','Color',lightblueish,'linewidth',1);
plot(tspk,Mcross,'o','Color',darkgray,'linewidth',1);
axis([0 Tmax 0  1.2]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
%legend('X_n','Z_n','\DeltaS_n');
legend('X_{1,n}','Z_{1,n}','\DeltaS_{1,n}','X_{2,n}','Z_{2,n}','\DeltaS_{2,n}','\DeltaS_n^#');


% Approximation of the Peak sequences by the toy model functions
% F = A+(1-A)*exp(-(t-tspk(1))/sigmadep)
% G = B*(1-C*exp(-(t-tspk(1))/sigmadep)

Mcross_depfp = ((1-etadep)*X1fp+etadep*X2fp);
Mcross_facfp = ((1-etafac)*Z1fp+etafac*Z2fp);

% X1, X2, Mcross_dep 

Fcut1 = 0.5*(X1(1)-X1fp)+X1fp;
for k=2:length(tspk)
    if X1(k)<Fcut1
        kaux = k;    
        break
    end
end
sgmadep1 = (tspk(kaux)-tspk(1))/log((1-X1fp)/(X1(kaux)-X1fp));
F1 = X1fp+(1-X1fp)*exp(-(tspk-tspk(1))/sgmadep1);
E_F1 = sum((X1-F1).^2);

Fcut2 = 0.5*(X2(1)-X2fp)+X2fp;
for k=2:length(tspk)
    if X2(k)<Fcut2
        kaux = k;    
        break
    end
end
sgmadep2 = (tspk(kaux)-tspk(1))/log((1-X2fp)/(X2(kaux)-X2fp));
F2 = X2fp+(1-X2fp)*exp(-(tspk-tspk(1))/sgmadep2);
E_F2 = sum((X2-F2).^2);


Fcut_cross = 0.5*(Mcross_dep(1)-Mcross_depfp)+Mcross_depfp;
for k=2:length(tspk)
    if Mcross_dep(k)<Fcut_cross
        kaux=k;
        break;
    end
end
sgmadepcross = (tspk(kaux)-tspk(1))/log((1-Mcross_depfp)/(Mcross_dep(kaux)-Mcross_depfp));
Fcross = Mcross_depfp+(1-Mcross_depfp)*exp(-(tspk-tspk(1))/sgmadepcross);
E_Fcross = sum((Mcross_dep-Fcross).^2);


% Z1, Z2, Mcross_fac 

Gcut1 = 0.5*(Z1fp-Z1(1))+Z1(1);
for k=2:length(tspk)
    if Z1(k)>Gcut1
        kaux = k;    
        break
    end
end
C1 = 1-afac/Z1fp;
sgmafac1 = (tspk(kaux)-tspk(1))/log(Z1fp*C1/(Z1fp-Z1(kaux)));
G1 = Z1fp*(1-C1*exp(-(tspk-tspk(1))/sgmafac1));
E_G1 = sum((Z1-G1).^2);

Gcut2 = 0.5*(Z2fp-Z2(1))+Z2(1);
for k=2:length(tspk)
    if Z2(k)>Gcut2
        kaux = k;    
        break
    end
end
C2 = 1-afac/Z2fp;
sgmafac2 = (tspk(kaux)-tspk(1))/log(Z2fp*C2/(Z2fp-Z2(kaux)));
G2 = Z2fp*(1-C2*exp(-(tspk-tspk(1))/sgmafac2));
E_G2 = sum((Z2-G2).^2);

Gcut_cross = 0.5*(Mcross_facfp-Mcross_fac(1))+Mcross_fac(1);
for k=2:length(tspk)
    if Mcross_fac(k)>Gcut_cross
        kaux = k;    
        break
    end
end
Ccross = 1-afac/Mcross_facfp;
sgmafaccross = (tspk(kaux)-tspk(1))/log(Mcross_facfp*Ccross/(Mcross_facfp-Mcross_fac(kaux)));
Gcross = Mcross_facfp*(1-Ccross*exp(-(tspk-tspk(1))/sgmafaccross));
E_Gcross = sum((Mcross_fac-Gcross).^2);

FGcross = ((1-etadep)*F1+etadep*F2).*((1-etafac)*G1+etafac*G2);

[E_F1 E_F2 E_Fcross E_G1 E_G2 E_Gcross]

figure
hold on
plot(-100,-100,'or','linewidth',2);
%plot(-100,-100,'og','linewidth',2);
%plot(-100,-100,'ob','linewidth',2);
plot(-100,-100,'o','Color',lightcoral,'linewidth',2);
%plot(-100,-100,'o','Color',mediumacquamarine,'linewidth',2);
%plot(-100,-100,'o','Color',lightblueish,'linewidth',2);
plot(-100,-100,'o','Color',lightgray,'linewidth',2);
plot(tspk,X1,'or','linewidth',1);
%plot(tspk,Z1,'og','linewidth',1);
%plot(tspk,M1,'ob','linewidth',1);
plot(tspk,X2,'o','Color',lightcoral,'linewidth',1);
%plot(tspk,Z2,'o','Color',mediumacquamarine,'linewidth',1);
%plot(tspk,M2,'o','Color',lightblueish,'linewidth',1);
plot(tspk,Mcross_dep,'o','Color',lightgray,'linewidth',1);
axis([0 Tmax 0  1.2]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
plot(tspk,F1,'-k','linewidth',1);
plot(tspk,F2,'-k','linewidth',1);
plot(tspk,Fcross,'-k','linewidth',1);
legend('X_{1,n}','X_{2,n}','\DeltaS_{dep,n}');

figure
hold on
%plot(-100,-100,'or','linewidth',2);
plot(-100,-100,'og','linewidth',2);
%plot(-100,-100,'ob','linewidth',2);
%plot(-100,-100,'o','Color',lightcoral,'linewidth',2);
plot(-100,-100,'o','Color',mediumacquamarine,'linewidth',2);
%plot(-100,-100,'o','Color',lightblueish,'linewidth',2);
plot(-100,-100,'o','Color',lightgray,'linewidth',2);
%plot(tspk,X1,'or','linewidth',1);
plot(tspk,Z1,'og','linewidth',1);
%plot(tspk,M1,'ob','linewidth',1);
%plot(tspk,X2,'o','Color',lightcoral,'linewidth',1);
plot(tspk,Z2,'o','Color',mediumacquamarine,'linewidth',1);
%plot(tspk,M2,'o','Color',lightblueish,'linewidth',1);
plot(tspk,Mcross_fac,'o','Color',lightgray,'linewidth',1);
axis([0 Tmax 0  1.2]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
plot(tspk,G1,'-k','linewidth',1);
plot(tspk,G2,'-k','linewidth',1);
plot(tspk,Gcross,'-k','linewidth',1);
legend('Z_{1,n}','Z_{2,n}','\DeltaS_{fac,n}');

figure
hold on
%plot(-100,-100,'or','linewidth',2);
%plot(-100,-100,'og','linewidth',2);
plot(-100,-100,'ob','linewidth',2);
%plot(-100,-100,'o','Color',lightcoral,'linewidth',2);
%plot(-100,-100,'o','Color',mediumacquamarine,'linewidth',2);
plot(-100,-100,'o','Color',lightblueish,'linewidth',2);
plot(-100,-100,'o','Color',lightgray,'linewidth',2);
%plot(tspk,X1,'or','linewidth',1);
%plot(tspk,Z1,'og','linewidth',1);
plot(tspk,M1,'ob','linewidth',1);
%plot(tspk,X2,'o','Color',lightcoral,'linewidth',1);
%plot(tspk,Z2,'o','Color',mediumacquamarine,'linewidth',1);
plot(tspk,M2,'o','Color',lightblueish,'linewidth',1);
plot(tspk,Mcross,'o','Color',lightgray,'linewidth',1);
axis([0 Tmax 0  1.2]);
set(gca,'fontsize',24);
xlabel('t  [ms]');
plot(tspk,F1.*G1,'-k','linewidth',1);
plot(tspk,F2.*G2,'-k','linewidth',1);
plot(tspk,FGcross,'-k','linewidth',1);
legend('\Delta S_{1,n} = X_{1,n} Z_{1,n}','\Delta S_{2,n} = X_{2,n} Z_{2,n}','\DeltaS_{n}^\times    = \Delta S_{dep,n} \Delta S_{fac,n}');

if KSE == 2
    
    SpkFreqin = 100;
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    tspk = (1:Nspk)*SpkPerin;

    Dtatspkin = SpkPerin;
    Tspkmax = floor(Tmax/SpkPerin);
    tspk = SpkPerin:SpkPerin:Tspkmax*SpkPerin;
     jspk = floor(tspk/dt);
  
    taudep1 = 500;
    taudep2vec = 10:10:1000;
    taufac1 = 100;
    taufac2vec = 10:10:1000;
    
    Pdep1 =  exp(-Dtatspkin/taudep1);
    Pfac1 =  exp(-Dtatspkin/taufac1);

    sgmadep2vec = zeros(1,length(taudep2vec));
    sgmafac2vec = zeros(1,length(taufac2vec));
    sgmadepcrossvec = zeros(1,length(taudep2vec));
    sgmafaccrossvec = zeros(1,length(taufac2vec));
    
    X1 = zeros(1,Tspkmax);
    Z1 = zeros(1,Tspkmax);
    Zaux1 = zeros(1,Tspkmax);
    X2 = zeros(1,Tspkmax);
    Z2 = zeros(1,Tspkmax);
    Zaux2 = zeros(1,Tspkmax);
    
    X1(1) = 1;
    Z1(1) = afac;
    X2(1) = 1;
    Z2(1) = afac;
    
    for j=1:Tspkmax-1
        X1(j+1) = Xinf+((1-adep)*X1(j)-Xinf)*Pdep1;       
        Z1(j+1) = (1-afac)*(Zinf+(Z1(j)-Zinf)*Pfac1)+afac;
    end
    X1fp = (1-Pdep1)*Xinf/(1-(1-adep)*Pdep1);
    Z1fp = ((1-Pfac1)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac1);
    
    Fcut1 = 0.5*(X1(1)-X1fp)+X1fp;
    for k=2:length(tspk)
        if X1(k)<Fcut1
            kaux = k;    
            break
        end
    end
    sgmadep1 = (tspk(kaux)-tspk(1))/log((1-X1fp)/(X1(kaux)-X1fp));
    F1 = X1fp+(1-X1fp)*exp(-(tspk-tspk(1))/sgmadep1);
    E_F1 = sum((X1-F1).^2);

    for l=1:length(taudep2vec)
        taudep2 = taudep2vec(l);      
        Pdep2 =  exp(-Dtatspkin/taudep2);
        for j=1:Tspkmax-1
            X2(j+1) = Xinf+((1-adep)*X2(j)-Xinf)*Pdep2;       
        end
        X2fp = (1-Pdep2)*Xinf/(1-(1-adep)*Pdep2);      
        Mcross_depfp = ((1-etadep)*X1fp+etadep*X2fp);
        Mcross_dep = ((1-etadep)*X1+etadep*X2);       
        
        Fcut2 = 0.5*(X2(1)-X2fp)+X2fp;
        for k=2:length(tspk)
            if X2(k)<Fcut2
                kaux = k;    
                break
            end
        end
        sgmadep2 = (tspk(kaux)-tspk(1))/log((1-X2fp)/(X2(kaux)-X2fp));
        F2 = X2fp+(1-X2fp)*exp(-(tspk-tspk(1))/sgmadep2);
        E_F2 = sum((X2-F2).^2);

        Fcut_cross = 0.5*(Mcross_dep(1)-Mcross_depfp)+Mcross_depfp;
        for k=2:length(tspk)
            if Mcross_dep(k)<Fcut_cross
                kaux=k;
                break;
            end
        end
        sgmadepcross = (tspk(kaux)-tspk(1))/log((1-Mcross_depfp)/(Mcross_dep(kaux)-Mcross_depfp));
        Fcross = Mcross_depfp+(1-Mcross_depfp)*exp(-(tspk-tspk(1))/sgmadepcross);
        E_Fcross = sum((Mcross_dep-Fcross).^2);
         
        sgmadep2vec(l) = sgmadep2;
        sgmadepcrossvec(l) = sgmadepcross;
        
        %[l E_F1 E_F2 E_Fcross]
        
    end
    
    for l=1:length(taudep2vec)
        taufac2 = taufac2vec(l);
        Pfac2 =  exp(-Dtatspkin/taufac2); 
        for j=1:Tspkmax-1      
            Z2(j+1) = (1-afac)*(Zinf+(Z2(j)-Zinf)*Pfac2)+afac;
        end
        Z2fp = ((1-Pfac2)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac2);
        Mcross_facfp = ((1-etafac)*Z1fp+etafac*Z2fp);
        Mcross_fac = ((1-etafac)*Z1+etafac*Z2);
        Gcut2 = 0.5*(Z2fp-Z2(1))+Z2(1);
        for k=2:length(tspk)
            if Z2(k)>Gcut2
                kaux = k;    
                break
            end
        end
        C2 = 1-afac/Z2fp;
        sgmafac2 = (tspk(kaux)-tspk(1))/log(Z2fp*C2/(Z2fp-Z2(kaux)));
        G2 = Z2fp*(1-C2*exp(-(tspk-tspk(1))/sgmafac2));
        E_G2 = sum((Z2-G2).^2);
        Gcut_cross = 0.5*(Mcross_facfp-Mcross_fac(1))+Mcross_fac(1);
        for k=2:length(tspk)
            if Mcross_fac(k)>Gcut_cross
                kaux = k;    
                break
            end
        end
        Ccross = 1-afac/Mcross_facfp;
        sgmafaccross = (tspk(kaux)-tspk(1))/log(Mcross_facfp*Ccross/(Mcross_facfp-Mcross_fac(kaux)));
        Gcross = Mcross_facfp*(1-Ccross*exp(-(tspk-tspk(1))/sgmafaccross));
        E_Gcross = sum((Mcross_fac-Gcross).^2);

        sgmafac2vec(l) = sgmafac2;
        sgmafaccrossvec(l) = sgmafaccross;

        % [E_G1 E_G2 E_Gcross]
end
    
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(-100,-100,'o','Color',lightgray,'linewidth',2);
    plot(taudep2vec,sgmadep1*ones(1,length(sgmadep2vec)),'o','Color',lightgray)
    plot(taudep2vec,sgmadepcrossvec,'ob')
    plot(taudep2vec,sgmadep2vec,'or')  
    axis([0 1000 0 100])
    set(gca,'fontsize',24);
    xlabel('\tau_{dep,2}');
    legend('\sigma_{d}^\times','\sigma_{d,2}','\sigma_{d,1}');
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(-100,-100,'o','Color',lightgray,'linewidth',2);
    plot(taufac2vec,sgmafac1*ones(1,length(sgmafac2vec)),'o','Color',lightgray)
    plot(taufac2vec,sgmafaccrossvec,'ob')
    plot(taufac2vec,sgmafac2vec,'or')   
    axis([0 1000 0 100])
    set(gca,'fontsize',24);
    xlabel('\tau_{dep,2}');
    legend('\sigma_{f}^\times','\sigma_{f,2}','\sigma_{f,1}');
    
    sgmadepplusfac2vec = (1./sgmadep2vec+1./sgmafac2vec).^-1;
    sgmadepplusfac1vec = (1./sgmadep1+1./sgmafac1).^-1;
    sgmadepplusfaccrossvec = (1./sgmadepcrossvec+1./sgmafaccrossvec).^-1;
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2); 
    plot(-100,-100,'o','Color',lightgray,'linewidth',2);
    plot(taufac2vec,sgmadepplusfac1vec*ones(1,length(sgmafac2vec)),'o','Color',lightgray)
    plot(taufac2vec,sgmadepplusfac2vec,'or')
    plot(taufac2vec,sgmadepplusfaccrossvec,'ob')   
    axis([0 1000 0 100])
    set(gca,'fontsize',24);
    xlabel('\tau_{dep,2}');
    legend('\sigma_{d+f}^\times','\sigma_{d+f,2}','\sigma_{d+f,1}');


end