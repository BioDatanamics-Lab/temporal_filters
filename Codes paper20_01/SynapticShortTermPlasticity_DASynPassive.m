% 2021-04-21

% DA model + synaptic dynamics + passive cell

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
mediumacquamarine = [0.4 0.8 0.6];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];

KSE = 2;
        % 1: single cell
        % 2: Changes in Deltaspk
        
        
% Functions

Q=@(a,tau,Delta) (1-a).*exp(-Delta/tau);
Xbar=@(a,tau,Delta,Xinf) (1-exp(-Delta/tau))*Xinf./(1-(1-a).*exp(-Delta/tau));
Zbar=@(a,tau,Delta,Zinf) ((1-exp(-Delta/tau))*(1-a)*Zinf+a)./(1-(1-a).*exp(-Delta/tau));

% Parameters

Xinf = 1;
Zinf = 0;
adep = 0.1;
afac = 0.1;

taudec = 5;

C = 1;
El = -60;
gL = 0.1;
Iapp = 0;

taudep = 250;
taufac = 250;


Gsyn = 0.1;
Esyn = 0;

% Time definitions

Tmax = 1000;
dt = 0.01;
t = 0:dt:Tmax;

% Input spike definitions

SpkFreqin = 100;
SpkPerin = 1000/SpkFreqin;
Nspk = floor(Tmax/SpkPerin);
tspk = (1:Nspk)*SpkPerin;

if KSE == 1
   
    % Single parameter set
    
     % Voltage response - synaptic dynamics: no short-term plasticity
    
    V = zeros(1,length(t));
    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    x(1) = 1;
    z(1) = 0;
    V(1) = El+Iapp/gL;
    spkcnt = 1;
    for j=1:length(t)-1
        k1v = (-gL*(V(j)-El)+Iapp-Gsyn*S(j)*(V(j)-Esyn))/C;
        k1s = -S(j)/taudec;   
        av = V(j)+k1v*dt;
        as = S(j)+k1s*dt;
        k2v = (-gL*(av-El)+Iapp-Gsyn*as*(av-Esyn))/C;
        k2s = -as/taudec;
        V(j+1) = V(j)+(k1v+k2v)*dt/2;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        if t(j)> tspk(spkcnt)
            DeltaS(spkcnt) = afac;
            S(j+1) = DeltaS(spkcnt) ;
            spkcnt = spkcnt+1;         
        end
    end
    
    Vb = V;
    Sb = S;
    
    % Sequence evolution (numerical)
    
    [Vpeakb,tpeakb,cntpeakb] = PeakOsc(V,t,tspk);
    [Vtroughb,ttroughb,cnttroughb] = TroughOsc(V,t,tspk);
    [Speakb,tpeaksb,cntpeaksb] = PeakOsc(S,t,tspk);
   
    
    % Voltage response - synaptic dynamics: short-term plasticity

    V = zeros(1,length(t));
    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    x(1) = 1;
    z(1) = 0;
    V(1) = El+Iapp/gL;
    spkcnt = 1;
    for j=1:length(t)-1
        k1v = (-gL*(V(j)-El)+Iapp-Gsyn*S(j)*(V(j)-Esyn))/C;
        k1s = -S(j)/taudec;   
        k1x = (Xinf-x(j))/taudep;
        k1z = (Zinf-z(j))/taufac;
        av = V(j)+k1v*dt;
        as = S(j)+k1s*dt;
        ax = x(j)+k1x*dt;
        az = z(j)+k1z*dt;
        k2v = (-gL*(av-El)+Iapp-Gsyn*as*(av-Esyn))/C;
        k2s = -as/taudec;
        k2x = (Xinf-ax)/taudep;
        k2z = (Zinf-az)/taufac;
        V(j+1) = V(j)+(k1v+k2v)*dt/2;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        x(j+1) = x(j)+(k1x+k2x)*dt/2;
        z(j+1) = z(j)+(k1z+k2z)*dt/2;
        if t(j)> tspk(spkcnt)
            x(j+1) = x(j+1)-adep*x(j);
            z(j+1) = z(j+1)+afac*(1-z(j+1));
            DeltaS(spkcnt) = x(j)*z(j+1);
            S(j+1) = DeltaS(spkcnt) ;
            spkcnt = spkcnt+1;         
        end
    end
   
    
    [Vpeak,tpeak,cntpeak] = PeakOsc(V,t,tspk);
    [Vtrough,ttrough,cnttrough] = TroughOsc(V,t,tspk);
    [Speak,tpeaks,cntpeaks] = PeakOsc(S,t,tspk);
    
    figure
    hold on
    plot(t,x,'-r','linewidth',2);
    plot(t,z,'-g','linewidth',2);
    plot(t,x.*z,'-b','linewidth',2);
    plot(t,S,'-','Color',lightblueish,'linewidth',2);
    axis([0 Tmax 0 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('x','z','\DeltaS = x z');
    
%     figure
%     hold on
%     plot(t,Vb,'-b','linewidth',2);
%     plot(t,Sb*10-65,'-','Color',lightblueish,'linewidth',2);
%     set(gca,'fontsize',24);
%     axis([0 Tmax -70 -50]);
%     xlabel('t  [ms]');
%     ylabel('V  [mV]');
    
    figure
    hold on
    plot(t,Vb,'-r','linewidth',2);
    plot(t,Sb*10-65,'-','Color',lightcoral,'linewidth',2);
    plot(t,V,'-b','linewidth',2);
    plot(t,S*10-65,'-','Color',lightblueish,'linewidth',2);
    set(gca,'fontsize',24);
    axis([0 Tmax -70 -50]);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    
    
    figure
    hold on
    plot(tpeak,Vpeak-El,'ob','linewidth',2);
    plot(tpeakb,Vpeakb-El,'or','linewidth',2);
    plot(tpeaks,Speak,'og','linewidth',2);
    plot(tpeak,Vpeak-Vtrough,'-','Color',lightblueish,'linewidth',2);    
    plot(tpeakb,Vpeakb-Vtroughb,'-','Color',lightcoral,'linewidth',2);
    plot(tpeaksb,Speakb,'--g','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 Tmax 0 3]);
    xlabel('t  [ms]');
    ylabel('');
    legend('V_n','V_0','S_n');
    
    Vpeakmax = max(Vpeak(1:cntpeak));
    Vpeakbmax = max(Vpeakb(1:cntpeakb));
    Speakmax = max(Speak(1:cntpeaks));
     
    
    figure
    hold on
    plot(tpeak,(Vpeak-El)/(Vpeakmax-El),'ob','linewidth',2);
    plot(tpeakb,(Vpeakb-El)/(Vpeakbmax-El),'or','linewidth',2);
    plot(tpeaks,Speak/Speakmax,'og','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 Tmax 0 1.3]);
    xlabel('t  [ms]');
    ylabel('');
    legend('V_n','V_0','S_n');  
    legend(['V_n       V_{n,max}    = ' num2str(Vpeakmax,'%2.1f')],['V^0_{n}       V^0_{n,max}    = ' num2str(Vpeakbmax,'%2.1f')],['S_n       S_{n,max}    = ' num2str(Speakmax,'%2.1f')],'Location','SouthEast');
    text(40,1.15,['\tau_m = '  num2str(C/gL), '          \tau_{dep} = '  num2str(taudep), '   \tau_{fac} = '  num2str(taufac)],'fontsize',24) ;
    
    
    figure
    hold on
    plot(t,V,'-b','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 Tmax -62 -50]);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
%     figure
%     hold on
%     plot(tpeak,Vpeak-El,'ob','linewidth',2);
%     plot(tpeakb,Vpeakb-El,'or','linewidth',2);
%     plot(tpeaks,Speak,'og','linewidth',2);
%     plot(tpeak,Vpeak-Vtrough,'-','Color',lightblueish,'linewidth',2);    
%     plot(tpeakb,Vpeakb-Vtroughb,'-','Color',lightcoral,'linewidth',2);
%     plot(tpeaksb,Speakb,'--g','linewidth',2);
%     set(gca,'fontsize',24);
%     axis([0 Tmax 0 3]);
%     xlabel('t  [ms]');
%     ylabel('');
%     legend('V_n','V_0','S_n');
    
elseif KSE == 2
    
    SpkFreqin = 40;
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    tspk = (1:Nspk)*SpkPerin;
    
    % Voltage response - synaptic dynamics: short-term plasticity

    V = zeros(1,length(t));
    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    x(1) = 1;
    z(1) = 0;
    V(1) = El+Iapp/gL;
    spkcnt = 1;
    for j=1:length(t)-1
        k1v = (-gL*(V(j)-El)+Iapp-Gsyn*S(j)*(V(j)-Esyn))/C;
        k1s = -S(j)/taudec;   
        k1x = (Xinf-x(j))/taudep;
        k1z = (Zinf-z(j))/taufac;
        av = V(j)+k1v*dt;
        as = S(j)+k1s*dt;
        ax = x(j)+k1x*dt;
        az = z(j)+k1z*dt;
        k2v = (-gL*(av-El)+Iapp-Gsyn*as*(av-Esyn))/C;
        k2s = -as/taudec;
        k2x = (Xinf-ax)/taudep;
        k2z = (Zinf-az)/taufac;
        V(j+1) = V(j)+(k1v+k2v)*dt/2;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        x(j+1) = x(j)+(k1x+k2x)*dt/2;
        z(j+1) = z(j)+(k1z+k2z)*dt/2;
        if t(j)> tspk(spkcnt)
            x(j+1) = x(j+1)-adep*x(j);
            z(j+1) = z(j+1)+afac*(1-z(j+1));
            DeltaS(spkcnt) = x(j)*z(j+1);
            S(j+1) = DeltaS(spkcnt) ;
            spkcnt = spkcnt+1;         
        end
    end
   
    
    [Vpeak1,tpeak1,cntpeak1] = PeakOsc(V,t,tspk);
   
    
    figure
    hold on
    plot(t,x,'-r','linewidth',2);
    plot(t,z,'-g','linewidth',2);
    plot(t,x.*z,'-b','linewidth',2);
    plot(t,S,'-','Color',lightblueish,'linewidth',2);
    axis([0 Tmax 0 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('x','z','\DeltaS = x z');
    
    figure
    hold on
    plot(t,V,'-b','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 Tmax -62 -50]);
    xlabel('t  [ms]');
    ylabel('V  [mV]');
    
    SpkFreqin = 80;
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    tspk = (1:Nspk)*SpkPerin;
    
    % Voltage response - synaptic dynamics: short-term plasticity

    V = zeros(1,length(t));
    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    x(1) = 1;
    z(1) = 0;
    V(1) = El+Iapp/gL;
    spkcnt = 1;
    for j=1:length(t)-1
        k1v = (-gL*(V(j)-El)+Iapp-Gsyn*S(j)*(V(j)-Esyn))/C;
        k1s = -S(j)/taudec;   
        k1x = (Xinf-x(j))/taudep;
        k1z = (Zinf-z(j))/taufac;
        av = V(j)+k1v*dt;
        as = S(j)+k1s*dt;
        ax = x(j)+k1x*dt;
        az = z(j)+k1z*dt;
        k2v = (-gL*(av-El)+Iapp-Gsyn*as*(av-Esyn))/C;
        k2s = -as/taudec;
        k2x = (Xinf-ax)/taudep;
        k2z = (Zinf-az)/taufac;
        V(j+1) = V(j)+(k1v+k2v)*dt/2;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        x(j+1) = x(j)+(k1x+k2x)*dt/2;
        z(j+1) = z(j)+(k1z+k2z)*dt/2;
        if t(j)> tspk(spkcnt)
            x(j+1) = x(j+1)-adep*x(j);
            z(j+1) = z(j+1)+afac*(1-z(j+1));
            DeltaS(spkcnt) = x(j)*z(j+1);
            S(j+1) = DeltaS(spkcnt) ;
            spkcnt = spkcnt+1;         
        end
    end
   
    [Vpeak2,tpeak2,cntpeak2] = PeakOsc(V,t,tspk);
    
    SpkFreqin = 121;
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    tspk = (1:Nspk)*SpkPerin;
    
    % Voltage response - synaptic dynamics: short-term plasticity

    V = zeros(1,length(t));
    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    x(1) = 1;
    z(1) = 0;
    V(1) = El+Iapp/gL;
    spkcnt = 1;
    for j=1:length(t)-1
        k1v = (-gL*(V(j)-El)+Iapp-Gsyn*S(j)*(V(j)-Esyn))/C;
        k1s = -S(j)/taudec;   
        k1x = (Xinf-x(j))/taudep;
        k1z = (Zinf-z(j))/taufac;
        av = V(j)+k1v*dt;
        as = S(j)+k1s*dt;
        ax = x(j)+k1x*dt;
        az = z(j)+k1z*dt;
        k2v = (-gL*(av-El)+Iapp-Gsyn*as*(av-Esyn))/C;
        k2s = -as/taudec;
        k2x = (Xinf-ax)/taudep;
        k2z = (Zinf-az)/taufac;
        V(j+1) = V(j)+(k1v+k2v)*dt/2;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        x(j+1) = x(j)+(k1x+k2x)*dt/2;
        z(j+1) = z(j)+(k1z+k2z)*dt/2;
        if t(j)> tspk(spkcnt)
            x(j+1) = x(j+1)-adep*x(j);
            z(j+1) = z(j+1)+afac*(1-z(j+1));
            DeltaS(spkcnt) = x(j)*z(j+1);
            S(j+1) = DeltaS(spkcnt) ;
            spkcnt = spkcnt+1;         
        end
    end
   
    
    [Vpeak3,tpeak3,cntpeak3] = PeakOsc(V,t,tspk);
    
    figure
    hold on
    plot(tpeak1,Vpeak1,'ob','linewidth',2);
    plot(tpeak2,Vpeak2,'or','linewidth',2);
    plot(tpeak3,Vpeak3,'og','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 Tmax -60 -50]);
    xlabel('t  [ms]');
    ylabel('V_n  [mV]');
    legend('f_{spk}=40','f_{spk}=80','f_{spk}=120');
    
    
end





function [Speak,tpeak,cntpeak] = PeakOsc(S,t,tspk) 

    dt = t(2)-t(1);
    jspk = floor(tspk/dt);
    
    
    Speak = zeros(1,length(tspk)-1);
    tpeak = zeros(1,length(tspk)-1);
    cntpeak = 0;
    for j=jspk(1):length(t)-1
        if S(j)>S(j-1) && S(j)>S(j+1)
            cntpeak = cntpeak+1;
            Speak(cntpeak) = S(j);
            tpeak(cntpeak) = t(j);
        end
    end
end

function [Strough,ttrough,cnttrough] = TroughOsc(S,t,tspk) 

    dt = t(2)-t(1);
    jspk = floor(tspk/dt);  
    

    Strough = zeros(1,length(tspk)-1);
    ttrough = zeros(1,length(tspk)-1);
    cnttrough = 0;
    for j=jspk(1):length(t)-1
        if S(j)<S(j-1) && S(j)<S(j+1)
            cnttrough = cnttrough+1;
            Strough(cnttrough) = S(j);
            ttrough(cnttrough) = t(j);
        end    
    end
    cnttrough = cnttrough+1;
    Strough(cnttrough) = S(end);
    ttrough(cnttrough) = t(end);
end


















