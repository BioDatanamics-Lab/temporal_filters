% 2021-04-14

% DA model + synaptic dynamics with time scales mimicking the membrane
% potential dynamics

clearvars;
close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
mediumacquamarine = [0.4 0.8 0.6];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];

KSE = 1;
        % 1: Synaptic update: additive S --> S+Delta S with different time
        % 2: Sn: Changes in taudec (imagesc)
        % 3: Sn: Changes in taudep (imagesc
        % 4: Sn: Changes in taudep (imagesc)
        % 5: Sn: Changes in taufac (imagesc)
        % 1-5: single depression-facilitation process
        % 1: 

% Functions

Q=@(a,tau,Delta) (1-a).*exp(-Delta/tau);
Xbar=@(a,tau,Delta,Xinf) (1-exp(-Delta/tau))*Xinf./(1-(1-a).*exp(-Delta/tau));
Zbar=@(a,tau,Delta,Zinf) ((1-exp(-Delta/tau))*(1-a)*Zinf+a)./(1-(1-a).*exp(-Delta/tau));

% Parameters

Xinf = 1;
Zinf = 0;
adep = 0.1;
afac = 0.2;

taudep = 250;
taufac = 250;
taudec = 100;

% Time definitions

Tmax = 1000;
dt = 0.01;
t = 0:dt:Tmax;

% Input spike definitions

SpkFreqin = 100;
SpkPerin = 1000/SpkFreqin;
Nspk = floor(Tmax/SpkPerin);
tspk = (1:Nspk)*SpkPerin;


% Sequences: depression and facilitation

Deltaspk = SpkPerin;

Nspk = 100;
nvec = 1:1:Nspk;
tvec = nvec*Deltaspk;

Qd = Q(adep,taudep,Deltaspk);
Qf = Q(afac,taufac,Deltaspk);
Qdn = Qd.^(nvec-1);
Qfn = Qf.^(nvec-1);
Xb = Xbar(adep,taudep,Deltaspk,Xinf);
Zb = Zbar(afac,taufac,Deltaspk,Zinf);

Xn = Xb*(1-Qdn)+Qdn*1; 
Zn = Zb*(1-Qfn)+Qfn*afac; 
DeltaSn = Xn.*Zn;



if KSE == 1
   
    % Single parameter set
    
    
    % Depression & facilitation long-term time scales (analytical)
    
    sgma_dep = Deltaspk./(Deltaspk./taudep-log(1-adep));
    sgma_fac = Deltaspk./(Deltaspk./taufac-log(1-afac));
    sgma_depfac = Deltaspk./(Deltaspk*(1/taudep+1/taufac)-log(1-adep)-log(1-afac));
    [sgma_dep sgma_fac sgma_depfac]

    % Synaptic dynamics: no short-term plasticity - baseline

    S = zeros(1,length(t));

    spkcnt = 1;
    for j=1:length(t)-1
        k1s = -S(j)/taudec;
        as = S(j)+k1s*dt;
        k2s = -as/taudec;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        if t(j)> tspk(spkcnt)
            S(j+1) = S(j+1)+afac;
            spkcnt = spkcnt+1;
        end
    end
    Sbase = S;
    
     % Sequence evolution (numerical)
     
    [Speakb,tpeakb,cntpeakb] = PeakOsc(Sbase,t,tspk);
    [Stroughb,ttroughb,cnttroughb] = TroughOsc(Sbase,t,tspk);
    
%     figure
%     hold on
%     plot(t,S,'-b','linewidth',2);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     ylabel('S');
%     axis([0 Tmax 0 1.2]);


    % Synaptic dynamics: short-term plasticity

    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    x(1) = 1;
    z(1) = 0;
    spkcnt = 1;
    for j=1:length(t)-1
        k1s = -S(j)/taudec;   
        k1x = (Xinf-x(j))/taudep;
        k1z = (Zinf-z(j))/taufac;
        as = S(j)+k1s*dt;
        ax = x(j)+k1x*dt;
        az = z(j)+k1z*dt;
        k2s = -as/taudec;
        k2x = (Xinf-ax)/taudep;
        k2z = (Zinf-az)/taufac;
        S(j+1) = S(j)+(k1s+k2s)*dt/2;
        x(j+1) = x(j)+(k1x+k2x)*dt/2;
        z(j+1) = z(j)+(k1z+k2z)*dt/2;
        if t(j)> tspk(spkcnt)
            x(j+1) = x(j+1)-adep*x(j);
            z(j+1) = z(j+1)+afac*(1-z(j+1));
            DeltaS(spkcnt) = x(j)*z(j+1);
            S(j+1) = S(j+1)+DeltaS(spkcnt);
            spkcnt = spkcnt+1;         
        end
    end
   
    % Sequence evolution (numerical)
    
    [Speak,tpeak,cntpeak] = PeakOsc(S,t,tspk);
    [Strough,ttrough,cnttrough] = TroughOsc(S,t,tspk);
  
    % Sequence evolution (analytical)
    
    Sn = zeros(1,length(nvec));
    Sn(1) = afac;
    A = exp(-Deltaspk/taudec);
    for k=2:length(nvec)
        Sn(k) = A^(k-1)*Sn(1);
        Sumaux = DeltaSn(k);
        for l=1:k-2
            Sumaux = Sumaux+A^l*DeltaSn(k-l);
        end
        Sn(k) = Sn(k)+Sumaux;
    end
    
    Snbbar = afac/(1-exp(-Deltaspk/taudec));
    Snb = Snbbar+exp(-(nvec-1)*Deltaspk/taudec)*(afac-Snbbar);
    
    
    % Sequence steady-state (analytical)
    
    Sumaux = 1;
    for l=1:10000
        Sumaux = Sumaux+A^l;
    end
    Snbar0 = Sumaux*DeltaSn(end);
    
    Snbar = DeltaSn(end)/(1-exp(-Deltaspk/taudec));
    
    for k=1:length(tspk)
        if Speak(k)>0.63*Snbar
            sgmadec_n = k*Deltaspk;
            break;
        end
    end
    for k=1:length(tspk)
        if Speakb(k)>0.63*Snbbar
            taudec_n = k*Deltaspk;
            break;
        end
    end
    

    figure
    hold on
    plot(t,x,'-r','linewidth',2);
    plot(t,z,'-g','linewidth',2);
    plot(t,x.*z,'-b','linewidth',2);
    axis([0 Tmax 0 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('x','z','\DeltaS = x z');

%     figure
%     hold on
%     plot(t,S,'-b','linewidth',2);
%     plot(t,Sbase,'-r','linewidth',2);
%     set(gca,'fontsize',24);
%     axis([0 Tmax 0 1.2]);
%     xlabel('t  [ms]');
%     xlabel('S  [ms]');
%     legend
%     
%     figure
%     hold on
%     plot(tpeak,Speak-Strough,'ob','linewidth',2);
%     plot(tpeakb,Speakb-Stroughb,'or','linewidth',2);
%     set(gca,'fontsize',24);
%     axis([0 Tmax 0 1.2]);
%     xlabel('t  [ms]');
%     xlabel('S  [ms]');

    Speakmax = max(Speak);
    Speakbmax = max(Speakb);
    DeltaSmax = max(DeltaS);
    
    figure
    hold on
    plot(t,S,'-b','linewidth',2);
    plot(t,Sbase,'-r','linewidth',2);
    plot(tpeak,Speak-Strough,'-','Color',lightblueish,'linewidth',2);
    plot(tpeakb,Speakb-Stroughb,'-','Color',lightcoral,'linewidth',2);
    plot(tvec,DeltaSn,'og','linewidth',2);
    %plot(tvec,DeltaSn,'o','Color',lightgray,'linewidth',2);
    set(gca,'fontsize',24);
    axis([0 Tmax 0 1.2]);
    xlabel('t  [ms]');
    ylabel('');
    legend('S','S_{base}','A','A_{base}','\Delta S_n');
    
    
    % Curves not normalized
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(-100,-100,'og','linewidth',2);
    plot(Speakb,'or','linewidth',2) 
    plot(DeltaS,'og','linewidth',2);  
    plot(Speak,'ob','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 spkcnt 0 2]);
    xlabel('Spike # (n)');
    legend('S_n','S_{b,n} ','\DeltaS_n','Location','SouthEast');
    
    
    % Curves normlized by their steady state values
    
%     figure
%     hold on
%     plot(-100,-100,'ob','linewidth',2);
%     plot(-100,-100,'or','linewidth',2);
%     plot(-100,-100,'og','linewidth',2);
%     plot(Speakb/Speakb(end),'or','linewidth',2) 
%     plot(DeltaS/DeltaS(end),'og','linewidth',2);  
%     plot(Speak/Speak(end),'ob','linewidth',2);
%     set(gca,'fontsize',24);
%     axis([0 spkcnt 0 2]);
%     xlabel('Spike # (n)');
%     legend('S_n','S_{b,n} ','\DeltaS_n','Location','SouthEast');
%     %text(10,0.1,['\tau_{dep} = '  num2str(taudep), '   \tau_{fac} = '  num2str(taufac) '   \tau_{dec} = '  num2str(taudec)],'fontsize',24) ;  
%     %text(25,0.1,['\tau_{dep} = '  num2str(taudep), '   \tau_{fac} = '  num2str(taufac)],'fontsize',24) ;  
%     text(2,1.1,['\tau_{dec} = '  num2str(taudec), '          \tau_{dep} = '  num2str(taudep), '   \tau_{fac} = '  num2str(taufac)],'fontsize',24) ;
%     text(50,0.6,['$\bar{S}$  = ' num2str(Speak(end),'%2.1f')],'Interpreter','Latex','fontsize',24)
%     text(50,0.4,['$\bar{S}_b$ = ' num2str(Speakb(end),'%2.1f')],'Interpreter','Latex','fontsize',24)
%     text(50,0.2,['$\bar{\Delta S}$ = ' num2str(DeltaS(end),'%2.1f')],'Interpreter','Latex','fontsize',24)
%    
    
    % Curves normalized by their maxima
    
    figure
    hold on
    plot(-100,-100,'ob','linewidth',2);
    plot(-100,-100,'or','linewidth',2);
    plot(-100,-100,'og','linewidth',2);
    plot(Speakb/Speakbmax,'or','linewidth',2) 
    plot(DeltaS/DeltaSmax,'og','linewidth',2);  
    plot(Speak/Speakmax,'ob','linewidth',2);
    set(gca,'fontsize',24);
    axis([0 spkcnt 0 1.3]);
    xlabel('Spike # (n)');
    %legend('S_n','S_{b,n} ','\DeltaS_n','Location','SouthEast');
    legend(['S_n       S_{n,max}    = ' num2str(Speakmax,'%2.1f')],['S^0_{n}       S^0_{n,max}    = ' num2str(Speakbmax,'%2.1f')],['\DeltaS_n    \DeltaS_{n,max} = ' num2str(DeltaSmax,'%2.1f')],'Location','SouthEast');
    axis([0 spkcnt 0 1.3]);
    text(5,1.15,['\tau_{dec} = '  num2str(taudec), '          \tau_{dep} = '  num2str(taudep), '   \tau_{fac} = '  num2str(taufac)],'fontsize',24) ;
    
    

    

    
    
%     figure
%     hold on
%     plot(Speak,'ob','linewidth',2)
%     %plot(Speakb,'or','linewidth',2)
%     plot(DeltaS,'og','linewidth',2);
%     %plot(Speak./(DeltaS/afac),'o','Color',lightblueish,'linewidth',2)
%     plot(Speakb.*(DeltaS/afac),'or','linewidth',2)
%     set(gca,'fontsize',24);
%     axis([0 spkcnt 0 1]);
%     xlabel('Spike #');
%     legend('S','S_{base}','\Delta S_n');
%     

elseif KSE == 2
    
    taudep = 500;
    taufac = 0.1;
    
    % Synaptic dynamics: short-term plasticity

    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    taudec_vec = 0:10:350;
    taudec_vec(1) = 0.1;
    Speak_vec = zeros(length(taudec_vec),1);
    
    for k=1:length(taudec_vec)
    
        k 
        taudec = taudec_vec(k);
        
        x(1) = 1;
        z(1) = 0;
        spkcnt = 1;
        for j=1:length(t)-1
            k1s = -S(j)/taudec;   
            k1x = (Xinf-x(j))/taudep;
            k1z = (Zinf-z(j))/taufac;
            as = S(j)+k1s*dt;
            ax = x(j)+k1x*dt;
            az = z(j)+k1z*dt;
            k2s = -as/taudec;
            k2x = (Xinf-ax)/taudep;
            k2z = (Zinf-az)/taufac;
            S(j+1) = S(j)+(k1s+k2s)*dt/2;
            x(j+1) = x(j)+(k1x+k2x)*dt/2;
            z(j+1) = z(j)+(k1z+k2z)*dt/2;
            if t(j)> tspk(spkcnt)
                x(j+1) = x(j+1)-adep*x(j);
                z(j+1) = z(j+1)+afac*(1-z(j+1));
                DeltaS(spkcnt) = x(j)*z(j+1);
                S(j+1) = S(j+1)+DeltaS(spkcnt);
                spkcnt = spkcnt+1;         
            end
        end

        % Sequence evolution (numerical)

        [Speak,tpeak,cntpeak] = PeakOsc(S,t,tspk);
        [Strough,ttrough,cnttrough] = TroughOsc(S,t,tspk);
        
        Speak_vec(k,1:spkcnt-1) = Speak;
        
        
    end
    
    figure
    imagesc(Speak_vec)
    colorbar
    set(gca,'Ydir','Normal');
    set(gca,'fontsize',24);
    xlabel('Spike # (n)');
    ylabel('\tau_{dec}');
    h = colorbar;
    %set(get(h,'label'),'string','S_n');
    set(get(h,'title'),'string','S_n');
    title([' \tau_{dep} = '  num2str(taudep), '      \tau_{fac} = '  num2str(taufac)],'fontsize',24) ;
    
    
elseif KSE == 3   

    taufac = 0.1;
    taudec = 100;
    
    % Synaptic dynamics: short-term plasticity

    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    taudep_vec = 10:10:300;
    if taudep_vec(1) == 0
        taudep_vec(1) = 0.1;
    end
    Speak_vec = zeros(length(taudep_vec),1);
    
    for k=1:length(taudep_vec)
    
        k 
        taudep = taudep_vec(k);
        
        x(1) = 1;
        z(1) = 0;
        spkcnt = 1;
        for j=1:length(t)-1
            k1s = -S(j)/taudec;   
            k1x = (Xinf-x(j))/taudep;
            k1z = (Zinf-z(j))/taufac;
            as = S(j)+k1s*dt;
            ax = x(j)+k1x*dt;
            az = z(j)+k1z*dt;
            k2s = -as/taudec;
            k2x = (Xinf-ax)/taudep;
            k2z = (Zinf-az)/taufac;
            S(j+1) = S(j)+(k1s+k2s)*dt/2;
            x(j+1) = x(j)+(k1x+k2x)*dt/2;
            z(j+1) = z(j)+(k1z+k2z)*dt/2;
            if t(j)> tspk(spkcnt)
                x(j+1) = x(j+1)-adep*x(j);
                z(j+1) = z(j+1)+afac*(1-z(j+1));
                DeltaS(spkcnt) = x(j)*z(j+1);
                S(j+1) = S(j+1)+DeltaS(spkcnt);
                spkcnt = spkcnt+1;         
            end
        end

        % Sequence evolution (numerical)

        [Speak,tpeak,cntpeak] = PeakOsc(S,t,tspk);
        [Strough,ttrough,cnttrough] = TroughOsc(S,t,tspk);
        
        Speak_vec(k,1:spkcnt-1) = Speak;
        
        
    end
    
    figure
    imagesc(Speak_vec)
    colorbar
    set(gca,'Ydir','Normal');
    set(gca,'fontsize',24);
    xlabel('Spike # (n)');
    ylabel('\tau_{dep}');
    h = colorbar;
    %set(get(h,'label'),'string','S_n');
    set(get(h,'title'),'string','S_n');
    title([' \tau_{dec} = '  num2str(taudec), '      \tau_{fac} = '  num2str(taufac)],'fontsize',24) ;
    
    
    
elseif KSE == 4   

    taufac = 250;
    taudec = 10;
    
    % Synaptic dynamics: short-term plasticity

    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    taudep_vec = 150:10:500;
    if taudep_vec(1) == 0
        taudep_vec(1) = 0.1;
    end
    Speak_vec = zeros(length(taudep_vec),1);
    
    for k=1:length(taudep_vec)
    
        k 
        taudep = taudep_vec(k);
        
        x(1) = 1;
        z(1) = 0;
        spkcnt = 1;
        for j=1:length(t)-1
            k1s = -S(j)/taudec;   
            k1x = (Xinf-x(j))/taudep;
            k1z = (Zinf-z(j))/taufac;
            as = S(j)+k1s*dt;
            ax = x(j)+k1x*dt;
            az = z(j)+k1z*dt;
            k2s = -as/taudec;
            k2x = (Xinf-ax)/taudep;
            k2z = (Zinf-az)/taufac;
            S(j+1) = S(j)+(k1s+k2s)*dt/2;
            x(j+1) = x(j)+(k1x+k2x)*dt/2;
            z(j+1) = z(j)+(k1z+k2z)*dt/2;
            if t(j)> tspk(spkcnt)
                x(j+1) = x(j+1)-adep*x(j);
                z(j+1) = z(j+1)+afac*(1-z(j+1));
                DeltaS(spkcnt) = x(j)*z(j+1);
                S(j+1) = S(j+1)+DeltaS(spkcnt);
                spkcnt = spkcnt+1;         
            end
        end

        % Sequence evolution (numerical)

        [Speak,tpeak,cntpeak] = PeakOsc(S,t,tspk);
        [Strough,ttrough,cnttrough] = TroughOsc(S,t,tspk);
        
        Speak_vec(k,1:spkcnt-1) = Speak;
        
        
    end
    
    figure
    imagesc(Speak_vec)
    colorbar
    set(gca,'Ydir','Normal');
    set(gca,'fontsize',24);
    xlabel('Spike # (n)');
    ylabel('\tau_{dep}');
    h = colorbar;
    %set(get(h,'label'),'string','S_n');
    set(get(h,'title'),'string','S_n');
    title([' \tau_{dec} = '  num2str(taudec), '      \tau_{fac} = '  num2str(taufac)],'fontsize',24) ;
    
    
elseif KSE == 5   

    taudep = 250;
    taudec = 50;
    
    % Synaptic dynamics: short-term plasticity

    S = zeros(1,length(t));
    x = zeros(1,length(t));
    z = zeros(1,length(t));

    DeltaS = zeros(1);
    
    taufac_vec = 0:10:300;
    if taufac_vec(1) == 0
        taufac_vec(1) = 0.1;
    end
    Speak_vec = zeros(length(taufac_vec),1);
    
    for k=1:length(taufac_vec)
    
        k 
        taufac = taufac_vec(k);
        
        x(1) = 1;
        z(1) = 0;
        spkcnt = 1;
        for j=1:length(t)-1
            k1s = -S(j)/taudec;   
            k1x = (Xinf-x(j))/taudep;
            k1z = (Zinf-z(j))/taufac;
            as = S(j)+k1s*dt;
            ax = x(j)+k1x*dt;
            az = z(j)+k1z*dt;
            k2s = -as/taudec;
            k2x = (Xinf-ax)/taudep;
            k2z = (Zinf-az)/taufac;
            S(j+1) = S(j)+(k1s+k2s)*dt/2;
            x(j+1) = x(j)+(k1x+k2x)*dt/2;
            z(j+1) = z(j)+(k1z+k2z)*dt/2;
            if t(j)> tspk(spkcnt)
                x(j+1) = x(j+1)-adep*x(j);
                z(j+1) = z(j+1)+afac*(1-z(j+1));
                DeltaS(spkcnt) = x(j)*z(j+1);
                S(j+1) = S(j+1)+DeltaS(spkcnt);
                spkcnt = spkcnt+1;         
            end
        end

        % Sequence evolution (numerical)

        [Speak,tpeak,cntpeak] = PeakOsc(S,t,tspk);
        [Strough,ttrough,cnttrough] = TroughOsc(S,t,tspk);
        
        Speak_vec(k,1:spkcnt-1) = Speak;
        
        
    end
    
    figure
    imagesc(Speak_vec)
    colorbar
    set(gca,'Ydir','Normal');
    set(gca,'fontsize',24);
    xlabel('Spike # (n)');
    ylabel('\tau_{fac}');
    h = colorbar;
    %set(get(h,'label'),'string','S_n');
    set(get(h,'title'),'string','S_n');
    title([' \tau_{dec} = '  num2str(taudec), '      \tau_{dep} = '  num2str(taudep)],'fontsize',24) ;
    
    
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




