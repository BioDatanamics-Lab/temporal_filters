% 2020-01-08

% Synaptic depression and facilitation models 
% (Dayan & Abbott, see Ermentrout & Terman book, Chapter 7) 

clearvars;
close all;


KSE = 5;
        % 1: Single input frequency
        %    As in SynapticShortTermPlasticity_ToyModels.m
        % 2: Burst: Two input frequencies
        % 3: Poisson inputs - Synaptic update (Delta S)
        % 4: Perturbations of periodic presynaptic input spike trains:
        %    normally distributed (random) - Synaptic update (Delta S)
        % 5: Poisson inputs - Synaptic response mimicking the voltage
        %    response (membrane potential time scales instad of synaptic
        %    time scales)
        % 
        % 301: Poisson inputs
        %      VarX, VarZ & VarM as a function of the input frequency
        % 302: Poisson inputs
        %      VarX, VarZ & VarM as a function of taudep & taufac
  
lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightsalmon = [0.9 0.5 0.4];
lightgray = [.7 .7 .7];
darkgray = [.3 .3 .3];
darkgray2 = [.1 .1 .1];
mediumacquamarine = [0.4 0.8 0.6];

% Functions

Fsqw=@(t,SpkFreqin,DC) (square(2*pi*SpkFreqin*t/1000-pi/2,DC)+1)/2;

% Parameters

taudec = 5;
taudep = 10;
taufac = 0.1;
adep = 0.1;
afac = 0.2;
Xinf = 1;
Zinf = 0;

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

    m = zeros(1,length(t));
    for j=2:length(t)
        m(j) = x(j-1)*z(j);
    end

    % x-, z- and DeltaS=x*z-traces

    figure
    hold on
    plot(t,x,'-r','linewidth',2);
    plot(t,z,'-g','linewidth',2);
    plot(t,x.*z,'-b','linewidth',2);
    axis([0 Tmax 0 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('x','z','x z');

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
    legend('x','z','x z');

    % Sequence of peaks computed from the analytical solutions to the
    % % differential equations using the previous sequence values to generate the
    % current ones (instead of the previous values in the differential
    % equation discretization). 

    Dtatspkin = SpkPerin;
    Tspkmax = floor(Tmax/SpkPerin);
    tspk = SpkPerin:SpkPerin:Tspkmax*SpkPerin;
    jspk = floor(tspk/dt);


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
    Xfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
    Zfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);


    plot(tspk,X,'o','Color',lightblueish,'linewidth',1);
    plot(tspk,Z,'o','Color',lightblueish,'linewidth',1);
    plot(tspk,M,'o','Color',lightblueish,'linewidth',1);
    legend('x','z','x z','M');



elseif KSE == 2
    
    taudep = 500;
    taufac = 500;
    
    % Time definitions

    Tmax = 1000;
    dt = 0.01;
    t = 0:dt:Tmax;

    
    % Input spike definitions
    
    SpkFreqin1 = 4;
    SpkFreqin2 = 200;
    SpkPerin1 = 1000/SpkFreqin1;
    DC = 50;
    SpkPerin2 = 1000/SpkFreqin2;
    Nspkaux = floor(Tmax/SpkPerin2);
    tspkaux = (1:Nspkaux)*SpkPerin2;
    tspkaux = tspkaux.*Fsqw(tspkaux,SpkFreqin1,DC);
    k = 0;
    tspk = zeros(1);
    for l=1:length(tspkaux)
        if tspkaux(l)>0
            k = k+1;
            tspk(k) = tspkaux(l);
        end
    end
    tspk = [tspk Tmax+dt];
    
   
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

    m = zeros(1,length(t));
    for j=2:length(t)
        m(j) = x(j-1)*z(j);
    end

%     % x-, z- and DeltaS=x*z-traces
% 
%     figure
%     hold on
%     plot(t,x,'-r','linewidth',2);
%     plot(t,z,'-g','linewidth',2);
%     plot(t,x.*z,'-b','linewidth',2);
%     axis([0 Tmax 0 1.2]);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     legend('x','z','x z');

    % x-, z- and m-traces (m=DeltaS) + peaks

    figure
    hold on
    plot(t,x,'-r','linewidth',2);
    plot(t,z,'-g','linewidth',2);
    plot(t,x.*z,'-b','linewidth',2);
    plot(t,m,'-','Color',lightblueish,'linewidth',1);
    plot(tspk(1:end-1),zeros(1,length(tspk(1:end-1))),'ok','linewidth',2);
    %plot(t,m,'-b','linewidth',2);
    axis([0 Tmax 0 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('x','z','x z');

    
    % Sequence of peaks computed from the analytical solutions to the
    % differential equations
   
  
    Tspkmax = length(tspk)-1;
    Dtaspkin = diff(tspk(1:end-1));
    X = zeros(1,Tspkmax);
    Z = zeros(1,Tspkmax);
    X(1) = 1;
    Z(1) = afac;
    for j=1:Tspkmax-1
        X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*exp(-Dtaspkin(j)/taudep);
        Z(j+1) = afac+(1-afac)*(Zinf+(Z(j)-Zinf)*exp(-Dtaspkin(j)/taufac));
    end
    M = X.*Z;
    plot(tspk(1:end-1),X,'o','Color',lightblueish,'linewidth',1);
    plot(tspk(1:end-1),Z,'o','Color',lightblueish,'linewidth',1);
    plot(tspk(1:end-1),M,'o','Color',lightblueish,'linewidth',1);
    legend('x','z','x z','M');
    
    if SpkFreqin1 == 4
        jaux = find(tspk>800,1);
        AmpX = max(X(jaux:end))-min(X(jaux:end));
        AmpZ = max(Z(jaux:end))-min(Z(jaux:end));
        AmpM = max(M(jaux:end))-min(M(jaux:end));
    end
    
    figure
    hold on
    plot(tspk(1:end-1),X,'or','linewidth',1);
    plot(tspk(1:end-1),Z,'og','linewidth',1);
    plot(tspk(1:end-1),M,'ob','linewidth',1);
    plot(tspk(1:end-1),-0.05*ones(1,length(tspk(1:end-1))),'ok','linewidth',2);
    axis([0 Tmax+400 -0.1 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    ylabel('\DeltaS');
    legend('X_n','Z_n','\DeltaS_n');
    text(20,1.12,['f_{spk}=' num2str(SpkFreqin2) '  \tau_{dep}=' num2str(taudep) '  \tau_{fac}=' num2str(taufac)],'fontsize',20)
    text(1000,0.35,['A(X_n)~' num2str(AmpX,2)],'Color','r','fontsize',20)
    text(1000,0.2,['A(Z_n)~' num2str(AmpZ,2)],'Color','g','fontsize',20)
    text(1000,0.05,['A(\DeltaS_n)~' num2str(AmpM,2)],'Color','b','fontsize',20)
    set(gca,'XTick',(0:250:1000))
    
    fprintf('AmpX=%f\t AmpZ=%f\t AmpM=%f\n',AmpX,AmpZ,AmpM);
        
elseif KSE == 3
    
    SPKSE = 2;
                % 1: Poisson
                % 2: Poisson (ISI>ISI_min)
                
    taudep = 1000;
    taufac = 1000;
               
    % Time definitions

    Tmax = 200000;
    dt = 0.01;
    t = 0:dt:Tmax;
                
    % Generation of Poisson-like PWC inputs
    
    
   if SPKSE == 1  
        Freq = 100;
        tr = 1;                                   % trials
        r = Freq/(length(t)-1)*Tmax/1000;          % [rate] = Spk/bin
        spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
        spk = spk(1:length(t));
        tspk = find(spk'>0)*dt;                 % Spike times
        ISI = diff(tspk);    
        CVspk = sqrt(var(ISI))/mean(ISI);
   elseif SPKSE == 2
        ISImin = 1;
        Freq = 100;
        tr = 1;                                 % trials
        r = Freq/(length(t)-1)*Tmax/1000;       % [rate] = Spk/bin
        spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
        spk = spk(1:length(t));
        tspkbase = find(spk'>0)*dt;                 % Spike times
        ISIbase = diff(tspkbase);
        ISI = zeros(1);
        i=0;
        for j=1:length(ISIbase)
            if ISIbase(j) > ISImin
                i=i+1;
                ISI(i) = ISIbase(j);
            end
        end    
        tspk = zeros(1);
        tspk(1) = tspkbase(1);
        for i=2:length(ISI)
            tspk(i) = tspk(i-1)+ISI(i);
        end
        CVspk = sqrt(var(ISI))/(mean(ISI)-ISImin);
   end
   
   tspk = [tspk Tmax+dt];
   
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

    m = zeros(1,length(t));
    for j=2:length(t)
        m(j) = x(j-1)*z(j);
    end
                    
    % x-, z- and m-traces (m=DeltaS) + peaks

%     figure
%     hold on
%     plot(t,x,'-r','linewidth',2);
%     plot(t,z,'-g','linewidth',2);
%     plot(t,x.*z,'-b','linewidth',2);
%     plot(t,m,'-','Color',lightblueish,'linewidth',1);
%     %plot(t,m,'-b','linewidth',2);
%     axis([0 Tmax 0 1.2]);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     legend('x','z','x z');
%     

% Sequences: Poisson spike train inputs
    
    Tspkmax = length(tspk)-1;
    Dtaspkin = diff(tspk(1:end-1));
    X = zeros(1,Tspkmax);
    Z = zeros(1,Tspkmax);
    X(1) = 1;
    Z(1) = afac;
    for j=1:Tspkmax-1
        X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*exp(-Dtaspkin(j)/taudep);
        Z(j+1) = afac+(1-afac)*(Zinf+(Z(j)-Zinf)*exp(-Dtaspkin(j)/taufac));
    end
    M = X.*Z;
    
    
%     figure
%     hold on
%     plot(t,x.*z,'-b','linewidth',2);
%     plot(t,m,'-','Color',lightblueish,'linewidth',1);
%     plot(tspk,zeros(1,length(tspk)),'ok','linewidth',1);
%     %plot(t,m,'-b','linewidth',2);
%     axis([0 Tmax 0 1.2]);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     legend('x z');
%     plot(tspk(1:end-1),X,'o','Color',lightblueish,'linewidth',1);
%     plot(tspk(1:end-1),Z,'o','Color',lightblueish,'linewidth',1);
%     plot(tspk(1:end-1),M,'o','Color',lightblueish,'linewidth',1);
%     legend('x','z','x z','M');
    
%     figure
%     hold on
%     plot(tspk(1:end-1),X,'or','linewidth',1);
%     plot(tspk(1:end-1),Z,'og','linewidth',1);
%     plot(tspk(1:end-1),M,'ob','linewidth',1);
%     axis([0 Tmax 0 1.2]);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     legend('X','Z','M');
    
    
    % Sequences: Periodic spike train inputs (constant frequency input
    % equal to the Poisson mean rate)
    
    SpkFreqin = Freq;
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    tspkc = (1:Nspk)*SpkPerin;
    
    Dtatspkin = SpkPerin;
    Tspkmax = floor(Tmax/SpkPerin);
    tspkc = SpkPerin:SpkPerin:Tspkmax*SpkPerin;
    jspkc = floor(tspk/dt);

    Pdep =  exp(-Dtatspkin/taudep);
    Pfac =  exp(-Dtatspkin/taufac);
    Xc = zeros(1,Tspkmax);
    Zc = zeros(1,Tspkmax);
    Zaux = zeros(1,Tspkmax);
    Xc(1) = 1;
    Zc(1) = afac;
    for j=1:Tspkmax-1
        Xc(j+1) = Xinf+((1-adep)*Xc(j)-Xinf)*Pdep;       
        Zc(j+1) = (1-afac)*(Zinf+(Zc(j)-Zinf)*Pfac)+afac;
    end
    Mc = Xc.*Zc;
    Xcfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
    Zcfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);
    Mcfp = Xcfp*Zcfp;
   
    sgmadep = 1/(1/taudep-log(1-adep)/SpkPerin);
    sgmafac = 1/(1/taufac-log(1-afac)/SpkPerin);
    sgmadepfac = 1/(1/taudep+1/taufac-log((1-adep)*(1-afac))/SpkPerin);
     
%     figure
%     hold on
%     plot(tspk(1:end-1),X,'ob','linewidth',1);
%     plot(tspkc,Xc,'or');
%     plot(tspk,-0.1*ones(1,length(tspk)),'ok','linewidth',1);
%     axis([0 Tmax -0.2 1.2]);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     ylabel('X')
%     legend('X  (Poisson)','X  (deterministic)') 
%     
%     
%     figure
%     hold on
%     plot(tspk(1:end-1),Z,'ob','linewidth',1);
%     plot(tspkc,Zc,'or');
%     plot(tspk,-0.1*ones(1,length(tspk)),'ok','linewidth',1);
%     axis([0 Tmax -0.2 1.2]);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     ylabel('Z')
%     legend('Z  (Poisson)','Z  (deterministic)') 
%     
%     figure
%     hold on
%     plot(tspk(1:end-1),M,'ob','linewidth',1);
%     plot(tspkc,Mc,'or');
%     plot(tspk,-0.1*ones(1,length(tspk)),'ok','linewidth',1);
%     axis([0 Tmax -0.2 1.2]);
%     set(gca,'fontsize',24);
%     xlabel('t  [ms]');
%     ylabel('\DeltaS')
%     legend('\DeltaS  (Poisson)','\DeltaS  (deterministic)')     
    
    figure
    hold on
    plot(tspk(1:end-1),X,'or','linewidth',2);
    plot(tspk(1:end-1),Z,'og','linewidth',2);
    plot(tspk(1:end-1),M,'ob','linewidth',2);
    plot(tspkc,Xc,'-r','linewidth',1);
    plot(tspkc,Zc,'-g','linewidth',1);
    plot(tspkc,Mc,'-b','linewidth',1);
    plot(tspk,-0.05*ones(1,length(tspk)),'ok','linewidth',1);
    axis([0 Tmax -0.1 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('X_n','Z_n','\DeltaS_n','X_{n,det}','Z_{n,det}','\DeltaS_{n,det}');
    text(20,1.12,['r_{spk}=' num2str(Freq) '  \tau_{dep}=' num2str(taudep) '  \tau_{fac}=' num2str(taufac)],'fontsize',20)
   
    
    
    
   
    
    fprintf('\n');
%     fprintf('mean(M)      = %f\n',mean(M));
%     fprintf('var(M)       = %f\n',var(M));
%     fprintf('var(M)       = %f\n',1000*var(M));
%     fprintf('CV(M)        = %f\n',CVDeltaS);

  %[mean(M) var(M)*1000 CVDeltaS CVspk Mcfp]
  
    CVDeltaS = sqrt(var(M))/mean(M);
    CVX = sqrt(var(X))/mean(X); 
    CVZ = sqrt(var(Z))/mean(Z);
   
    fprintf('CV(DeltaSpk)   = %f\n',CVspk);
    fprintf('mean(DeltaSpk) = %f\n',mean(ISI)-ISImin);
    fprintf('var(DeltaSpk)  = %f\n',var(ISI));
    
%     fprintf('\n')
%     COV = cov(X,Z);
%     fprintf('COV(X,Z):\t%f\t%f\t%f\t%f\n',COV(1,1),COV(1,2),COV(2,1),COV(2,2));
%     Rho = mean((X-mean(X)).*(Z-mean(Z)))/sqrt(var(X)*var(Z));
%     fprintf('Rho(X,Z):\t%f\n',Rho);
    
%     fprintf('\n')
%     fprintf('mean(X) = %f\t mean(Z) = %f\t mean(M) = %f\n',mean(X),mean(Z),mean(M));
%     fprintf('X_{fp}  = %f\t Z_{fp}  = %f\t M_{fp}  = %f\n',Xcfp,Zcfp,Mcfp);
%     fprintf('var(X)  = %f\t var(Z)  = %f\t var(M)  = %f\n',var(X),var(Z),var(M));
%     fprintf('CV(X)   = %f\t CV(Z)   = %f\t CV(M)   = %f\n',CVX,CVZ,CVDeltaS);
    
    
    jaux = find(tspk>1000,1);
    Xaux = X(jaux:end);
    Zaux = Z(jaux:end);
    Maux = M(jaux:end);
    CVDeltaS = sqrt(var(Maux))/mean(Maux);
    CVX = sqrt(var(Xaux))/mean(Xaux); 
    CVZ = sqrt(var(Zaux))/mean(Zaux);
    Rho = mean((X-mean(X)).*(Z-mean(Z)))/sqrt(var(X)*var(Z));
    fprintf('\n')
    fprintf('mean(X) = %f\t mean(Z) = %f\t mean(M) = %f\n',mean(Xaux),mean(Zaux),mean(Maux));
    fprintf('X_{fp}  = %f\t Z_{fp}  = %f\t M_{fp}  = %f\n',Xcfp,Zcfp,Mcfp);
    fprintf('\n')
    fprintf('var(X)  = %f\t var(Z)  = %f\t var(M)  = %f\n',var(Xaux),var(Zaux),var(Maux));
    fprintf('CV(X)   = %f\t CV(Z)   = %f\t CV(M)   = %f\n',CVX,CVZ,CVDeltaS);
    fprintf('sgmad   = %f\t sgmaf   = %f\t sgmadf  = %f\n',sgmadep,sgmafac,sgmadepfac);
    fprintf('\n')
    fprintf('Rho(X,Z):\t%f\n',Rho);
    fprintf('\n')
    
%     mean(X)*mean(Z)
%     (var(X)+mean(X)^2)*(var(Z)+mean(Z)^2)-mean(X)^2*mean(Z)^2

    toff = 3000;
    joff = find(tspk<toff,1,'last');
    tspka = tspk(1:joff);
    Xa = X(1:joff);
    Za = Z(1:joff);
    Ma = M(1:joff);
    joffc = find(tspkc<toff,1,'last');
    tspkca = tspkc(1:joffc);
    Xca = Xc(1:joffc);
    Zca = Zc(1:joffc);
    Mca = Mc(1:joffc);
    
    figure
    hold on
    plot(tspka,Xa,'or','linewidth',2);
    plot(tspka,Za,'og','linewidth',2);
    plot(tspka,Ma,'ob','linewidth',2);
    plot(tspkca,Xca,'-r','linewidth',1);
    plot(tspkca,Zca,'-g','linewidth',1);
    plot(tspkca,Mca,'-b','linewidth',1);
    plot(tspka,-0.05*ones(1,length(tspka)),'ok','linewidth',1);
    axis([0 toff+1000 -0.1 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    text(20,1.12,['f_{spk}=' num2str(Freq) '  \tau_{dep}=' num2str(taudep) '  \tau_{fac}=' num2str(taufac)],'fontsize',20)
    text(3050,0.35,['V(X_n)~' num2str(var(Xaux)*1000,2)],'Color','r','fontsize',20)
    text(3050,0.2,['V(Z_n)~' num2str(var(Zaux)*1000,2)],'Color','g','fontsize',20)
    text(3050,0.05,['V(\DeltaS_n)~' num2str(var(Maux)*1000,2)],'Color','b','fontsize',20)
    legend('X_n','Z_n','\DeltaS_n','X_{n,det}','Z_{n,det}','\DeltaS_{n,det}','NumColumns',2);

    
elseif KSE == 4
    
    % Perturbations of periodic presynaptic input spike trains: normally
    % distributed (random)
    
    PRTYPE = 1;
        % 1: Uniform value of D (variance) for all unperturbed ISIs
        % 2: D (variance) proportional to the unperturbed ISI
        
    
    adep = 0.1;
    afac = 0.2;
    
    taudep = 100;
    taufac = 100;
                
   % Time definitions

    Tmax = 100000;
    dt = 0.01;
    t = 0:dt:Tmax;
    
    % Input spike definitions

    SpkFreqin = 200;
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    tspk = (1:Nspk)*SpkPerin;
    
    if PRTYPE == 1
        D = 1;
        delta = randn(1,Nspk);
        ISIpert = sqrt(D)*delta+SpkPerin;
    elseif PRTYPE == 2
        Dbaseline = 250;
        D = Dbaseline*SpkPerin/1000;
        delta = randn(1,Nspk);
        ISIpert = sqrt(D)*delta+SpkPerin;
    end
    
    tspkpert = zeros(1,length(tspk));
    tspkpert(1) = tspk(1);
    for k=2:length(ISIpert)
        tspkpert(k) = tspkpert(k-1)+ISIpert(k);
    end
    
    
    % Synaptic short-term plasticit Sequences - regular presynaptic spike trains
    
    
    Dtaspkin = diff(tspk);
    X = zeros(1,length(tspk));
    Z = zeros(1,length(tspk));
    X(1) = 1;
    Z(1) = afac;
    for j=1:length(tspk)-1
        X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*exp(-Dtaspkin(j)/taudep);
        Z(j+1) = afac+(1-afac)*(Zinf+(Z(j)-Zinf)*exp(-Dtaspkin(j)/taufac));
    end
    M = X.*Z;
    Pdep =  exp(-Dtaspkin/taudep);
    Pfac =  exp(-Dtaspkin/taufac);
    Xfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
    Zfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);
    
    % Synaptic short-term plasticit Sequences - small perturbations to
    % regular presynaptic spike trains
    
    Dtaspkinpert = ISIpert;
    Xp = zeros(1,length(tspkpert));
    Zp = zeros(1,length(tspkpert));
    Xp(1) = 1;
    Zp(1) = afac;
    for j=1:length(tspkpert)-1
        Xp(j+1) = Xinf+((1-adep)*Xp(j)-Xinf)*exp(-Dtaspkinpert(j)/taudep);
        Zp(j+1) = afac+(1-afac)*(Zinf+(Zp(j)-Zinf)*exp(-Dtaspkinpert(j)/taufac));
    end
    Mp = Xp.*Zp;
    
    
    % Approximation of the Peak sequences by the toy model functions
    % Unperturbed system
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
    E_F = sum((X-F).^2);
    
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
    E_G = sum((Z-G).^2);
    
    % Computation of the toy model functions at the perturbed spike times
    
    F_p = Xfp+(1-Xfp)*exp(-(tspkpert-tspk(1))/sgmadep);
    G_p = Zfp*(1-C*exp(-(tspkpert-tspk(1))/sgmafac));
    
    
    figure
    hold on
    plot(tspk,X,'or','linewidth',2);
    plot(tspk,Z,'og','linewidth',2);
    plot(tspk,M,'ob','linewidth',1);
    plot(tspkpert,Xp,'o','Color',lightcoral,'linewidth',2);
    plot(tspkpert,Zp,'o','Color',mediumacquamarine,'linewidth',2);
    plot(tspkpert,Mp,'o','Color',lightblueish,'linewidth',1);
    axis([0 Tmax 0 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('X_n','Z_n','\Delta S_n');
    
    figure
    hold on
    plot(tspk,X,'-r','linewidth',2);
    plot(tspk,Z,'-g','linewidth',2);
    plot(tspk,M,'-b','linewidth',2);
    plot(tspkpert,Xp,'or','linewidth',1);
    plot(tspkpert,Zp,'og','linewidth',1);
    plot(tspkpert,Mp,'ob','linewidth',1);
    plot(tspkpert,-0.05*ones(1,length(tspkpert)),'ok','linewidth',1);
    axis([0 Tmax -0.1 1.1]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('X_n(0)','Z_n(0)','\DeltaS_n(0)','X_n(\delta_p)','Z_n(\delta_p)','\DeltaS_n(\delta_p)');
    
    
    
    figure
    hold on
    plot(tspk,X,'or','linewidth',2);
    plot(tspk,Z,'og','linewidth',2);
    plot(tspk,M,'ob','linewidth',1);
    plot(tspkpert,Xp,'o','Color',lightcoral,'linewidth',2);
    plot(tspkpert,Zp,'o','Color',mediumacquamarine,'linewidth',2);
    plot(tspkpert,Mp,'o','Color',lightblueish,'linewidth',1);
    axis([0 Tmax 0 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    plot(tspk,F,'-r','linewidth',2);
    plot(tspk,G,'-g','linewidth',2);
    plot(tspkpert,F_p,'-k','linewidth',2);
    plot(tspkpert,G_p,'-k','linewidth',2);
    legend('X_n','Z_n','\Delta S_n');
   
 
 
    figure
    hold on
    plot(tspkpert,F_p-Xp,'or','linewidth',2);
    plot(tspkpert,G_p-Zp,'og','linewidth',2);
    plot(tspkpert,F_p.*G_p-Xp.*Zp,'ob','linewidth',2);
    axis([0 Tmax -0.2 0.2]);
    set(gca,'fontsize',24);
    legend('X_n','Z_n','\Delta S_n');
    
    
    Xc = Xp(floor(length(Xp)/2):floor(length(Xp))); 
    Zc = Zp(floor(length(Zp)/2):floor(length(Zp)));
    Mc = Xc.*Zc;
    
    
    fprintf('\n');
    fprintf('mean(Xc) = %f\t X_{fp} = %f\n',mean(Xc),Xfp);
    fprintf('mean(Zc) = %f\t Z_{fp} = %f\n',mean(Zc),Zfp);
    fprintf('\n');
    COV = cov(Xc,Zc)
    fprintf('\n');
    fprintf('var(Xc) = %f\n',var(Xc));
    fprintf('var(Zc) = %f\n',var(Zc));
    fprintf('var(Xc.*Zc) = %f\n',var(Mc));
    fprintf('\n');
    fprintf('Min(ISIpert) = %f\t Max(ISIpert) = %f\n',min(ISIpert),max(ISIpert)); 
    fprintf('\n');
    fprintf('var(Xc)=%f\t var(Zc)=%f\t var(Mc)=%f\n',var(Xc),var(Zc),var(Mc));
     
    
%     fprintf('COV(X,X) = %f\t COV(X,Z) = %f\n',COV(1,1),COV(1,2));
%     fprintf('COV(Y,X) = %f\t COV(Z,Z) = %f\n',COV(2,1),COV(2,2));

  text(20,1.12,['\Delta_{spk}=' num2str(SpkPerin) '  \tau_{dep}=' num2str(taudep) '  \tau_{fac}=' num2str(taufac)],'fontsize',20)
   
elseif KSE == 5 

    SPKSE = 2;
                % 1: Poisson
                % 2: Poisson (ISI>ISI_min)
                
    taudep = 1000;
    taufac = 1000;
    taudec = 10;           
    
    % Time definitions

    Tmax = 500000;
    dt = 0.01;
    t = 0:dt:Tmax;
                
    % Generation of Poisson-like PWC inputs
    
    
   if SPKSE == 1  
        Freq = 100;
        tr = 1;                                   % trials
        r = Freq/(length(t)-1)*Tmax/1000;          % [rate] = Spk/bin
        spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
        spk = spk(1:length(t));
        tspk = find(spk'>0)*dt;                 % Spike times
        ISI = diff(tspk);    
        CVspk = sqrt(var(ISI))/mean(ISI);
   elseif SPKSE == 2
        ISImin = 1;
        Freq = 200;
        tr = 1;                                 % trials
        r = Freq/(length(t)-1)*Tmax/1000;       % [rate] = Spk/bin
        spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
        spk = spk(1:length(t));
        tspkbase = find(spk'>0)*dt;                 % Spike times
        ISIbase = diff(tspkbase);
        ISI = zeros(1);
        i=0;
        for j=1:length(ISIbase)
            if ISIbase(j) > ISImin
                i=i+1;
                ISI(i) = ISIbase(j);
            end
        end    
        tspk = zeros(1);
        tspk(1) = tspkbase(1);
        for i=2:length(ISI)
            tspk(i) = tspk(i-1)+ISI(i);
        end
        CVspk = sqrt(var(ISI))/(mean(ISI)-ISImin);
   end
   
   tspk = [tspk Tmax+dt];
    
   % Sequences: Poisson spike train inputs
    
    Tspkmax = length(tspk)-1;
    Dtaspkin = diff(tspk(1:end-1));
    X = zeros(1,Tspkmax);
    Z = zeros(1,Tspkmax);
    X(1) = 1;
    Z(1) = afac;
    for j=1:Tspkmax-1
        X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*exp(-Dtaspkin(j)/taudep);
        Z(j+1) = afac+(1-afac)*(Zinf+(Z(j)-Zinf)*exp(-Dtaspkin(j)/taufac));
    end
    M = X.*Z;
    
    S = zeros(1,length(tspk)-1);
    S(1) = afac;
    for j=1:Tspkmax-1
        S(j+1) = exp(-Dtaspkin(j)/taudec)*S(j)+M(j+1);
    end
    
   
    
    % Sequences: Periodic spike train inputs (constant frequency input
    % equal to the Poisson mean rate)
    
    SpkFreqin = Freq;
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    tspkc = (1:Nspk)*SpkPerin;
    
    Dtatspkin = SpkPerin;
    Tspkmax = floor(Tmax/SpkPerin);
    tspkc = SpkPerin:SpkPerin:Tspkmax*SpkPerin;
    jspkc = floor(tspk/dt);

    Pdep =  exp(-Dtatspkin/taudep);
    Pfac =  exp(-Dtatspkin/taufac);
    Xc = zeros(1,Tspkmax);
    Zc = zeros(1,Tspkmax);
    Zaux = zeros(1,Tspkmax);
    Xc(1) = 1;
    Zc(1) = afac;
    for j=1:Tspkmax-1
        Xc(j+1) = Xinf+((1-adep)*Xc(j)-Xinf)*Pdep;       
        Zc(j+1) = (1-afac)*(Zinf+(Zc(j)-Zinf)*Pfac)+afac;
    end
    Mc = Xc.*Zc;
    Xcfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
    Zcfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);
    Mcfp = Xcfp*Zcfp;
    
    Sc = zeros(1,length(tspkc));
    Sc(1) = afac;
    for j=1:length(tspkc)-1
        Sc(j+1) = exp(-Dtatspkin/taudec)*Sc(j)+Mc(j+1);
    end
        

    figure
    hold on
    plot(tspk(1:end-1),X,'or','linewidth',2);
    plot(tspk(1:end-1),Z,'og','linewidth',2);
    plot(tspk(1:end-1),M,'ob','linewidth',2);
    plot(tspkc,Xc,'-r','linewidth',1);
    plot(tspkc,Zc,'-g','linewidth',1);
    plot(tspkc,Mc,'-b','linewidth',1);
    plot(tspk,-0.05*ones(1,length(tspk)),'ok','linewidth',1);
    axis([0 Tmax -0.1 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('X_n','Z_n','\DeltaS_n','X_{n,det}','Z_{n,det}','\DeltaS_{n,det}');
    text(20,1.12,['r_{spk}=' num2str(Freq) '  \tau_{dep}=' num2str(taudep) '  \tau_{fac}=' num2str(taufac)],'fontsize',20)
   
    
    figure
    hold on
    plot(tspk(1:end-1),S,'ob','linewidth',2);
    plot(tspkc,Sc,'o','Color',lightgray,'linewidth',2);
    axis([0 Tmax -0.1 1.2]);
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
    legend('S_n','S_{n,det}');
    
elseif KSE == 301
    
    SPKSE = 2;
    
    % Poisson inputs
    % VarX, VarZ & VarM as a function of the input frequency
    
    taudep = 1000;
    taufac = 1000;

    % VarX, VarZ & VarM as a function of the input frequency

    Tmax = 500000;
    dt = 0.01;
    t = 0:dt:Tmax;

    Freqvec = 10:10:250;
    
    VarXvec = zeros(1,length(Freqvec));
    VarZvec = zeros(1,length(Freqvec));
    VarMvec = zeros(1,length(Freqvec));

    for k=1:length(Freqvec)

        if SPKSE == 1  
            Freq = Freqvec(k);
            tr = 1;                                   % trials
            r = Freq/(length(t)-1)*Tmax/1000;          % [rate] = Spk/bin
            spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
            spk = spk(1:length(t));
            tspk = find(spk'>0)*dt;                 % Spike times
            ISI = diff(tspk);    
            CVspk = sqrt(var(ISI))/mean(ISI);
       elseif SPKSE == 2
            ISImin = 1;
            Freq = Freqvec(k);
            tr = 1;                                 % trials
            r = Freq/(length(t)-1)*Tmax/1000;       % [rate] = Spk/bin
            spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
            spk = spk(1:length(t));
            tspkbase = find(spk'>0)*dt;             % Spike times
            ISIbase = diff(tspkbase);
            ISI = zeros(1);
            i=0;
            for j=1:length(ISIbase)
                if ISIbase(j) > ISImin
                    i=i+1;
                    ISI(i) = ISIbase(j);
                end
            end    
            tspk = zeros(1);
            tspk(1) = tspkbase(1);
            for i=2:length(ISI)
                tspk(i) = tspk(i-1)+ISI(i);
            end
            CVspk = sqrt(var(ISI))/(mean(ISI)-ISImin);
        end

        % Sequences: Poisson spike train inputs

        Tspkmax = length(tspk)-1;
        Dtaspkin = diff(tspk(1:end-1));
        X = zeros(1,Tspkmax);
        Z = zeros(1,Tspkmax);
        X(1) = 1;
        Z(1) = afac;
        for j=1:Tspkmax-1
            X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*exp(-Dtaspkin(j)/taudep);
            Z(j+1) = afac+(1-afac)*(Zinf+(Z(j)-Zinf)*exp(-Dtaspkin(j)/taufac));
        end
        M = X.*Z;

        jaux = find(tspk>1000,1);
        Xaux = X(jaux:end);
        Zaux = Z(jaux:end);
        Maux = M(jaux:end);


        VarXvec(k) = var(Xaux);
        VarZvec(k) = var(Zaux);
        VarMvec(k) = var(Maux);

        [Freq CVspk VarXvec(k) VarZvec(k) VarMvec(k)]
        
        % Sequences: Uniform spike train inputs (constant spiking frequency
        % equal to the Poisson mean firing rate
        
        SpkFreqin = Freq;
        SpkPerin = 1000/SpkFreqin;
        Nspk = floor(Tmax/SpkPerin);
        tspkc = (1:Nspk)*SpkPerin;

        Dtatspkin = SpkPerin;
        Tspkmax = floor(Tmax/SpkPerin);
        tspkc = SpkPerin:SpkPerin:Tspkmax*SpkPerin;
        jspkc = floor(tspk/dt);

        Pdep =  exp(-Dtatspkin/taudep);
        Pfac =  exp(-Dtatspkin/taufac);
        Xc = zeros(1,Tspkmax);
        Zc = zeros(1,Tspkmax);
        Zaux = zeros(1,Tspkmax);
        Xc(1) = 1;
        Zc(1) = afac;
        for j=1:Tspkmax-1
            Xc(j+1) = Xinf+((1-adep)*Xc(j)-Xinf)*Pdep;       
            Zc(j+1) = (1-afac)*(Zinf+(Zc(j)-Zinf)*Pfac)+afac;
        end
        Mc = Xc.*Zc;
        Xcfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
        Zcfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);
        Mcfp = Xcfp*Zcfp;

    end

    figure
    hold on
    plot(Freqvec,VarXvec,'or','linewidth',2);
    plot(Freqvec,VarZvec,'og','linewidth',2);
    plot(Freqvec,VarMvec,'ob','linewidth',2);
    plot(Freqvec,VarXvec,'-r','linewidth',1);
    plot(Freqvec,VarZvec,'-g','linewidth',1);
    plot(Freqvec,VarMvec,'-b','linewidth',1);
    axis([0 taudepvec(end) 0 0.01])
    set(gca,'fontsize',24);
    xlabel('r_{spk}');
    legend('Var(X_n)','Var(Z_n)','Var(\DeltaS_n)')
    text(50,0.0092,['\tau_{dep}=' num2str(taudep) '   \tau_{fac}=' num2str(taufac)],'fontsize',20)

elseif KSE == 302
    
    %   Poisson inputs
    %   VarX, VarZ & VarM as a function of taudep & taufac (taudep=taufac)
    
    SPKSE = 2;
    
    taudepvec = 100:100:2500;
    taufacvec = taudepvec;
    
    Tmax = 500000;
    dt = 0.01;
    t = 0:dt:Tmax;

    if SPKSE == 1  
        Freq = 100;
        tr = 1;                                   % trials
        r = Freq/(length(t)-1)*Tmax/1000;          % [rate] = Spk/bin
        spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
        spk = spk(1:length(t));
        tspk = find(spk'>0)*dt;                 % Spike times
        ISI = diff(tspk);    
        CVspk = sqrt(var(ISI))/mean(ISI);
   elseif SPKSE == 2
        ISImin = 1;
        Freq = 200;
        tr = 1;                                 % trials
        r = Freq/(length(t)-1)*Tmax/1000;       % [rate] = Spk/bin
        spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
        spk = spk(1:length(t));
        tspkbase = find(spk'>0)*dt;             % Spike times
        ISIbase = diff(tspkbase);
        ISI = zeros(1);
        i=0;
        for j=1:length(ISIbase)
            if ISIbase(j) > ISImin
                i=i+1;
                ISI(i) = ISIbase(j);
            end
        end    
        tspk = zeros(1);
        tspk(1) = tspkbase(1);
        for i=2:length(ISI)
            tspk(i) = tspk(i-1)+ISI(i);
        end
        CVspk = sqrt(var(ISI))/(mean(ISI)-ISImin);
    end
    
    VarXvec = zeros(1,length(taudepvec));
    VarZvec = zeros(1,length(taudepvec));
    VarMvec = zeros(1,length(taudepvec));
    
    MeanXvec = zeros(1,length(taudepvec));
    MeanZvec = zeros(1,length(taudepvec));
    MeanMvec = zeros(1,length(taudepvec));
    
    Xcfpvec = zeros(1,length(taudepvec));
    Zcfpvec = zeros(1,length(taudepvec));
    Mcfpvec = zeros(1,length(taudepvec));
    
    
    for k=1:length(taudepvec)
    
        taudep = taudepvec(k);
        taufac = taufacvec(k);
    
        % Sequences: Poisson spike train inputs

        Tspkmax = length(tspk)-1;
        Dtaspkin = diff(tspk(1:end-1));
        X = zeros(1,Tspkmax);
        Z = zeros(1,Tspkmax);
        X(1) = 1;
        Z(1) = afac;
        for j=1:Tspkmax-1
            X(j+1) = Xinf+((1-adep)*X(j)-Xinf)*exp(-Dtaspkin(j)/taudep);
            Z(j+1) = afac+(1-afac)*(Zinf+(Z(j)-Zinf)*exp(-Dtaspkin(j)/taufac));
        end
        M = X.*Z;

        jaux = find(tspk>1000,1);
        Xaux = X(jaux:end);
        Zaux = Z(jaux:end);
        Maux = M(jaux:end);


        VarXvec(k) = var(Xaux);
        VarZvec(k) = var(Zaux);
        VarMvec(k) = var(Maux);
        
        MeanXvec(k) = mean(Xaux);
        MeanZvec(k) = mean(Zaux);
        MeanMvec(k) = mean(Maux);
        
        % Sequences: Uniform spike train inputs (constant spiking frequency
        % equal to the Poisson mean firing rate
        
        SpkFreqin = Freq;
        SpkPerin = 1000/SpkFreqin;
        Nspk = floor(Tmax/SpkPerin);
        tspkc = (1:Nspk)*SpkPerin;

        Dtatspkin = SpkPerin;
        Tspkmax = floor(Tmax/SpkPerin);
        tspkc = SpkPerin:SpkPerin:Tspkmax*SpkPerin;
        jspkc = floor(tspk/dt);

        Pdep =  exp(-Dtatspkin/taudep);
        Pfac =  exp(-Dtatspkin/taufac);
        Xc = zeros(1,Tspkmax);
        Zc = zeros(1,Tspkmax);
        Zaux = zeros(1,Tspkmax);
        Xc(1) = 1;
        Zc(1) = afac;
        for j=1:Tspkmax-1
            Xc(j+1) = Xinf+((1-adep)*Xc(j)-Xinf)*Pdep;       
            Zc(j+1) = (1-afac)*(Zinf+(Zc(j)-Zinf)*Pfac)+afac;
        end
        Mc = Xc.*Zc;
        Xcfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
        Zcfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);
        Mcfp = Xcfp*Zcfp;
        
        [taudep taufac VarXvec(k) VarZvec(k) VarMvec(k)]  
        
        Xcfpvec(k) = Xcfp;
        Zcfpvec(k) = Zcfp;
        Mcfpvec(k) = Mcfp;
        
    end
    
    figure
    hold on
    plot(taudepvec,VarXvec,'or','linewidth',2);
    plot(taudepvec,VarZvec,'og','linewidth',2);
    plot(taudepvec,VarMvec,'ob','linewidth',2);
    plot(taudepvec,VarXvec,'-r','linewidth',1);
    plot(taudepvec,VarZvec,'-g','linewidth',1);
    plot(taudepvec,VarMvec,'-b','linewidth',1);
    set(gca,'fontsize',24);
    axis([0 taudepvec(end) 0 0.01])
    xlabel('\tau_{dep} = \tau_{fac}');
    legend('Var(X_n)','Var(Z_n)','Var(\DeltaS_n)')
    text(80,0.0092,['r_{spk}=' num2str(Freq)],'fontsize',20)
    
    figure
    hold on
    plot(taudepvec,MeanXvec,'or','linewidth',2);
    plot(taudepvec,MeanZvec,'og','linewidth',2);
    plot(taudepvec,MeanMvec,'ob','linewidth',2);
    plot(taudepvec,MeanXvec,'-r','linewidth',1);
    plot(taudepvec,MeanZvec,'-g','linewidth',1);
    plot(taudepvec,MeanMvec,'-b','linewidth',1);
    plot(taudepvec,Xcfpvec,'--r','linewidth',1);
    plot(taudepvec,Zcfpvec,'--g','linewidth',1);
    plot(taudepvec,Mcfpvec,'--b','linewidth',1);
    set(gca,'fontsize',24);
    axis([0 taudepvec(end) 0 1])
    xlabel('\tau_{dep} = \tau_{fac}');
    legend('Mean(X_n)','Mean(Z_n)','Mean(\DeltaS_n)')
    text(80,0.92,['r_{spk}=' num2str(Freq)],'fontsize',20)
   
end

