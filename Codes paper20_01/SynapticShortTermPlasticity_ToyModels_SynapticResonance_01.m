% 2020-01-08

% Synaptic depression and facilitation models 
% (Dayan & Abbott, see Ermentrout & Terman book, Chapter 7) 

% clearvars;
% close all;

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightgray = [.7 .7 .7];

% Parameters

taudec = 5;
taudep = 100;
taufac = 100;
adep = 0.1;
afac = 0.2;
Xinf = 1;
Zinf = 0;

KSE = 2;
           % 1: Single presynaptic input frequency. Temporal evolutions
           % 2: Steady-state, frequency-dependent response 

if KSE == 1          
           
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
    xlabel('t  [ms]');

    figure
    hold on
    plot(tspk,F,'-r','linewidth',2);
    plot(tspk,G,'-g','linewidth',2);
    plot(tspk,H,'-b','linewidth',2);
    plot(tspk,X,'or','linewidth',1);
    plot(tspk,Z,'og','linewidth',1);
    plot(tspk,M,'ob','linewidth',1);
    axis([0 Tmax 0  1.1]);
    legend('X_n','Z_n','\DeltaS_n')
    set(gca,'fontsize',24);
    xlabel('t  [ms]');

    
elseif KSE == 2
    
    % Frequency definitions
    
    SpkFreqinmin = 0;
    SpkFreqinmax = 1000;
    DeltaSpkFreqin = 0.1;
    SpkFreqin = SpkFreqinmin:DeltaSpkFreqin:SpkFreqinmax;
    SpkFreqin(1) = 0.1;
    SpkPerin = 1000./SpkFreqin;
    Dtatspkin = SpkPerin;
    
    
    % Frequency-dependent steady-state peak envelope responses
    
   
              
    Pdep =  exp(-Dtatspkin/taudep);
    Pfac =  exp(-Dtatspkin/taufac);

    Xfp = (1-Pdep)*Xinf./(1-(1-adep)*Pdep);
    Zfp = ((1-Pfac)*(1-afac)*Zinf+afac)./(1-(1-afac)*Pfac);
    DeltaSfp = Xfp.*Zfp;
    
    
    
    
    % Approximation of Xfp and Zfp by the sigmoid function of the
    % ln(SpkFreqin) (rescaling the input frequency by its natural
    % logarithm). 
    
    f = log(SpkFreqin);
    for j=1:length(SpkFreqin)
        if Xfp(j)<=Xinf/2;
            fhlf_d = SpkFreqin(j);
            jaux = j;
            break;          
        end
    end
    Xfp1 = (Xfp(jaux+1)-Xfp(jaux-1))/(f(jaux+1)-f(jaux-1));
    fslp_d = -1/(4*Xfp1);
    F = Xinf./(1+exp((log(SpkFreqin)-log(fhlf_d))/fslp_d));
    E_F = sum((F-Xfp).^2);
    Faux = 1./(1+exp((f-log(fhlf_d))/fslp_d)); 
    for j=1:length(SpkFreqin)
        if Zfp(j)>=afac+(1-afac)/2;
            fhlf_f = SpkFreqin(j);
            jaux = j;
            break;
        end
    end
    Zfp1 = (Zfp(jaux+1)-Zfp(jaux-1))/(f(jaux+1)-f(jaux-1));
    fslp_f = (1-afac)/(4*Zfp1);
    G = afac+(1-afac)./(1+exp(-(log(SpkFreqin)-log(fhlf_f))/fslp_f));
    E_G = sum((G-Zfp).^2);
    Gaux = afac+(1-afac)./(1+exp(-(f-log(fhlf_f))/fslp_f));
    
    % Computation of the characteristic frequencies using Xfp and Zfp
    decaux = 0.37;
    betaux = (1-afac)*Zinf;
    kappa_d = 1000/(taudep*log((1-decaux*(1-adep))/(1-decaux)));
    kappa_f = 1000/(-taufac*log((1-betaux-afac)*(1-decaux)/(((1-betaux-afac)*(1-decaux)+betaux+afac)*(1-afac)-betaux)))
    
    % Computation of the characteristic frequencies using the approximating
    % functions F and G.
    
    kappa_d_app = exp(log(fhlf_d))*(Xinf/decaux-1)^fslp_d;
    kappa_f_app = exp(log(fhlf_f))*((1-decaux)/decaux)^fslp_f;
    
    [kappa_d kappa_f kappa_d_app kappa_f_app E_F E_G]
    
    figure
    hold on
    plot(f,Xfp,'or','linewidth',2);
    plot(f,Zfp,'og','linewidth',2);
    plot(f,Faux,'--k','linewidth',1)
    plot(f,Gaux,'--k','linewidth',1)
    set(gca,'fontsize',24);
    
   
    figure
    hold on
    plot(SpkFreqin,Xfp,'-r','linewidth',2);
    plot(SpkFreqin,Zfp,'-g','linewidth',2);
    plot(SpkFreqin,DeltaSfp,'-b','linewidth',2);
%     plot(SpkFreqin,F,'--r','linewidth',1);
%     plot(SpkFreqin,G,'--g','linewidth',1);
%     plot(SpkFreqin,F.*G,'--b','linewidth',1);
    axis([0 SpkFreqinmax 0  1]);
    h = legend('$\bar{X}$','$\bar{Z}$','$\bar{\Delta S}$');
    set(h,'interpreter','Latex')        
    set(gca,'fontsize',24);
    xlabel('f_{spk}  [Hz]');
    ylabel('');
    
    
    
    
end
