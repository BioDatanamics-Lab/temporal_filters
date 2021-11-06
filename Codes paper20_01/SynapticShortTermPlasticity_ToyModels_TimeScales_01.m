% 2020-01-18

% Synaptic depression and facilitation models 
% (Dayan & Abbott, see Ermentrout & Terman book, Chapter 7) 

% Inherited from "SynapticShortTermPlasticity_ToyModels.m"

% Computes the recurrent sequences for X, Z and DeltaS=X*Z
% Analitical estimation of the parameters for the approximating functions
% F, G and H = F*G, in particular sgma_d and sgma_f (long-term depression
% and facilitation time scales.
% Expansion of H in sume of exponents and approximation of H by Hcut
% consisting of the three first terms of the sum of exponents disregarding
% the term combining sgma_d and sigma_fac, and leaving only the constant
% term and the terms including sgmadep (sgma_dep) and sgmafac (sgma_fac)
% separately.

% clearvars;
% close all;

KSE = 2;
           % 1: Single input frequency and set of parameters
           % 2: Range of input frequencies for representative values of 
           %    taudep and taufac (all other parameters fixed)

lightblueish = [.4 .6 .9];
lightcoral = [0.94 0.5 0.5];
lightgray = [.7 .7 .7];

% Parameters

taudec = 5;
taudep = 300;
taufac = 300;
adep = 0.1;
afac = 0.2;
Xinf = 1;
Zinf = 0;

% Time definitions

Tmax = 10000;
dt = 0.01;
t = 0:dt:Tmax;

if KSE == 1
    
    % Input spike definitions

    SpkFreqin = 80;
    SpkPerin = 1000/SpkFreqin;
    Nspk = floor(Tmax/SpkPerin);
    

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
    DeltaS = X.*Z;
    jspk = floor(tspk/dt);
    Xfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
    Zfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);

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

    H = F.*G;

    Hexpanded = Zfp*(Xfp-Xfp*C*exp(-(tspk-tspk(1))/sgmafac)+(1-Xfp)*exp(-(tspk-tspk(1))/sgmadep)-(1-Xfp)*C*exp(-(tspk-tspk(1))*(1/sgmadep+1/sgmafac)));
    Hcut = Zfp*(Xfp-Xfp*C*exp(-(tspk-tspk(1))/sgmafac)+(1-Xfp)*exp(-(tspk-tspk(1))/sgmadep));

    [taudep sgmadep taufac sgmafac E_F E_G]

    figure
    hold on
    plot(tspk,F,'-r','linewidth',2);
    plot(tspk,G,'-g','linewidth',2);
    plot(tspk,H,'-b','linewidth',2);
    plot(tspk,X,'or','linewidth',1);
    plot(tspk,Z,'og','linewidth',1);
    plot(tspk,DeltaS,'ob','linewidth',1);
    axis([0 Tmax 0  1.1]);
    legend('X_n','Z_n','\DeltaS_n')
    set(gca,'fontsize',24);
    xlabel('t  [ms]');
elseif KSE == 2
    
    % Input spike definitions

    SpkFreqinmin = 10;
    SpkFreqinmax = 200;
    DeltaSpkFreqin = 10;
    SpkFreqin = SpkFreqinmin:DeltaSpkFreqin:SpkFreqinmax;
    
    sgmadep = zeros(1,length(SpkFreqin));
    sgmafac = zeros(1,length(SpkFreqin));
    E_F = zeros(1,length(SpkFreqin));
    E_G = zeros(1,length(SpkFreqin));
    
    auxerror = 0;
    
    for l=1:length(SpkFreqin)
        
        SpkPerin = 1000/SpkFreqin(l);
        Nspk = floor(Tmax/SpkPerin);
        
        % Sequence of peaks computed from the analytical solutions to the
        % differential equations using the previous sequence values to generate the
        % current ones (instead of the previous values in the differential
        % equation discretization). 

        Dtatspkin = SpkPerin;
        Tspkmax = floor(Tmax/SpkPerin);
        tspk = SpkPerin:SpkPerin:Tspkmax*SpkPerin;



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
        DeltaS = X.*Z;
        jspk = floor(tspk/dt);
        Xfp = (1-Pdep)*Xinf/(1-(1-adep)*Pdep);
        Zfp = ((1-Pfac)*(1-afac)*Zinf+afac)/(1-(1-afac)*Pfac);

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
        sgmadep(l) = (tspk(kaux)-tspk(1))/log((1-Xfp)/(X(kaux)-Xfp));
        F = Xfp+(1-Xfp)*exp(-(tspk-tspk(1))/sgmadep(l));
        E_F(l) = sum((X-F).^2);

        Gcut = 0.5*(Zfp-Z(1))+Z(1);
        for k=2:length(tspk)
            if Z(k)>Gcut
                kaux = k;    
                break
            end
        end

        C = 1-afac/Zfp;
        sgmafac(l) = (tspk(kaux)-tspk(1))/log(Zfp*C/(Zfp-Z(kaux)));
        G = Zfp*(1-C*exp(-(tspk-tspk(1))/sgmafac(l)));
        E_G(l) = sum((Z-G).^2);

        if E_F(l) > 1.0e-20
            fprintf('E_F too big \n')
            auxerror = 1;
        end
        if E_G(l) > 1.0e-20
            fprintf('E_G too big \n')     
            auxerror = 1;
        end
        
    end
    
    sgmadepfac = 1./(1./sgmadep+1./sgmafac);
    
    [taudep taufac adep afac]
    
    [SpkFreqin' sgmadep' sgmafac' E_F' E_G']
    
    if auxerror == 1
        fprintf('E_F or E_G too big: check! \n')
    end
    
    
    figure
    hold on
    plot(SpkFreqin,sgmadep,'ob','linewidth',2);    
    plot(SpkFreqin,sgmadep,'-b','linewidth',1);    
    axis([0 SpkFreqinmax 0  210]);    
    set(gca,'fontsize',24);
    xlabel('f_{spk}  [Hz]');
    ylabel('\sigma_d  [ms]');
    
    figure
    hold on
    plot(SpkFreqin,sgmafac,'ob','linewidth',2);    
    plot(SpkFreqin,sgmafac,'-b','linewidth',1);    
    axis([0 SpkFreqinmax 0  210]);    
    set(gca,'fontsize',24);
    xlabel('f_{spk}  [Hz]');
    ylabel('\sigma_f  [ms]');
    
    figure
    hold on
    plot(SpkFreqin,sgmadepfac,'ob','linewidth',2);    
    plot(SpkFreqin,sgmadepfac,'-b','linewidth',1);    
    axis([0 SpkFreqinmax 0  210]);    
    set(gca,'fontsize',24);
    xlabel('f_{spk}  [Hz]');
    ylabel('\sigma_{d+f}  [ms]');
    
    figure
    hold on
    plot(SpkFreqin,sgmadep,'or','linewidth',2);
    plot(SpkFreqin,sgmafac,'og','linewidth',2);
    plot(SpkFreqin,sgmadepfac,'ob','linewidth',2);
    plot(SpkFreqin,sgmadep,'-r','linewidth',1);
    plot(SpkFreqin,sgmafac,'-g','linewidth',1);
    plot(SpkFreqin,sgmadepfac,'-b','linewidth',1);
    axis([0 SpkFreqinmax 0  210]);
    legend('\sigma_{d}','\sigma_f','\sigma_{d+f}')
    set(gca,'fontsize',24);
    xlabel('f_{spk}  [Hz]');
    ylabel('[ms]');
    
end


