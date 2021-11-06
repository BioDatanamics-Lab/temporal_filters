clear all; close all; clc

f1 = figure(1); set(gcf,'Position', [10 10 1450 650])
axFnt = 12; ttlFnt = 20; lblFnt = 13;

%% \Delta S TFs MT Model
dep = [150 1 150 150 50]; fac = [1 150 150 50 150];
frq = [1 20 100 150];


clr = {'--k','r','g','b'};
tFinForError = 301; tFin = 500;
ODEparams.tau_dec = .01; % kill no summation
ODEparams.C = 1; ODEparams.E_L = -60; ODEparams.g_L = 1; ODEparams.I_app = 0; ODEparams.G_ex = 1; ODEparams.E_ex = 0; ODEparams.U = .1;

for k = 1:5 %dep/fac pairs
    h4 = subplot(2,3,k);
    ODEparams.tau_dep = dep(k); ODEparams.tau_fac = fac(k);
    for j = 1:length(frq)

        % setting simulation lengths
        spikePeriod = 1000./frq(j);

        %%%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [rawSig.sig, rawSig.time, rawSig.isFired] = numerical(ODEparams, spikePeriod, tFin);        
        [pk, pkLoc] = findpeaks(rawSig.sig(3,:),rawSig.time);
        plot(pkLoc,pk,clr{j},'Linewidth',1.3); hold on;

    end
    xlabel('time [ms]','FontSize',lblFnt); ylabel('\Delta S [units]','FontSize',lblFnt);
    xlim([0 300]); %ylim([0 .4]);
    ttl = title(strcat('A',num2str(k)),'FontSize',ttlFnt);
    ttl.Units = 'Normalize'; 
    ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
    ttl.HorizontalAlignment = 'left';
    if k == 1
        ylim([0 .33])
        legend({strcat('\Delta S (',num2str(frq(1)),' Hz)'),...
           strcat('\Delta S (',num2str(frq(2)),' Hz)'), ...
           strcat('\Delta S (',num2str(frq(3)),' Hz)'), ...
           strcat('\Delta S (',num2str(frq(4)),' Hz)')}, ...
            'FontSize',axFnt,'Location','northeast','NumColumns',2)
    end
end

%saveas(f1,'../img/TF_DeltaS_types','epsc');

pause; close

%%%% add summation // S TFs
f1 = figure(1); set(gcf,'Position', [10 10 1450 350])
dep = [150 1 150]; fac = [1 50 150];
frq = [100]; dec = [.01 5 15];

clr = {'--k','r','b'};
tFinForError = 301; tFin = 500; 
ODEparams.C = 1; ODEparams.E_L = -60; ODEparams.g_L = 1; ODEparams.I_app = 0; ODEparams.G_ex = 1; ODEparams.E_ex = 0; ODEparams.U = .1;

for k = 1:3 %dep fac pairs
    hAx(k)=subplot(1,3,k); 
    ODEparams.tau_dep = dep(k); ODEparams.tau_fac = fac(k);
    for j = 1:length(dec)
        ODEparams.tau_dec = dec(j); % synthetic -- no summation
        
        % setting simulation lengths
        spikePeriod = 1000./frq(1);

        %%%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [rawSig.sig, rawSig.time, rawSig.isFired] = numerical(ODEparams, spikePeriod, tFin);
        [pk, pkLoc] = findpeaks(rawSig.sig(3,:),rawSig.time);
        plot(pkLoc,pk,clr{j},'Linewidth',1.3); hold on;

    end
    xlabel('time [ms]','FontSize',lblFnt); ylabel('\Delta S [units]','FontSize',lblFnt); %ylim([0 1]);
    hAx(k).Position(4)=hAx(k).Position(4)*0.94;
    dmy = ylim; xlim([0 300]);
    hTtl(k)=title(strcat('A',num2str(k)),'FontSize',ttlFnt,'Position',[0 dmy(2)],'HorizontalAlignment','left');
    hAx(k).Position(2)=hAx(k).Position(2)+hAx(k).Position(4)*0.05;
    
    if k == 1
        legend({'\Delta S ',...
           strcat('S (\tau_{dec} = ',num2str(dec(2)),')'), ...
           strcat('S (\tau_{dec} = ',num2str(dec(3)),')')}, ...
            'FontSize',axFnt,'Location','northeast','NumColumns',1)
    end
    hAx=gca;
    hAx.Position=hAx.Position.*[1 1 1 0.95];
end

%saveas(f1,'../img/TF_S_types','epsc');

pause; close

% post synaptic summation // Voltage TFs
f1 = figure(1); set(gcf,'Position', [10 10 1450 350])
dep = [150 1 150]; fac = [1 50 150];
frq = [100]; dec = 5;
gL = [0.1, 0.5, 1];

clr = {'r','b','g'};
tFinForError = 301; tFin = 500; 
ODEparams.C = 1; ODEparams.E_L = -60; ODEparams.I_app = 0; ODEparams.G_ex = 1; ODEparams.E_ex = 0; ODEparams.U = .1;

for k = 1:3 %dep fac pairs
    h4 = subplot(1,3,k);
    ODEparams.tau_dep = dep(k); ODEparams.tau_fac = fac(k); ODEparams.tau_dec = dec(1);
    for j = 1:length(gL)
        ODEparams.g_L = gL(j);
        
        % setting simulation lengths
        spikePeriod = 1000./frq(1);

        %%%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [rawSig.sig, rawSig.time, rawSig.isFired] = numerical(ODEparams, spikePeriod, tFin);
        [pk, pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
        plot(pkLoc,pk,clr{j},'Linewidth',1.3); hold on;

    end
    [pk, pkLoc] = findpeaks(rawSig.sig(3,:),rawSig.time);
    plot(pkLoc,(pk*60)-60,'linestyle','--','color',[.5 .5 .5]); hold on;
    xlabel('time [ms]','FontSize',lblFnt); ylabel('V [mV]','FontSize',lblFnt);
    xlim([0 300]); %ylim([0 .4]);

    if k == 1
        ylim([-60 -20]);
        legend({
           strcat('V (\tau_{mem} = ',num2str(1/gL(1)),')'), ...
           strcat('V (\tau_{mem} = ',num2str(1/gL(2)),')'), ...
           strcat('V (\tau_{mem} = ',num2str(1/gL(3)),')'),'S (rescaled)'}, ...
            'FontSize',axFnt-1,'Location','northeast','NumColumns',2)
    end
    h4.Position(4)=h4.Position(4)*0.94;
    dmy = ylim;
    hTtl=title(strcat('A',num2str(k)),'FontSize',ttlFnt,'Position',[0 dmy(2)],'HorizontalAlignment','left');
    h4.Position(2)=h4.Position(2)+h4.Position(4)*0.05;end

%saveas(f1,'../img/TF_V_types','epsc');
