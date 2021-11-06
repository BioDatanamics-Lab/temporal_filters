clear all; close all; clc;
f1 = figure(1); set(gcf,'Position', [10 10 1450 700])
tFinForError = 1000;

axFnt = 12; ttlFnt = 23; lblFnt = 14;

%%% Fitting High Pass Temporal Filters
figure(1);
ODEparams.tau_fac = 500;
dp = [0 50]; dpIx = length(dp);
frq = [10:10:150]; frqIx = length(frq);
for i = 1:dpIx
    ODEparams.tau_dep = dp(i);
    for j = 1:frqIx
        
        % setting simulation lengths
		spikePeriod = 1000./frq(j);

		% Modify Simulation Length Depending on Frequency
		tFin = 1000;
	 	% Other settings.
		ODEparams.C = 1; ODEparams.E_L = -60; ODEparams.g_L = 1; ODEparams.I_app = 0; ODEparams.G_ex = .5; ODEparams.E_ex = 0; ODEparams.tau_dec = 3; ODEparams.U = .1;

		%%%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[rawSig.sig, rawSig.time, rawSig.isFired] = numerical(ODEparams, spikePeriod, tFin);
        
        isDeltaS = 0; thresh = .01;
		%%%%% TEMPORAL FILTER ANALYZER %%%%%%%%%%%%%%%%%%%%%%%
		dmy = tfCompute(rawSig, thresh, 3, tFinForError,ODEparams, spikePeriod, tFin, isDeltaS);
       
		if dmy.class == 2
			% 1 - cvx fitted value of sigmaD; 2 - MSE; 3 - computed sigmaD
            sigmaD(i,j,4) = dmy.class;
			sigmaD(i,j,1) = dmy.sigma1.cvx;
			sigmaD(i,j,2) = dmy.sigma1.MSE^.5;
            sigmaD(i,j,3) = dmy.sigma1.analytic;
			%{
		    %checking plots
		    figure;
		    [pk, pkLoc] = findpeaks(rawSig.sig(3,:),rawSig.time);
		    fitS = pk(end)-(pk(end) - pk(1))*exp(-pkLoc/dmy.sigma1.cvx);
		    plot(rawSig.time, rawSig.sig(3,:)); hold on; plot(pkLoc, pk); plot(pkLoc, fitS);
		    legend('signal','pk trace','pk fit'); %xlim([0 tFinForError]);
		    title(num2str(frq(j)));
			%}
		else
            sigmaD(i,j,4) = dmy.class;
			sigmaD(i,j,1) = NaN;
			sigmaD(i,j,2) = NaN;
            sigmaD(i,j,3) = NaN;
        end

        if frq(j) == 40
            if i == 1
                h = subplot(2,3,2);
				ttlStr = 'A2';
            else
                h = subplot(2,3,3);
				ttlStr = 'A3';
            end
            [pk, pkLoc] = findpeaks(rawSig.sig(3,:),rawSig.time);
            fitS = pk(end)-(pk(end) - pk(1))*exp(-pkLoc/dmy.sigma1.cvx);
            plot(rawSig.time, rawSig.sig(3,:),'-b'); hold on; plot(pkLoc, fitS,'-k');
            h.FontSize = axFnt;
			ttl = title(h,ttlStr,'FontSize',ttlFnt);
			ttl.Units = 'Normalize'; 
			ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
			ttl.HorizontalAlignment = 'left';
			ylim([0 1]); xlim([0 tFinForError]);
			xlabel('t [ms]','FontSize',lblFnt); ylabel('S [units]','FontSize',lblFnt);
			legend('Synaptic Response','Fitted Temporal Filter','FontSize',axFnt)
        end
        

    end
end
k = subplot(2,3,1);
plot(k, frq, squeeze(sigmaD(1,:,1)), '-bo','MarkerSize',3); hold on;
plot(k, frq, squeeze(sigmaD(2,:,1)), '-ro','MarkerSize',3);
k.FontSize = axFnt;
ttl = title(k,'A1','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
legend('\tau_{dep} = 0','\tau_{dep} = 50','FontSize',axFnt);
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\sigma_{f,S}','FontSize',lblFnt);

'high pass error max'
max(max(sigmaD(:,:,2)))

%%% Fitting low Pass Temporal Filters
ODEparams.tau_dep = 500;
fc = [0 50]; dpIx = length(fc);
frq = [10:10:150]; frqIx = length(frq);
for i = 1:dpIx
    ODEparams.tau_fac = fc(i);
    for j = 1:frqIx
        
        % setting simulation lengths
		spikePeriod = 1000./frq(j);

		% Modify Simulation Length Depending on Frequency
		tFin = 500;
	 	% Other settings.
		ODEparams.C = 1; ODEparams.E_L = -60; ODEparams.g_L = 1; ODEparams.I_app = 0; ODEparams.G_ex = .5; ODEparams.E_ex = 0; ODEparams.tau_dec = 3; ODEparams.U = .1;

		%%%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[rawSig.sig, rawSig.time, rawSig.isFired] = numerical(ODEparams, spikePeriod, tFin);
        
        isDeltaS = 0; thresh = .01;
		%%%%% TEMPORAL FILTER ANALYZER %%%%%%%%%%%%%%%%%%%%%%%
		dmy = tfCompute(rawSig, thresh, 3, tFinForError,ODEparams, spikePeriod, tFin, isDeltaS);
       
		if dmy.class == 4
			% 1 - cvx fitted value of sigmaD; 2 - MSE; 3 - computed sigmaD
            sigmaD(i,j,4) = dmy.class;
			sigmaD(i,j,1) = dmy.sigma1.cvx;
			sigmaD(i,j,2) = dmy.sigma1.MSE^.5;
            sigmaD(i,j,3) = dmy.sigma1.analytic;
			%{
		    %checking plots
		    figure;
		    [pk, pkLoc] = findpeaks(rawSig.sig(3,:),rawSig.time);
		    fitS = pk(end)-(pk(end) - pk(1))*exp(-pkLoc/dmy.sigma1.cvx);
		    plot(rawSig.time, rawSig.sig(3,:)); hold on; plot(pkLoc, pk); plot(pkLoc, fitS);
		    legend('signal','pk trace','pk fit'); % xlim([0 tFinForError]);
		    title(num2str(frq(j)));
			%}
		else
            sigmaD(i,j,4) = dmy.class;
			sigmaD(i,j,1) = NaN;
			sigmaD(i,j,2) = NaN;
            sigmaD(i,j,3) = NaN;
        end

        if frq(j) == 20
            if i == 1
                h = subplot(2,3,5);
				ttlStr = 'B2';
            else
                h = subplot(2,3,6);
				ttlStr = 'B3';
            end
            [pk, pkLoc] = findpeaks(rawSig.sig(3,:),rawSig.time);
            fitS = pk(end)-(pk(end) - pk(1))*exp(-pkLoc/dmy.sigma1.cvx);
            plot(rawSig.time, rawSig.sig(3,:),'-b'); hold on; plot(pkLoc, fitS,'-k');
            h.FontSize = axFnt;
			ttl = title(h,ttlStr,'FontSize',ttlFnt);
			ttl.Units = 'Normalize'; 
			ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
			ttl.HorizontalAlignment = 'left';
			ylim([0 .3]); xlim([0 tFinForError]);
			xlabel('t [ms]','FontSize',lblFnt); ylabel('S [units]','FontSize',lblFnt);
			legend('Synaptic Response','Fitted Temporal Filter','FontSize',axFnt)
        end
    end
end
k = subplot(2,3,4);
plot(frq, squeeze(sigmaD(1,:,1)), '-bo','MarkerSize',3); hold on;
plot(frq, squeeze(sigmaD(2,:,1)), '-ro','MarkerSize',3);
k.FontSize = axFnt;
ttl = title(k,'B1','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';
legend('\tau_{fac} = 0','\tau_{fac} = 50','FontSize',axFnt);
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\sigma_{d,S}','FontSize',lblFnt);

'low pass error max'
max(max(sigmaD(:,:,2)))

saveas(f1,'../img/fig_Squant1','epsc');

