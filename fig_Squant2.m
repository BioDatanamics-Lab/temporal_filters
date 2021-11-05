clear all; close all; clc;
ODEparams.tau_dec = 3; ODEparams.U = .1;
f1 = figure(1); set(gcf,'Position', [10 10 1450 700])
tFinForError = 301;

axFnt = 12; ttlFnt = 23; lblFnt = 14;

%% USED IN PART B and markers in A
valPlt = [200]; frqPlt = [30 90 150]; %frqPlt2 = [40 100];

%%% Note: A1 and A2 TF time scales are also obtained in A3's computation.
%%% A1 and A2 can also be plotted using information obtained in A3's computation.
%%% Errors reported in A3's computation verify this.

strtPlotFrq = 2;
%%% Extract low pass component
figure(1);
h(1) = subplot(2,3,1);
ODEparams.tau_fac = 0;
dp = [100 200 300]; dpIx = length(dp);
frq = [10:10:150]; frqIx = length(frq);
for i = 1:dpIx
    ODEparams.tau_dep = dp(i);
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
			sigmaD(i,j,1) = dmy.sigma1.cvx;
			sigmaD(i,j,2) = dmy.sigma1.MSE^.5;
            sigmaD(i,j,3) = dmy.sigma1.analytic;
		else
			'Stalling computation: expecting decay (class 4) instead getting a different class:'
			dmy.class
			keyboard
		end
        
    end
end
plot(h(1), frq(strtPlotFrq:end), squeeze(sigmaD(1,strtPlotFrq:end,1)), '-bo','MarkerSize',3); hold on;
plot(h(1), frq(strtPlotFrq:end), squeeze(sigmaD(2,strtPlotFrq:end,1)), '-ro','MarkerSize',3);
plot(h(1), frq(strtPlotFrq:end), squeeze(sigmaD(3,strtPlotFrq:end,1)), '-go','MarkerSize',3);
dmy = [squeeze(sigmaD(2,3,1)) squeeze(sigmaD(2,9,1)) squeeze(sigmaD(2,15,1))];
plot(h(1), frqPlt(1), dmy(1), 'ko','MarkerSize',5,'LineWidth',2);
plot(h(1), frqPlt(2), dmy(2), 'ks','MarkerSize',5,'LineWidth',2);
plot(h(1), frqPlt(3), dmy(3), 'kd','MarkerSize',5,'LineWidth',2);
h(1).FontSize = axFnt;
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\sigma_{d,S} [ms]','FontSize',lblFnt); ylim([0 250]); xlim([0 153]);
legend('\tau_{dep} = 100','\tau_{dep} = 200','\tau_{dep} = 300','FontSize',axFnt)
ttl = title('A1','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

'max error on low pass'
max(max(sigmaD(:,:,2)))

%%% Extract high pass component
figure(1);
h(2) = subplot(2,3,2);
ODEparams.tau_dep = 0;
fc = [100 200 300]; fcIx = length(fc);
frq = [10:10:150]; frqIx = length(frq);
for i = 1:fcIx
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

		if dmy.class == 2
			% 1 - cvx fitted value of sigmaF; 2 - MSE; 3 - computed sigmaF
			sigmaF(i,j,1) = dmy.sigma1.cvx;
			sigmaF(i,j,2) = dmy.sigma1.MSE^.5;
            sigmaF(i,j,3) = dmy.sigma1.analytic;
		else
			'Stalling computation: expecting rising (class 2) instead getting a different class:'
			dmy.class
			keyboard
		end

    end
end
plot(h(2), frq(strtPlotFrq:end), squeeze(sigmaF(1,strtPlotFrq:end,1)), '-bo','MarkerSize',3); hold on;
plot(h(2), frq(strtPlotFrq:end), squeeze(sigmaF(2,strtPlotFrq:end,1)), '-ro','MarkerSize',3);
plot(h(2), frq(strtPlotFrq:end), squeeze(sigmaF(3,strtPlotFrq:end,1)), '-go','MarkerSize',3);
dmy = [squeeze(sigmaF(2,3,1)) squeeze(sigmaF(2,9,1)) squeeze(sigmaF(2,15,1))];
plot(h(2), frqPlt(1), dmy(1), 'ko','MarkerSize',5,'LineWidth',2);
plot(h(2), frqPlt(2), dmy(2), 'ks','MarkerSize',5,'LineWidth',2);
plot(h(2), frqPlt(3), dmy(3), 'kd','MarkerSize',5,'LineWidth',2);
h(2).FontSize = axFnt;
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\sigma_{f,S} [ms]','FontSize',lblFnt); ylim([0 250]); xlim([0 153]);
legend('\tau_{fac} = 100','\tau_{fac} = 200','\tau_{fac} = 300','FontSize',axFnt)
ttl = title('A2','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

'max error on high pass'
max(max(sigmaF(:,:,2)))

%%% A3 - Fitting Band Pass Temporal Filters 
figure(1);
h(3) = subplot(2,3,3);
vals = [100 200 300]; valsIx = length(vals);
%vals = [100]; valsIx = length(vals);
frq = [10:10:150]; frqIx = length(frq);
%frq = [60]; frqIx = length(frq);
cnt = 3;
for i = 1:valsIx
	for kkk = (1:1:frqIx)

		% setting simulation lengths
		spikePeriod = 1000./frq(kkk);

		% Modify Simulation Length Depending on Frequency
		tFin = 500;
		ODEparams.tau_dep = vals(i); ODEparams.tau_fac = vals(i);
	 	% Other settings.
		ODEparams.C = 1; ODEparams.E_L = -60; ODEparams.g_L = 1; ODEparams.I_app = 0; ODEparams.G_ex = .5; ODEparams.E_ex = 0; ODEparams.tau_dec = 3; ODEparams.U = .1;

		%%%%% SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		[rawSig.sig, rawSig.time, rawSig.isFired] = numerical(ODEparams, spikePeriod, tFin);

		%% For figure 3 plotting later
		%% 1 - R- ; 2- u+
		%% soln_r(3).peaks - soln_r(3).traughs
		pkR = findpeaks(-1*rawSig.sig(1,:))*-1; ssHold(i,kkk,1) = pkR(end);
		pku = findpeaks(rawSig.sig(2,:)); ssHold(i,kkk,2) = pku(end);
	   
        
        isDeltaS = 0; thresh = .01;
		%%%%% TEMPORAL FILTER ANALYZER %%%%%%%%%%%%%%%%%%%%%%%
		dmy = tfCompute(rawSig, thresh, 3, tFinForError,ODEparams, spikePeriod, tFin, isDeltaS);

		% 1- sigmaR, 2- sigmaU, 3-sigmaR+U
        % // row 1: value;  row 2: err; row 3: analytic. 
		simgaStore(i,kkk,1,1) = dmy.sigma1.cvx; simgaStore(i,kkk,2,1) = dmy.sigma2.cvx; simgaStore(i,kkk,3,1) = dmy.sigma3.cvx(1);
		simgaStore(i,kkk,1,2) = dmy.sigma1.MSE^.5; simgaStore(i,kkk,2,2) = dmy.sigma2.MSE^.5; simgaStore(i,kkk,3,2) = dmy.sigma3.MSE^.5;
        simgaStore(i,kkk,1,3) = dmy.sigma1.analytic; simgaStore(i,kkk,2,3) = dmy.sigma2.analytic; simgaStore(i,kkk,3,3) = dmy.sigma3.analytic;


		%%% PLOTTING B1-B3
		if ~isempty(find(vals(i) == valPlt)) && ~isempty(find(frq(kkk) == frqPlt))
			cnt = cnt + 1;
			% rerun analysis for S
	        isDeltaS = 0; thresh = .01;
			%%%%%%%%% TEMPORAL FILTER ANALYZER %%%%%%%%%%%%%%%%%%%%%%%
			dmy = tfCompute(rawSig, thresh, 3, tFinForError,ODEparams, spikePeriod, tFin, isDeltaS);

			%% plot with S
			k = subplot(2,3,cnt);
			[pk, pkLoc] = findpeaks(rawSig.sig(3,:),rawSig.time);
			plot(pkLoc, pk, 'bo','MarkerSize',3); hold on;
			sigmaR = dmy.sigma1.cvx; sigmaU = dmy.sigma2.cvx; x = dmy.sigma3.cvx;
			plot(pkLoc, pk(end) + x(2)*exp(-pkLoc/sigmaR) + x(3)*exp(-pkLoc/sigmaU) + (pk(1)-x(2)-x(3)-pk(end))*exp(-pkLoc/x(1)), 'Color','blue')
			plot(pkLoc, pk(end) + x(2)*exp(-pkLoc/sigmaR) + x(3)*exp(-pkLoc/sigmaU), 'Color', [.7 .7 .7], 'LineStyle','--')
            k.FontSize = axFnt;
			legend('S_n^+','S^{d+f}_{fit}','S^{d+f}_{cut}'); ylim([0 1]); xlim([0 300]); xlabel('t [ms]','FontSize',lblFnt); text(10,.9,strcat('f_{spk} = ',num2str(frq(kkk))),'FontSize',axFnt)
            ylabel('S [Units]','FontSize',lblFnt);
			ttl = title(strcat('B',num2str(cnt-3)),'FontSize',ttlFnt);
            ttl.Units = 'Normalize'; 
            ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
            ttl.HorizontalAlignment = 'left';
		end

	end
end
plot(subplot(2,3,3), frq(strtPlotFrq:end), squeeze(simgaStore(1,strtPlotFrq:end,3,1)), '-bo','MarkerSize',3); hold on;
plot(subplot(2,3,3), frq(strtPlotFrq:end), squeeze(simgaStore(2,strtPlotFrq:end,3,1)), '-ro','MarkerSize',3);
plot(subplot(2,3,3), frq(strtPlotFrq:end), squeeze(simgaStore(3,strtPlotFrq:end,3,1)), '-go','MarkerSize',3);
plot(subplot(2,3,3), frq(strtPlotFrq:end), squeeze(simgaStore(1,strtPlotFrq:end,3,3)), '--b','MarkerSize',3); hold on;
plot(subplot(2,3,3), frq(strtPlotFrq:end), squeeze(simgaStore(2,strtPlotFrq:end,3,3)), '--r','MarkerSize',3);
plot(subplot(2,3,3), frq(strtPlotFrq:end), squeeze(simgaStore(3,strtPlotFrq:end,3,3)), '--g','MarkerSize',3);
dmy = [squeeze(simgaStore(2,3,3,1)) squeeze(simgaStore(2,9,3,1)) squeeze(simgaStore(2,15,3,1))];
plot(subplot(2,3,3), frqPlt(1), dmy(1), 'ko','MarkerSize',5,'LineWidth',2);
plot(subplot(2,3,3), frqPlt(2), dmy(2), 'ks','MarkerSize',5,'LineWidth',2);
plot(subplot(2,3,3), frqPlt(3), dmy(3), 'kd','MarkerSize',5,'LineWidth',2);
h(3).FontSize = axFnt;
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\sigma_{d+f,S} [ms]','FontSize',lblFnt); ylim([0 250]); xlim([0 frq(end)]); xlim([0 153]);
%legend('cvx \tau_{fac} = \tau_{dep}  = 100','cvx \tau_{fac} = \tau_{dep} = 200','cvx \tau_{fac} = \tau_{dep} = 300','cmptd \tau_{fac} = \tau_{dep}  = 100','cmptd \tau_{fac} = \tau_{dep} = 200','cmptd \tau_{fac} = \tau_{dep} = 300')
legend('\tau_{fac} = \tau_{dep}  = 100','\tau_{fac} = \tau_{dep} = 200','\tau_{fac} = \tau_{dep} = 300','FontSize',axFnt)
ttl = title('A3','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

'low pass error (should match abv)'
max(max(simgaStore(:,:,1,2)))

'high pass error (should match abv)'
max(max(simgaStore(:,:,2,2)))

'band pass error'
max(max(simgaStore(:,:,3,2)))

%% SigmaD Markers
plot(subplot(1,3,1), frq, squeeze(simgaStore(1,:,1,2)), '-bo','MarkerSize',3); hold on;
plot(subplot(1,3,1), frq, squeeze(simgaStore(2,:,1,2)), '-ro','MarkerSize',3);
plot(subplot(1,3,1), frq, squeeze(simgaStore(3,:,1,2)), '-go','MarkerSize',3);
xlabel('f_{spk} [Hz]'); ylabel('Error [S units]'); xlim([0 frq(end)]);
legend('\tau_{fac} = \tau_{dep}  = 100','\tau_{fac} = \tau_{dep} = 200','\tau_{fac} = \tau_{dep} = 300')
title(subplot(1,3,1),'A1','position',[30 .0006]);

%% SigmaF Markers
plot(subplot(1,3,2), frq, squeeze(simgaStore(1,:,2,2)), '-bo','MarkerSize',3); hold on;
plot(subplot(1,3,2), frq, squeeze(simgaStore(2,:,2,2)), '-ro','MarkerSize',3);
plot(subplot(1,3,2), frq, squeeze(simgaStore(3,:,2,2)), '-go','MarkerSize',3);
xlabel('f_{spk} [Hz]'); ylabel('Error [S units]'); xlim([0 frq(end)]);
legend('\tau_{fac} = \tau_{dep}  = 100','\tau_{fac} = \tau_{dep} = 200','\tau_{fac} = \tau_{dep} = 300')
title(subplot(1,3,2),'A2','position',[30 .00004]);

%% SigmaD+F Markers
plot(subplot(1,3,3), frq, squeeze(simgaStore(1,:,3,2)), '-bo','MarkerSize',3); hold on;
plot(subplot(1,3,3), frq, squeeze(simgaStore(2,:,3,2)), '-ro','MarkerSize',3);
plot(subplot(1,3,3), frq, squeeze(simgaStore(3,:,3,2)), '-go','MarkerSize',3);
xlabel('f_{spk} [Hz]'); ylabel('Error [S units]'); xlim([0 frq(end)]);
legend('\tau_{fac} = \tau_{dep}  = 100','\tau_{fac} = \tau_{dep} = 200','\tau_{fac} = \tau_{dep} = 300')
title(subplot(1,3,3),'A3','position',[25 .007]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%saveas(f1,'../img/fig_Squant2','epsc');



