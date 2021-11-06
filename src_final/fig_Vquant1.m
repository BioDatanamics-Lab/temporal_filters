clear all; tic; close all;
% note: this takes about an hour to run
tic
axFnt = 12; ttlFnt = 23; lblFnt = 14;

ODEparams.G_ex = 1; ODEparams.tau_dec = 3; ODEparams.I_app = 0;
ODEparams.C = 1; ODEparams.U = .1; ODEparams.E_L = -60; ODEparams.E_ex = 0;
tFinForError = 500;

%dc = [.1 .2 .3 .5 .6 .7 .8 1 1.5 2];
dc = [.5 1.5];
dcIx = length(dc);
frq = [10:10:150]; frqIx = length(frq);
fc = [200]; fcIx = length(fc);

% Compute and Store High Pass Fitlers of V
for ii = 1:dcIx
    for jj = 1:frqIx
        for jjj = 1:fcIx

        ODEparams.g_L = dc(ii);
        
        ODEparams.tau_dep = 0; ODEparams.tau_fac = fc(jjj);
        
        isDeltaS = 0; thresh = .01; tFin = tFinForError;
        [rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/frq(jj), tFin);
        dmy = tfCompute_classOnly(rawSig, thresh, 4, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
 
        if dmy.class ~= 3
            dmy = tfCompute(rawSig, thresh, 4, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
            dmyS = tfCompute(rawSig, thresh, 3, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
            store3a(ii,jj,1) = dmy.sigma1.cvx;
            store3a(ii,jj,2) = dmy.sigma1.MSE^.5;
            store3a(ii,jj,3) = dmy.sigma1.analytic;
            store3a(ii,jj,4) = dmyS.sigma1.cvx;
            
            [pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            fitV = pk(end) - (pk(end) - pk(1)).*exp(-(pkLoc-pkLoc(1))/store3a(ii,jj,1));
            [diffs3a(ii,jj) diffsIx3a(ii,jj)] = max(abs(fitV-pk));
            if diffs3a(ii,jj) > 2
                fail3a(ii,jj) = 1;
            else
                fail3a(ii,jj) = 0;
            end
            
            %{
            figure;
            plot(rawSig.time,rawSig.sig(4,:)); hold on;
            [pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            fitV = pk(end) - (pk(end) - pk(1)).*exp(-(pkLoc-pkLoc(1))/store3a(ii,jj,1));
            %diffs = fitV-pk; find(diffs > 2);
            plot(pkLoc, fitV,'-k');
            title(frq(jj))
            %keyboard
            %}
        else
            store3a(ii,jj,1) = NaN;
            store3a(ii,jj,2) = NaN;
            store3a(ii,jj,3) = NaN;    
            store3a(ii,jj,4) = NaN;
            fail3a(ii,jj) = NaN;
            
        end  
        
    end
    end
end


% Compute and Store Low Pass Fitlers of V
for ii = 1:dcIx
    for jj = 1:frqIx
        for jjj = 1:fcIx

        ODEparams.g_L = dc(ii);
        
        ODEparams.tau_fac = 0; ODEparams.tau_dep = fc(jjj);
        
        isDeltaS = 0; thresh = .01; tFin = tFinForError;
        [rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/frq(jj), tFin);
        dmy = tfCompute_classOnly(rawSig, thresh, 4, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
 
        if dmy.class ~= 3
            dmy = tfCompute(rawSig, thresh, 4, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
            dmyS = tfCompute(rawSig, thresh, 3, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
            store3b(ii,jj,1) = dmy.sigma1.cvx;
            store3b(ii,jj,2) = dmy.sigma1.MSE^.5;
            store3b(ii,jj,3) = dmy.sigma1.analytic;
            store3b(ii,jj,4) = dmyS.sigma1.cvx;
            
            [pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            fitV = pk(end) - (pk(end) - pk(1)).*exp(-(pkLoc-pkLoc(1))/store3b(ii,jj,1));
            [diffs3b(ii,jj) diffsIx3b(ii,jj)] = max(abs(fitV-pk));
            if diffs3b(ii,jj) > 2
                fail3b(ii,jj) = 1;
            else
                fail3b(ii,jj) = 0;
            end
            
            %{
            figure;
            plot(rawSig.time,rawSig.sig(4,:)); hold on;
            [pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            fitV = pk(end) - (pk(end) - pk(1)).*exp(-(pkLoc-pkLoc(1))/store3b(ii,jj,1));
            plot(pkLoc, fitV, '-k');
            title(frq(jj))
            %keyboard
            %}
        else
            store3b(ii,jj,1) = NaN;
            store3b(ii,jj,2) = NaN;
            store3b(ii,jj,3) = NaN;    
            store3b(ii,jj,4) = NaN;
            fail3b(ii,jj) = NaN;
            
        end  
        
    end
    end
end

for ii = 1:dcIx
    for jj = 1:frqIx
        
        ODEparams.g_L = dc(ii);
        
        ODEparams.tau_dep = fc; ODEparams.tau_fac = fc;
        
        isDeltaS = 0; thresh = .01; tFin = tFinForError;
        [rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/frq(jj), tFin);
        dmy = tfCompute_classOnly(rawSig, thresh, 4, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
 
        if dmy.class == 3
      
            dmyS = tfCompute(rawSig, thresh, 3, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
            store(ii,jj,1,:) = [dmyS.sigma1.cvx dmyS.sigma2.cvx dmyS.sigma3.cvx'];
                       
            error = 9999;
            winners = [dmyS.sigma1.cvx dmyS.sigma2.cvx dmyS.sigma3.cvx'];

            [pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            test = @(y, t1, t2, t3, A, B) pk(end) + A*exp(-(y-pkLoc(1))/t1) + B*exp(-(y-pkLoc(1))/t2) + (pk(1)-A-B-pk(end))*exp(-(y-pkLoc(1))/t3);
            testY = @(y) test(y, winners(1),winners(2),winners(3),winners(4),winners(5));
            error = sum([testY(pkLoc)-pk].^2);
            for B1 = [-30:3:30] %30 
            for A1 = [-30:3:30] %1 %(-3:.3:3)
            for sigma3 = linspace(max(0,dmyS.sigma3.cvx(1)-30),dmyS.sigma3.cvx(1)+30,20) %winners(3) %
            for sigma2 = linspace(max(0,dmyS.sigma2.cvx(1)-30),dmyS.sigma2.cvx+30,20) %winners(2)
            for sigma1 = linspace(max(0,dmyS.sigma1.cvx(1)-30),dmyS.sigma1.cvx+30,20) %winners(1)
                testY = @(y) test(y, sigma1, sigma2, sigma3,A1,B1);
                %figure; set(gcf,'Position', [10 10 1000 600]); plot(rawSig.time,rawSig.sig(4,:)); hold on; fplot(testY,[0 tFin]);        
                newerr = sum([testY(pkLoc)-pk].^2);
                if  newerr < error
                winners = [sigma1 sigma2 sigma3 A1 B1];
                error = newerr;
                end
            %title(strcat(num2str(A),',',num2str(B)));
            %pause
            end
            end
            end
            end
            end
            
            %
            delta = @(x) pk - (pk(end) + x(4)*exp(-(pkLoc-pkLoc(1))/x(1)) + x(5)*exp(-(pkLoc-pkLoc(1))/x(2)) + (pk(1)-pk(end)-x(4)-x(5))*exp(-(pkLoc-pkLoc(1))/x(3)) );
            grad = @(x) [sum( (delta(x)*-2) .* (1/x(1)*x(4)*exp(-(pkLoc-pkLoc(1))/x(1))) ),...
                         sum( (delta(x)*-2) .* (1/x(2)*x(5)*exp(-(pkLoc-pkLoc(1))/x(2))) ),...
                         sum( (delta(x)*-2) .* (1/x(3)*(pk(1)-pk(end)-x(4)-x(5))*exp(-(pkLoc-pkLoc(1))/x(3))) ),...
                         sum( (delta(x)*2) .* (-exp(-(pkLoc-pkLoc(1))/x(1)) + exp(-(pkLoc-pkLoc(1))/x(3)))   ),...
                         sum( (delta(x)*2) .* (-exp(-(pkLoc-pkLoc(1))/x(2)) + exp(-(pkLoc-pkLoc(1))/x(3)))   )];
            % time step; intial conditions;
            alpha = [2, 2, 2, .5, .5]; x0 = winners'; %[dmyS.sigma1.cvx, dmyS.sigma2.cvx, dmyS.sigma3.cvx']';
            [xopt5,fopt5,niter,gnorm,dx] = gradDFit_V(delta, grad, alpha, x0);
            store(ii,jj,2,:) = xopt5;
            fitV = @(y) pk(end) + xopt5(4)*exp(-(y-pkLoc(1))/xopt5(1)) + xopt5(5)*exp(-(y-pkLoc(1))/xopt5(2)) + (pk(1)-pk(end)-xopt5(4)-xopt5(5))*exp(-(y-pkLoc(1))/xopt5(3)); 
            
            % plot testing.
            %{
            fitVStart = @(y) pk(end) + winners(4)*exp(-(y-pkLoc(1))/winners(1)) + winners(5)*exp(-(y-pkLoc(1))/winners(2)) + (pk(1)-pk(end)-winners(4)-winners(5))*exp(-(y-pkLoc(1))/winners(3));
            figure; set(gcf,'Position', [10 10 1000 600]); plot(rawSig.time,rawSig.sig(4,:)); hold on; fplot(fitVStart,[0 tFin]); title(frq(jj));
            %
            fplot(fitV,[0 tFin]); legend('signal','end fit','start fit');
            %}
            [pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            [diffs(ii,jj) diffsIx(ii,jj)] = max(abs(fitV(pkLoc)-pk));
            
            if diffs(ii,jj) > 1
                failV1(ii,jj) = 1;
            else
                failV1(ii,jj) = 0;
            end
            
            errIdx = find(pkLoc<tFinForError);
            errorV1(ii,jj) = (sum([fitV(pkLoc(errIdx))-pk(errIdx)].^2)/length(errIdx))^.5;
                
            %{
            figure;
            plot(rawSig.time,rawSig.sig(4,:)); hold on;
            [pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            fitV = pk(end) - (pk(end) - pk(1)).*exp(-(pkLoc-pkLoc(1))/store(ii,jj,1));
            diffs = fitV-pk; find(diffs > 2);
            plot(pkLoc, fitV,'-k');
            title(frq(jj))
            %keyboard
            %}
            
        else
            store(ii,jj,1) = NaN;
            store(ii,jj,2) = NaN;
            store(ii,jj,3) = NaN;    
            store(ii,jj,4) = NaN;
            failV1(ii,jj) = NaN;
            
        end  
        
    end
end
%save('Vquant1_small.mat')
%clear all; close all; clc
%load('Vquant1_big.mat')

dcPlotIx = [1 2];
lgdStrCat = {strcat('\tau_{mem} = ',num2str((1/dc(dcPlotIx(1))))),strcat('\tau_{mem} = ',num2str(round(1/dc(dcPlotIx(2)),2))),'S'};
f1 = figure(1); set(gcf,'Position', [10 10 1450 950])
h = subplot(3,3,1);
% frq = 80; gl = .5; dep fac are  200 200 when appropropriate.
ODEparams.g_L = .5; ODEparams.tau_dep = 0; ODEparams.tau_fac = 200;
[rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/80, tFin);
plot(rawSig.time,rawSig.sig(4,:),'-b'); hold on;
[pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
fitV = pk(end) - (pk(end) - pk(1)).*exp(-(pkLoc-pkLoc(1))/store3a(dcPlotIx(1),8,1));
plot(pkLoc,fitV,'-k');
h.FontSize = axFnt;
legend({'V','Temporal Filter'},'FontSize',axFnt,'Location','northwest')
xlabel('t [ms]','FontSize',lblFnt); ylabel('V [mV]','FontSize',lblFnt); xlim([0 150]); ylim([-60 -30]);
ttl = title('A1','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(3,3,2);
ODEparams.g_L = .5; ODEparams.tau_dep = 200; ODEparams.tau_fac = 0;
[rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/80, tFin);
plot(rawSig.time,rawSig.sig(4,:),'-b'); hold on;
[pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
fitV = pk(end) - (pk(end) - pk(1)).*exp(-(pkLoc-pkLoc(1))/store3b(dcPlotIx(1),8,1));
plot(pkLoc,fitV,'-k');
h.FontSize = axFnt;
legend({'V','Temporal Filter'},'FontSize',axFnt)
xlabel('t [ms]','FontSize',lblFnt); ylabel('V [mV]','FontSize',lblFnt); xlim([0 150]); ylim([-60 -45]);
ttl = title('A2','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(3,3,3);
ODEparams.g_L = .5; ODEparams.tau_dep = 200; ODEparams.tau_fac = 200;
[rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/80, tFin);
plot(rawSig.time,rawSig.sig(4,:),'-b'); hold on;
[pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
xopt5 = store(dcPlotIx(1),8,2,:);
fitV = pk(end) + xopt5(4)*exp(-(pkLoc-pkLoc(1))/xopt5(1)) + xopt5(5)*exp(-(pkLoc-pkLoc(1))/xopt5(2)) + (pk(1)-pk(end)-xopt5(4)-xopt5(5))*exp(-(pkLoc-pkLoc(1))/xopt5(3));
plot(pkLoc,fitV,'-k');
fitVcut = pk(end) + xopt5(4)*exp(-(pkLoc-pkLoc(1))/xopt5(1)) + xopt5(5)*exp(-(pkLoc-pkLoc(1))/xopt5(2));
plot(pkLoc,fitVcut, 'Color', [.7 .7 .7], 'LineStyle','--')
fitVcut = pk(end) + xopt5(4)*exp(-(pkLoc-pkLoc(1))/xopt5(1)) + (pk(1)-pk(end)-xopt5(4)-xopt5(5))*exp(-(pkLoc-pkLoc(1))/xopt5(3));
plot(pkLoc,fitVcut, 'Color', [.7 .7 .7], 'LineStyle','-.')
fitVcut = pk(end) + xopt5(5)*exp(-(pkLoc-pkLoc(1))/xopt5(2)) + (pk(1)-pk(end)-xopt5(4)-xopt5(5))*exp(-(pkLoc-pkLoc(1))/xopt5(3));
plot(pkLoc,fitVcut, 'Color', [.7 .7 .7], 'LineStyle',':')
h.FontSize = axFnt;
legend({'V','Temporal Filter'},'FontSize',axFnt)
xlabel('t [ms]','FontSize',lblFnt); ylabel('V [mV]','FontSize',lblFnt); xlim([0 150]); ylim([-60 -45]);
ttl = title('A3','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(3,3,4); strt = 2;
plot(frq(strt:end),squeeze(store3a(dcPlotIx(1),strt:end,1)),'-sb'); hold on;
plot(frq(strt:end),squeeze(store3a(dcPlotIx(2),strt:end,1)),'-sr'); hold on;
plot(frq(strt:end),squeeze(store3a(dcPlotIx(2),strt:end,4)),'-sk'); hold on;
h.FontSize = axFnt;
legend(lgdStrCat,'FontSize',axFnt)
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\sigma_{f,V} [ms]','FontSize',lblFnt); ylim([0 160]);  xlim([10 160]);
ttl = title('B1','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(3,3,5);
plot(frq(strt:end),squeeze(store3b(dcPlotIx(1),strt:end,1)),'-sb'); hold on;
plot(frq(strt:end),squeeze(store3b(dcPlotIx(2),strt:end,1)),'-sr'); hold on;
plot(frq(strt:end),squeeze(store3b(dcPlotIx(2),strt:end,4)),'-sk'); hold on;
plot(frq(12),squeeze(store3b(dcPlotIx(1),12,1)),'-xb','MarkerSize',10,'LineWidth',2); hold on;
h.FontSize = axFnt;
legend(lgdStrCat,'FontSize',axFnt)
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\sigma_{d,V} [ms]','FontSize',lblFnt); ylim([0 160]); xlim([10 160]);;
ttl = title('B2','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

 
for i = 1:2
    for j = 1:15
        gLsort(i,j,:) = sort(squeeze(store(i,j,2,1:3)))';
    end
end
 
h = subplot(3,3,8); strt = 2;
plot(frq(strt:end),squeeze(gLsort(dcPlotIx(1),strt:end,2)),'-sb'); hold on;
plot(frq(strt:end),squeeze(gLsort(dcPlotIx(2),strt:end,2)),'-sr'); hold on;
plot(frq(strt:end),squeeze(store(1,strt:end,1,1)),'-sk'); hold on;
h.FontSize = axFnt;
legend(lgdStrCat,'FontSize',axFnt)
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\rho_{1} [ms]','FontSize',lblFnt);  ylim([0 160]); xlim([10 160]);
ttl = title('C2','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(3,3,7);
plot(frq(strt:end),squeeze(gLsort(dcPlotIx(1),strt:end,3)),'-sb'); hold on;
plot(frq(strt:end),squeeze(gLsort(dcPlotIx(2),strt:end,3)),'-sr'); hold on;
plot(frq(strt:end),squeeze(store(1,strt:end,1,2)),'-sk'); hold on;
h.FontSize = axFnt;
legend(lgdStrCat,'FontSize',axFnt)
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\rho_{2} [ms]','FontSize',lblFnt); ylim([0 160]); xlim([10 160]);
ttl = title('C1','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(3,3,9);
plot(frq(strt:end),squeeze(gLsort(dcPlotIx(1),strt:end,1)),'-sb'); hold on;
plot(frq(strt:end),squeeze(gLsort(dcPlotIx(2),strt:end,1)),'-sr'); hold on;
plot(frq(strt:end),squeeze(store(1,strt:end,1,3)),'-sk'); hold on;
h.FontSize = axFnt;
legend(lgdStrCat,'FontSize',axFnt)
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\rho_{3} [ms]','FontSize',lblFnt); ylim([0 160]); xlim([10 160]);
ttl = title('C3','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

%load('Vquant1_small.mat','errorV1')
dmy1 = store3b(dcPlotIx,:,2); dmy2 = store3a(dcPlotIx,:,2); dmy3 = errorV1(:); dmy = [dmy1(:)' dmy2(:)' dmy3(:)'];
'MSE error Upper bound'
max(dmy)

toc

%saveas(f1,'../img/fig_Vquant1','epsc');
%}