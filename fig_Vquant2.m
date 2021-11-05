clear all; tic; close all;

axFnt = 12; ttlFnt = 23; lblFnt = 14;

ODEparams.G_ex = 1; ODEparams.tau_dec = 3; ODEparams.I_app = 0;
ODEparams.C = 1; ODEparams.U = .1; ODEparams.E_L = -60; ODEparams.E_ex = 0;
tFinForError = 500;

%dc = [.5 1.5]; dcIx = length(dc);
dc = [.5]; dcIx = length(dc);
frq = [10:10:150]; frqIx = length(frq);
fc = [200]; fcIx = length(fc);

for ii = 1:dcIx
    for jj = 1:frqIx
tic
        ODEparams.g_L = dc(ii);
        
        ODEparams.tau_fac = 0; ODEparams.tau_dep = fc;
        
        isDeltaS = 0; thresh = .01; tFin = tFinForError;
        [rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/frq(jj), tFin);
        dmy = tfCompute_classOnly(rawSig, thresh, 4, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
        
        store2b(ii,jj,:) = [NaN NaN NaN NaN NaN];
        
        if dmy.class == 3
            
            dmyS = tfCompute(rawSig, thresh, 3, tFinForError,ODEparams, 1000/frq(jj) , tFin, isDeltaS);
            store2a(ii,jj,:) = [dmyS.sigma1.cvx dmyS.sigma2.cvx dmyS.sigma3.cvx'];

            error = 9999;
            winners = [dmyS.sigma1.cvx 50 20 30 1];

            [pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            test = @(y, t1, t2, t3, A, B) pk(end) + A*exp(-(y-pkLoc(1))/t1) + B*exp(-(y-pkLoc(1))/t2) + (pk(1)-A-B-pk(end))*exp(-(y-pkLoc(1))/t3);
            testY = @(y) test(y, winners(1),winners(2),winners(3),winners(4),winners(5));
            error = sum([testY(pkLoc)-pk].^2)^.5;
            for B1 = [-15:1:-3] %[-15.1:1:-3.1] 
            for A1 = [3:1:15] %[3.1:3:15.1] 
            for sigma3 = [.1:.1:1] %[.1:.1:1]
            for sigma2 = [3:1:6] %[2:1:6] 
            for sigma1 = [24:2:40] %linspace(max(0,dmyS.sigma1.cvx(1)-30),dmyS.sigma1.cvx+30,20) %
                testY = @(y) test(y, sigma1, sigma2, sigma3,A1,B1);
                %figure; set(gcf,'Position', [10 10 1000 600]); plot(rawSig.time,rawSig.sig(4,:)); hold on; fplot(testY,[0 tFin]);        
                newerr = sum([testY(pkLoc)-pk].^2)^.5;
                if  newerr < error
                winners = [sigma1 sigma2 sigma3 A1 B1];
                error = newerr;
                end
            end
            end
            end
            end
            end
            %winners
			%error

            delta = @(x) pk - (pk(end) + x(4)*exp(-(pkLoc-pkLoc(1))/x(1)) + x(5)*exp(-(pkLoc-pkLoc(1))/x(2)) + (pk(1)-pk(end)-x(4)-x(5))*exp(-(pkLoc-pkLoc(1))/x(3)) );
            grad = @(x) [sum( (delta(x)*-2) .* (1/x(1)*x(4)*exp(-(pkLoc-pkLoc(1))/x(1))) ),...
                         sum( (delta(x)*-2) .* (1/x(2)*x(5)*exp(-(pkLoc-pkLoc(1))/x(2))) ),...
                         sum( (delta(x)*-2) .* (1/x(3)*(pk(1)-pk(end)-x(4)-x(5))*exp(-(pkLoc-pkLoc(1))/x(3))) ),...
                         sum( (delta(x)*2) .* (-exp(-(pkLoc-pkLoc(1))/x(1)) + exp(-(pkLoc-pkLoc(1))/x(3)))   ),...
                         sum( (delta(x)*2) .* (-exp(-(pkLoc-pkLoc(1))/x(2)) + exp(-(pkLoc-pkLoc(1))/x(3)))   )];
            % time step; intial conditions;
            alpha = [1, 1, 1, .25, .25]; x0 = winners'; %[dmyS.sigma1.cvx, dmyS.sigma2.cvx, dmyS.sigma3.cvx']';
            [xopt5,fopt5,niter,gnorm,dx] = gradDFit_V(delta, grad, alpha, x0);
            errorV2(ii,jj) = fopt5^.5;
            store2b(ii,jj,:) = xopt5;
            fitV = @(y) pk(end) + xopt5(4)*exp(-(y-pkLoc(1))/xopt5(1)) + xopt5(5)*exp(-(y-pkLoc(1))/xopt5(2)) + (pk(1)-pk(end)-xopt5(4)-xopt5(5))*exp(-(y-pkLoc(1))/xopt5(3));
			sigma3coeff(ii,jj) = pk(1)-pk(end)-xopt5(4)-xopt5(5);
            %
            %plotting fits
            %fitVStart = @(y) pk(end) + winners(4)*exp(-y/winners(1)) + winners(5)*exp(-y/winners(2)) + (pk(1)-pk(end)-winners(4)-winners(5))*exp(-y/winners(3));
			%fitV = fitVStart;
			%{
            winners = xopt5;
            fitV1 = @(y) pk(end) + winners(5)*exp(-(y-pkLoc(1))/winners(2)) + (pk(1)-pk(end)-winners(4)-winners(5))*exp(-(y-pkLoc(1))/winners(3));
            fitV2 = @(y) pk(end) + winners(4)*exp(-(y-pkLoc(1))/winners(1)) + (pk(1)-pk(end)-winners(4)-winners(5))*exp(-(y-pkLoc(1))/winners(3));
            fitV3 = @(y) pk(end) + winners(4)*exp(-(y-pkLoc(1))/winners(1)) + winners(5)*exp(-(y-pkLoc(1))/winners(2));

            figure; set(gcf,'Position', [10 10 1000 600]); plot(rawSig.time,rawSig.sig(4,:)); hold on; fplot(fitV,[0 tFin],'-k');
            fplot(fitV1,[0 tFin],'--k','linewidth',2); hold on; fplot(fitV2,[0 tFin],'--r','linewidth',2); fplot(fitV3,[0 tFin],'--b','linewidth',2);
            title(frq(jj)); legend('signal','start fit','exclude 1','exclude 2','exclude 3');
            %}
            

			[pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
            [diffs2(ii,jj) diffsIx2(ii,jj)] = max(abs(fitV(pkLoc)-pk));
            if diffs2(ii,jj) > 2
                failV2(ii,jj) = 1;
            else
                failV2(ii,jj) = 0;
            end
            %

        else
            store2b(ii,jj,1) = NaN;
            store2b(ii,jj,2) = NaN;
            store2b(ii,jj,3) = NaN;    
            store2b(ii,jj,4) = NaN;
			store2a(ii,jj,1:3) = [NaN NaN NaN];
            failV2(ii,jj) = NaN;
			sigma3coeff(ii,jj) = NaN;
            
        end  
  toc      
    end
end

f1 = figure(1); set(gcf,'Position', [10 10 1450 700])
choice = squeeze(store2b(1,1:end,2)); %(squeeze(store2b(1,1:end,2))+squeeze(store2b(1,1:end,3)))/2;
newCoeff = squeeze(store2b(1,:,5)); %sigma3coeff + squeeze(store2b(1,:,5));
ODEparams.g_L = dc(1);

h = subplot(2,3,4); j = 13;
[rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/130, tFin);
[pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
fitV3 = pk(end) + store2b(1,j,4)*exp(-(pkLoc-pkLoc(1))/store2b(1,j,1)) + newCoeff(1,j)*exp(-(pkLoc-pkLoc(1))/choice(j));
plot(rawSig.time,rawSig.sig(4,:),'-b'); hold on;
plot(pkLoc,fitV3,'--k','Linewidth',1.3);
xopt5 = store2b(1,j,:);
fitV = pk(end) + xopt5(4)*exp(-(pkLoc-pkLoc(1))/xopt5(1)) + xopt5(5)*exp(-(pkLoc-pkLoc(1))/xopt5(2)) + (pk(1)-pk(end)-xopt5(4)-xopt5(5))*exp(-(pkLoc-pkLoc(1))/xopt5(3));
plot(pkLoc,fitV,':k','Linewidth',1.3);
h.FontSize = axFnt;
legend({'V','2 Sigma','3 Sigma'},'FontSize',axFnt)
xlabel('t [ms]','FontSize',lblFnt); ylabel('V [mV]','FontSize',lblFnt); xlim([0 150]); %xlim([0 dc(end)]);
ttl = title('B1','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(2,3,5); j = 14;
[rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/140, tFin);
[pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
fitV3 = pk(end) + store2b(1,j,4)*exp(-(pkLoc-pkLoc(1))/store2b(1,j,1)) + newCoeff(1,j)*exp(-(pkLoc-pkLoc(1))/choice(j));
plot(rawSig.time,rawSig.sig(4,:),'-b'); hold on;
plot(pkLoc,fitV3,'--k','Linewidth',1.3);
xopt5 = store2b(1,j,:);
fitV = pk(end) + xopt5(4)*exp(-(pkLoc-pkLoc(1))/xopt5(1)) + xopt5(5)*exp(-(pkLoc-pkLoc(1))/xopt5(2)) + (pk(1)-pk(end)-xopt5(4)-xopt5(5))*exp(-(pkLoc-pkLoc(1))/xopt5(3));
plot(pkLoc,fitV,':k','Linewidth',1.3);
h.FontSize = axFnt;
legend({'V','2 Sigma','3 Sigma'},'FontSize',axFnt)
xlabel('t [ms]','FontSize',lblFnt); ylabel('V [mV]','FontSize',lblFnt); xlim([0 150]); %xlim([0 dc(end)]);
ttl = title('B2','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(2,3,6); j = 15;
[rawSig.sig, rawSig.time, isFired] = numerical(ODEparams, 1000/150, tFin);
[pk pkLoc] = findpeaks(rawSig.sig(4,:),rawSig.time);
winners = store2b(1,15,:);
fitV3 = pk(end) + winners(4)*exp(-(pkLoc-pkLoc(1))/winners(1)) + winners(5)*exp(-(pkLoc-pkLoc(1))/winners(2));
plot(rawSig.time,rawSig.sig(4,:),'-b'); hold on;
plot(pkLoc, fitV3,'--k','Linewidth',1.3);
xopt5 = store2b(1,j,:);
fitV = pk(end) + xopt5(4)*exp(-(pkLoc-pkLoc(1))/xopt5(1)) + xopt5(5)*exp(-(pkLoc-pkLoc(1))/xopt5(2)) + (pk(1)-pk(end)-xopt5(4)-xopt5(5))*exp(-(pkLoc-pkLoc(1))/xopt5(3));
plot(pkLoc,fitV,':k','Linewidth',1.3);h.FontSize = axFnt;
legend({'V','2 Sigma','3 Sigma'},'FontSize',axFnt)
xlabel('t [ms]','FontSize',lblFnt); ylabel('V [mV]','FontSize',lblFnt); xlim([0 150]); %xlim([0 dc(end)]);
ttl = title('B3','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(2,3,1);
plot(frq,squeeze(store2b(1,:,1)),'-sb'); hold on;
plot(frq,squeeze(store2a(1,:,1)),'-sk'); hold on;
h.FontSize = axFnt;
legend({'\tau_{mem} = 2','S'},'FontSize',axFnt,'Location','northwest')
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\rho_{1} [ms]','FontSize',lblFnt); xlim([0 160]); ylim([0 50]);
ttl = title('A1','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

h = subplot(2,3,2);
a = plot(frq,squeeze(store2b(1,:,2)),'-sb'); hold on;
b = plot(frq,squeeze(store2b(1,:,3)),':sb');
h.FontSize = axFnt;
legend([a b], {'\rho_{2}','\rho_{3}'},'FontSize',axFnt)
xlabel('f_{spk} [Hz]','FontSize',lblFnt); ylabel('\rho [ms]','FontSize',lblFnt); xlim([0 160]); ylim([0 10]);
ttl = title('A2','FontSize',ttlFnt);
ttl.Units = 'Normalize'; 
ttl.Position(1) = 0; % use negative values (ie, -0.1) to move further left
ttl.HorizontalAlignment = 'left';

saveas(f1,'../img/fig_Vquant2','epsc');


