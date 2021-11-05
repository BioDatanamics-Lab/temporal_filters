function [tfStuff] = tfCompute(rawSig, thresh, sigIdx, tFinForError, ODEparams, spikePeriod, tFin, isDeltaS)

    % intialize TF selector
	if isDeltaS == 1
        'this section is depricated'
        asdf
		%frq = 1000/spikePeriod;
		%[soln_r] = recurrent(ODEparams,frq);
		%pk = soln_r(3).peaks - soln_r(3).traughs;
		%pkLoc = soln_r(3).peakTimes;
	else
		[pk, pkLoc] = findpeaks(rawSig.sig(sigIdx,:),rawSig.time);
	end
	firstPk = pk(1); [maxPk, midx] = max(pk); ssPk = pk(end);
	xq = [pkLoc(1):.01:pkLoc(end)]; sS = spline(pkLoc,pk,xq);

    if abs(firstPk - maxPk) < thresh && abs(ssPk - maxPk) < thresh   %even across the board
        tfStuff.class = 1;
    elseif (maxPk - firstPk) >= thresh && abs(ssPk - maxPk) < (maxPk - firstPk)*.05 %|| check.fac == 1 %increases to plateu
        tfStuff.class = 2;
        [tfStuff.sigma1.cvx, tfStuff.sigma1.MSE, tfStuff.sigma1.analytic] = sigmaSFacOrDep(rawSig,sS,0,sigIdx, tFinForError,spikePeriod,ODEparams,isDeltaS);
        tfStuff.sigma2.cvx = NaN; tfStuff.sigma3.cvx = NaN;
        tfStuff.sigma2.MSE = NaN; tfStuff.sigma3.MSE = NaN;
        tfStuff.sigma2.analytic = NaN; tfStuff.sigma3.analytic = NaN;
    elseif abs(firstPk - maxPk) < thresh && (maxPk - ssPk) >= thresh %|| check.dep == 1%decay
        tfStuff.class = 4;
        [tfStuff.sigma1.cvx, tfStuff.sigma1.MSE, tfStuff.sigma1.analytic] = sigmaSFacOrDep(rawSig,sS,1,sigIdx,tFinForError,spikePeriod,ODEparams,isDeltaS);
        tfStuff.sigma2.cvx = NaN; tfStuff.sigma3.cvx = NaN;
        tfStuff.sigma2.MSE = NaN; tfStuff.sigma3.MSE = NaN;
        tfStuff.sigma2.analytic = NaN; tfStuff.sigma3.analytic = NaN;
    elseif (maxPk - firstPk) >= thresh && (maxPk - ssPk) >= (maxPk - firstPk)*.05 %temporary peak
        tfStuff.class = 3;

		facStore = ODEparams.tau_fac; depStore = ODEparams.tau_dep;

		% high pass TF
		ODEparams.tau_dep = 0;
		[rawSig.sig, rawSig.time, rawSig.isFired] = numerical(ODEparams, spikePeriod, tFin);
        %chk.dep = 1; chk.fac = 0;
		tfStuff_low = tfCompute(rawSig, .01, sigIdx, tFinForError, ODEparams, spikePeriod, tFin, isDeltaS);
        tfStuff.sigma2 = tfStuff_low.sigma1;

		% low pass TF
		ODEparams.tau_fac = 0; ODEparams.tau_dep = depStore;
		[rawSig.sig, rawSig.time, rawSig.isFired] = numerical(ODEparams, spikePeriod, tFin);
        %chk.dep = 0; chk.fac = 1;
		tfStuff_high = tfCompute(rawSig, .01, sigIdx, tFinForError, ODEparams, spikePeriod, tFin, isDeltaS);
		tfStuff.sigma1 = tfStuff_high.sigma1;

		sigmaR = tfStuff.sigma1.cvx; sigmaU = tfStuff.sigma2.cvx;
		analyticR = tfStuff.sigma1.analytic; analyticU = tfStuff.sigma2.analytic;

        
        % CVX fit
        % intialize two fixed/one free sigma fit
        error = 999; winners = [9999 9999 9999];
        for A1 = (-3:.3:3)
        for B1 = (-3:.3:3)
        for sigma1 = [1:5:100] %[1:2:10,20:25:100,200,300,400]
                test = @(y, A, B, sigma) pk(end) + A*exp(-y/sigmaR) + B*exp(-y/sigmaU) + (pk(1)-A-B-pk(end))*exp(-y/sigma);
                testY = @(y) test(y, A1, B1, sigma1);
                newerr = sum([testY(pkLoc)-pk].^2);
            if  newerr < error
                winners = [A1 B1 sigma1];
                error = newerr;
            end
            %title(strcat(num2str(A),',',num2str(B)));
            %pause
        end
        end
        end
        
        testY = @(y) test(y, (1/sigmaR+1/sigmaU)^-1, -.3, .1);
        check = sum([testY(pkLoc)-pk].^2);
        if check < error
            winners = [(1/sigmaR+1/sigmaU)^-1, -.3 .1];
        end
    
        % sigma - x(1), A - x(2), B - x(3)
        delta = @(x) pk - (pk(end) + x(2)*exp(-pkLoc/sigmaR) + x(3)*exp(-pkLoc/sigmaU) + (pk(1)-x(2)-x(3)-pk(end))*exp(-pkLoc/x(1)) );
        grad = @(x) [sum( (delta(x)*-2) .* (1/x(1)*(pk(1)-x(2)-x(3)-pk(end))*exp(-pkLoc/x(1))) ),...
                     sum( (delta(x)*2) .* (-exp(-pkLoc/sigmaR) + exp(-pkLoc/x(1)))   ),...
                     sum( (delta(x)*2) .* (-exp(-pkLoc/sigmaU) + exp(-pkLoc/x(1)))   )];
        % time step; intial conditions;
        alpha = [.25, .075, .075]; x0 = [winners(3) winners(1) winners(2)]'; %x0 = [(1/sigmaR+1/sigmaU)^-1, -.3, .1]';
        [xopt3,fopt3,niter,gnorm,dx] = gradDFit(delta, grad, alpha, x0);
		tfStuff.sigma3.cvx = xopt3;

		% error computation.
		fitS1 = @(x) pk(end) + xopt3(2)*exp(-x/sigmaR) + xopt3(3)*exp(-x/sigmaU) + (pk(1)-xopt3(2)-xopt3(3)-pk(end))*exp(-x/xopt3(1));
		errIdx = find(pkLoc<tFinForError);
		tfStuff.sigma3.MSE = sum([fitS1(pkLoc(errIdx))-pk(errIdx)].^2)/length(errIdx);
        
        tfStuff.sigma3.analytic = (analyticR^-1 + analyticU^-1)^-1;

		% two sigma fit.
        % [tfStuff.twoSigmaFit.initConditionsWarning, tfStuff.twoSigmaFit.fits, tfStuff.twoSigmaFit.MSE] = twoSigmaFit(rawSig,xq,3,tFinForError);
    else
        'TF doesnt appear to fit any classification' 
        keyboard
    end

end

function [sigmas,MSEs,sigma0] = sigmaSFacOrDep(rawSig,sS,isDecay,sigIdx,tFinForError,spikePeriod,ODEparams,isDeltaS)

	if isDeltaS == 1
		%frq = 1000/spikePeriod;
		%[soln_r] = recurrent(ODEparams,frq);
		%pk = soln_r(3).peaks - soln_r(3).traughs;
		%pkLoc = soln_r(3).peakTimes;
        'this section is depricated'
        asdf
	else
		[pk, pkLoc] = findpeaks(rawSig.sig(sigIdx,:),rawSig.time);
	end
    firstPk = pk(1); [maxPk, midx] = max(pk); ssPk = pk(end);

    if isDecay == 1
        xinf = ssPk; 
    elseif isDecay == 0
        xinf = maxPk;
    else
        'huh?'
        keyboard
    end

    % manual fit
    line63 = xinf -((xinf-firstPk)*(exp(-1)));
    if isDecay == 1
        idxFind = find([sS < line63]);
    elseif isDecay == 0
        idxFind = find([sS > line63]);
    end
    sigma0 = (idxFind(1)/10^2-pkLoc(1)/10^2)/1;

    %err = [pk - (pk(end)-(pk(end) - pk(1))*exp(-pkLoc/sigma0))];
    %MSEs0 = err*err';

    %% CVX Fit
    delta = @(x) pk - (pk(end)-(pk(end) - pk(1))*exp(-(pkLoc-pkLoc(1))/x(1)));
    grad = @(x) sum( (delta(x)*2) .* (1/x(1)*(pk(end)-pk(1))*exp(-(pkLoc-pkLoc(1))/x(1))) );
    % time step; intial conditions;
    alpha = 1; x0 = sigma0;
    [xopt1,fopt1,niter,gnorm,dx] = gradDFit(delta, grad, alpha, x0);
    fitS1 = @(x) pk(end)-(pk(end) - pk(1))*exp(-(x-pkLoc(1))/xopt1);
	sigmas = xopt1;

	%% error computation.
	errIdx = find(pkLoc<tFinForError);
	MSEs = sum([fitS1(pkLoc(errIdx))-pk(errIdx)].^2)/length(errIdx);

end