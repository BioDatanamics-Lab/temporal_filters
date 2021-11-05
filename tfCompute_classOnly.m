function [tfStuff] = tfCompute_classOnly(rawSig, thresh, sigIdx, tFinForError, ODEparams, spikePeriod, tFin, isDeltaS)

    [pk, pkLoc] = findpeaks(rawSig.sig(sigIdx,:),rawSig.time);
    firstPk = pk(1); [maxPk, midx] = max(pk); ssPk = pk(end);
    xq = [pkLoc(1):.01:pkLoc(end)]; sS = spline(pkLoc,pk,xq);

    if abs(firstPk - maxPk) < thresh && abs(ssPk - maxPk) < thresh   %even across the board
        tfStuff.class = 1;
    elseif (maxPk - firstPk) >= thresh && abs(ssPk - maxPk) < (maxPk - firstPk)*.05 %|| check.fac == 1 %increases to plateu
        tfStuff.class = 2;
    elseif abs(firstPk - maxPk) < thresh && (maxPk - ssPk) >= thresh %|| check.dep == 1%decay
        tfStuff.class = 4;
    elseif (maxPk - firstPk) >= thresh && (maxPk - ssPk) >= (maxPk - firstPk)*.05 %temporary peak
        tfStuff.class = 3;
    else
        'da foq? courseclass'
        keyboard
    end

end