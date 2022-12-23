function [DT, counterTime, IterationCounter, stableCounter, ...
          maxStableDT] = updateDT(DT, minDT, maxDT, ...
                                  ITER, IterationCounter, ...
                                  counterTime, maxCounterTime, ...
                                  stableCounter, maxStableDT, ...
                                  nstep)
%% Update time step increment
if nstep > 2*maxCounterTime
    IterationCounter(counterTime) = ITER(2);
    counterTime = counterTime + 1;
    if ITER(2) > 5
        if DT == maxStableDT
            maxStableDT = 0;
        end
        DT = DT * 1/4;
        if DT < minDT
            DT = minDT;
        end
        stableCounter = 1;
        counterTime = 1;
        IterationCounter = zeros(maxCounterTime,1);
    end
    if counterTime > maxCounterTime
        if ~any(IterationCounter > 5)
            if maxStableDT > DT 
                DT = maxStableDT;
                stableCounter = 1;
            elseif DT > maxStableDT
                if stableCounter > 5
                    maxStableDT = DT;
                    DT = DT*1.1;
                else
                    maxStableDT = DT;
                    DT = DT*1.5;
                    stableCounter = 1;
                end
            elseif DT == maxStableDT
                stableCounter = stableCounter + 1;
                if stableCounter > 5
                    DT = DT*1.1;
                end
            end
            if DT > maxDT
                DT = maxDT;
            end
        end
        counterTime = 1;
        IterationCounter = zeros(maxCounterTime,1);
    end
end