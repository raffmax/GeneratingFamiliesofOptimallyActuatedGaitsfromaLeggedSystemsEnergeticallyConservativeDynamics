function trajectory = getCompleteTrajectory(obj,trajectory)
    nTime = length(trajectory); %obj.systemRef.nTime;
    flipOperator = obj.systemRef.flipOperator;
    T = trajectory(end).t(end);
    % get periodicStates
    if isa(obj.systemRef,'LimitCycle')
        periodicStates = obj.systemRef.Sol.rfpData.periodicStates;
    else
        periodicStates = obj.systemRef.periodicStates;
    end
    stateIdx  = (1:obj.systemRef.nState);
    offsetIdx = stateIdx(~ismember(stateIdx,periodicStates));
    for iDomain = 1:nTime
        trajectory(iDomain+nTime).t = trajectory(iDomain).t-trajectory(iDomain).t(1)+trajectory(nTime+iDomain-1).t(end);
        trajectory(iDomain+nTime).x = trajectory(iDomain).x*flipOperator';
        % add offset
        trajectory(iDomain+nTime).x(:,offsetIdx) = trajectory(iDomain).x(:,offsetIdx)-trajectory(iDomain).x(1,offsetIdx)+ trajectory(iDomain+nTime-1).x(end,offsetIdx);
        if size(trajectory(iDomain).z,2) == 4
            trajectory(iDomain+nTime).z = trajectory(iDomain).z(:,[2,1,4,3]);
        elseif size(trajectory(iDomain).z,2) == 5
            trajectory(iDomain+nTime).z = trajectory(iDomain).z(:,[2,1,4,3,5]);
        else
            trajectory(iDomain+nTime).z = flip(trajectory(iDomain).z,2);
        end
        % get f, lambda, tau
        t = trajectory(iDomain).t;
        x = trajectory(iDomain+nTime).x;
        z = trajectory(iDomain+nTime).z;
        for i = 1:length(t)
            [f, lambda,~,~,tau] = obj.systemRef.model.FlowMap(t(i),x(i,:)',z(i,:)',T, obj.systemRef.controller, false, false);
            trajectory(iDomain+nTime).tau(i,:) = -tau';
            trajectory(iDomain+nTime).f(i,:) = f';
            lCount = 0;
            for iZ = 1:obj.systemRef.model.nZ
                if z(i,iZ) 
                    nLambda = obj.systemRef.model.nLambda(iZ);
                    trajectory(iDomain+nTime).lambda(i,1:nLambda,iZ) = lambda(lCount+(1:nLambda))';
                    lCount = lCount + nLambda;
                end
            end
        end
    end
end