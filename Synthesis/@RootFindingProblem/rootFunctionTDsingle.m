function [fun,jacobian,objective] = rootFunctionTDsingle(obj,init,x,computeLxx)
%rootFunctionTDsingle(obj,init,decVar)
%   This function uses time-driven single shooting in simulation to obtain
%   the rootFunction value

% computeLxx == 0 % don't compute second derivative in optimization problem
%                 % use scaling
% computeLxx == 1 % compute second derivative in optimization problem
%                 % use scaling
% computeLxx == 2 % don't compute second derivative in optimization problem
%                 % don't use scaling

if nargin<4
    % don't compute second derivative in optimization problem
    computeLxx = 0;
end

% get initial state, inital domain durations and gradient flags 
[x0,tDomain,~,epsilon,decVar] = obj.getStateTimeCtrlEpsFromDecVar(x,init);
z0              = obj.sequence.FlowMapOrder(:,1);
zT              = obj.sequence.JumpMapOrder(:,end);
T               = sum(tDomain);
Grad1           = obj.options.Grad1;
Grad2           = obj.options.Grad2;

%% run simulation
[xT,simData,~]    = obj.simulation.singleShootingTD(x0,tDomain,Grad1,Grad2);

% define structure of jacobian
if Grad1
    % allocate space for Jacobian
    jacobian = zeros(obj.nCon,obj.nDecVar);
    % compose gradient information
    [xE,intE] = obj.getGradientsFromSimData(simData);
    % compute discrete map to xT
    xT_xE = obj.flipOperator*simData.jumpMapData.xGrad1(:,:,end);
else
    jacobian = [];
end

%% define periodicity constraints
% compute selection matrix in periodicty of x0
selectMatrix = eye(obj.nState);
for index = obj.periodicStatesHomotopy
    selectMatrix(index,index) = epsilon;
end

conPeriod = xT(obj.periodicStates)-selectMatrix(obj.periodicStates,:)*x0;

% add information to jacobian
add2JacobianPeriodicity();

%% define event constraints
nEvents   = length(obj.idxCon.Events);
conEvent  = zeros(nEvents,1);
eventData = cell(nEvents,1);
for iEvent = 1:nEvents
    % get x at event
    xe_i = simData.flowMapData.x_E(:,iEvent);
    te_i = sum(tDomain(1:iEvent));
    [conEvent(iEvent),eventData{iEvent}] = obj.model.getEventFunctional(te_i,xe_i,[],T,obj.sequence.events(iEvent),obj.controller,Grad1,Grad2);
    add2JacobianEvent(iEvent);
end
%% define functional constraints
nFun       = length(obj.idxCon.Additional);
conFun     = zeros(nFun,1);
fncs       = obj.functionals;

if ~isempty(fncs)
    idxInt = find(cellfun(@(x) x.int,fncs)); % index of integral-functionals
    
    % loop over all functionals!
    iConFun = 0; % counter of equality constraint functional
    for iFNC = 1:length(fncs)
        if obj.fixedParameter && iFNC==1
            % user-provided fixedParameter c must be in first functional!
            % get constant parameter from Solution-struct init
            c   = init.param.(fncs{iFNC}.c);
        else
            c   = 0;
        end
        
        if ~fncs{iFNC}.cost
            iConFun = iConFun +1; % increase equality constraint counter
            intPos = find(idxInt==iFNC); % find position of user-defined integral
            [conFun(iConFun),tau,tau_t,tau_x,tau_xi,tau_T,intT] = getFunctionalValue(iFNC);
            add2JacobianFunctional(iFNC,iConFun)
        end
    end
end

%% assamble constraints
fun = [conPeriod;...
       conEvent;...
       conFun];
   
%% get objective/cost data
if ~isempty(obj.functionals)
    idxCost = find(cellfun(@(x) x.cost,obj.functionals));
    if isempty(idxCost)
        % there is no objective/cost information for this problem
        objective = [];
    elseif length(idxCost)>1
        error('Multi-Cost-Functionals are not implemented yet!')
    else
        intPos = find(idxInt==idxCost); %find position of user-defined integral
        [objective.cost,tau,tau_t,tau_x,tau_xi,tau_T,intT] = getFunctionalValue(idxCost); 
        if Grad1
            objective.grad1 = getFunctionalGrad1(idxCost);
        end
    end
else
    objective = [];
end

% compute jacobian including Lagrange multipliers
if computeLxx == 1
    [fun,jacobian] = obj.getOptimJacobian2rootFunction(decVar,fun,jacobian,objective,init);
end

% if ~isempty(jacobian) && size(jacobian,2)~=length(x)
%     % this is not a continuation, thus a square problem
%     jacobian = jacobian(:,[obj.idxVarInEqCon,obj.idx.Multiplier]);
% end
   

%% nested functions 
%                   getFunctionalValue, getFunctionalGrad1
%                   add2JacobianPeriodicity, add2JacobianEvent, add2JacobianFunctional    
%
    function [F,tau,tau_t,tau_x,tau_xi,tau_T,intT] = getFunctionalValue(iFNC)
        if ~isempty(intPos)
            intT       = simData.flowMapData.integral_E(intPos,end);
        else
            intT = [];
        end
        
        switch fncs{iFNC}.type
            case 'start'
                [tau, tau_t, tau_x, tau_xi, tau_T] = obj.controller.inputTau(0, x0, z0, obj.model, T);
                F                                  = fncs{iFNC}.F(x0,tau,c,obj);
            case 'end'
                [tau, tau_t, tau_x, tau_xi, tau_T] = obj.controller.inputTau(T, xT, zT, obj.model, T);
                F                                  = fncs{iFNC}.F(xT,intT,tau,T,c,obj);  
        end
    end

    function grad = getFunctionalGrad1(iFNC)
        %allocate gradient
        grad = zeros(1,size(jacobian,2));
        functional = fncs{iFNC};
        switch functional.type
            case 'start'
                F_tau = functional.F_tau(x0,tau,c,obj);
                F_x   = functional.F_x(x0,tau,c,obj) + F_tau*tau_x;
                grad(obj.idx.FreeState) = F_x(obj.freeStates);
                grad(obj.idx.FreeCtrl)  = F_tau*tau_xi(:,obj.freeCtrls);
                % compute F_eps
                F_eps = functional.F_eps(x0,tau,c,obj);
                grad(obj.idx.Eps) = F_eps(obj.freeEps);

            case 'end'
                F_tau = functional.F_tau(xT,intT,tau,T,c,obj);
                F_x   = functional.F_x(xT,intT,tau,T,c,obj) + F_tau*tau_x;
                if isempty(functional.F_int)
                    F_int = [];
                else
                    F_int = functional.F_int(xT,intT,tau,T,c,obj);
                end
                
                F_x0  = F_x*xT_xE*xE.x0(:,:,end);
                if ~isempty(F_int)
                    F_x0 = F_x0 + F_int*intE.x0(intPos,:,end);
                end
                grad(obj.idx.FreeState) = F_x0(obj.freeStates);
                % get time derivatives
                % compute dF/dt_i
                F_T = functional.F_T(xT,intT,tau,T,c,obj);
                if obj.controller.timeBased
                    F_T = F_T + F_tau*(tau_T+tau_t);
                end
                for iTime = 1:obj.nTime
                   if ismember(iTime,obj.freeStates)
                       F_ti = F_T+F_x*xT_xE*xE.t(:,iTime,end);
                       if ~isempty(F_int)
                           F_ti = F_ti +F_int*intE.t(intPos,iTime,end);
                       end
                       grad(obj.idx.FreeTime(iTime)) = F_ti;
                   end
                end
                % compute dF/dxi
                F_xi  = F_tau*tau_xi + F_x*xT_xE*xE.xi(:,:,end);
                if ~isempty(F_int)
                    F_xi = F_xi + F_int*intE.xi(intPos,:,end);
                end
                grad(obj.idx.FreeCtrl) = F_xi(obj.freeCtrls);
               % compute dF/deps
               F_eps = functional.F_eps(xT,intT,tau,T,c,obj)...
                       + F_x*xT_xE*xE.eps(:,:,end)...
                       + F_x*simData.jumpMapData.epsGrad1(:,:,end);
               if ~isempty(F_int)
                   F_eps = F_eps + F_int*intE.eps(intPos,:,end);
               end
               grad(obj.idx.Eps) = F_eps(obj.freeEps);
        end
    end

    function add2JacobianPeriodicity()
        if ~Grad1
            return
        end
        
        % get size of periodic states
        nPer = length(obj.periodicStates);
        % compute dxT/dx0
        xT_x0 = xT_xE*xE.x0(:,:,end);
        jacobian(1:nPer,obj.idx.FreeState) = xT_x0(obj.periodicStates,obj.freeStates)-selectMatrix(obj.periodicStates,obj.freeStates);
        % compute dxT/dt_i
        for iTime = 1:obj.nTime
           if ismember(iTime,obj.freeTimes)
               xT_ti = xT_xE*xE.t(:,iTime,end);
               jacobian(1:nPer,obj.idx.FreeTime(iTime)) = xT_ti(obj.periodicStates);
           end
        end
        % compute dxT/dxi
        xT_xi = xT_xE*xE.xi(:,:,end);
        jacobian(1:nPer,obj.idx.FreeCtrl) = xT_xi(obj.periodicStates,obj.freeCtrls);
        
       % compute dxT/deps
       xT_eps = xT_xE*xE.eps(:,:,end)...
                + obj.flipOperator*simData.jumpMapData.epsGrad1(:,:,end);

       % d/deps selectMatrix*x0 = (d/deps selectMatrix) * x0
       % = diag(obj.periodicStatesHomotopy) * x0
       % or equivalently -> x0(~obj.periodicStatesHomotopy) = 0
       selMat_epsx0 = x0;
       selMat_epsx0(~ismember((1:obj.nState),obj.periodicStatesHomotopy)) = 0;
       % d xT / d eps - d (selectMatrix*x0) / d eps
       if ~isempty(obj.idx.Eps)
           jacobian(1:nPer,obj.idx.Eps) = xT_eps(obj.periodicStates)-selMat_epsx0(obj.periodicStates);
       end
    end

    function add2JacobianEvent(iEvent)
        if ~obj.options.Grad1
            return
        end

        % get size of periodic states
        nPer = length(obj.periodicStates);
        % compute de_i/dx0
        e_x0 = eventData{iEvent}.grad_x*xE.x0(:,:,iEvent);
        idx  = nPer+(1:(obj.nCon-nPer));
        jacobian(idx(iEvent),obj.idx.FreeState) = e_x0(obj.freeStates); 
        % compute de_i/dt_j
        for iTime = 1:obj.nTime
           if ismember(iTime,obj.freeTimes)
               e_ti = eventData{iEvent}.grad_x*xE.t(:,iTime,iEvent);
               if obj.controller.timeBased
                   e_ti = e_ti + eventData{iEvent}.grad_T;
                   if iEvent >= iTime
                       e_ti = e_ti + eventData{iEvent}.grad_t;
                   end
               end
               jacobian(idx(iEvent),obj.idx.FreeTime(iTime)) = e_ti;
           end
        end
        % compute de_i/dxi
        jacobian(idx(iEvent),obj.idx.FreeCtrl) = eventData{iEvent}.grad_x*xE.xi(:,obj.freeCtrls,iEvent) +...
                                             eventData{iEvent}.grad_xi(obj.freeCtrls);
             
        if ~isempty(obj.idx.Eps)
            % compute de_i/deps
            jacobian(idx(iEvent),obj.idx.Eps) = eventData{iEvent}.grad_x*xE.eps(:,:,iEvent) +...
                                                eventData{iEvent}.grad_eps;
        end
    end

    function add2JacobianFunctional(iFNC,iConFun)
        if ~obj.options.Grad1
            return
        end
        
        grad = getFunctionalGrad1(iFNC);
        jacobian(end-nFun+iConFun,:) = grad; % fill jacobian
    end

end