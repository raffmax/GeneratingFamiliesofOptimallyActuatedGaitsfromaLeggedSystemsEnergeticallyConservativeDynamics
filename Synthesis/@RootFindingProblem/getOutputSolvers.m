function [sol,output] = getOutputSolvers(obj,init,x,output,jacobian,exitflag)
%getOutputSolvers 

[x0,tDomain,xi,epsilon,decVar] = obj.getStateTimeCtrlEpsFromDecVar(x,init);
T               = sum(tDomain); % Period Time

sol                        = Solution(x0,tDomain,xi,epsilon,init);
sol.rfpData.decVar         = x;
sol.rfpData.periodicStates = obj.periodicStates;

if ~isequal(obj.flipOperator, eye(obj.nState))
    sol.rfpData.flipOperator = obj.flipOperator;
else
    sol.rfpData.flipOperator = [];
end

if obj.fixedParameter
    field             = obj.functionals{1}.c;
    sol.param.(field) = init.param.(field);
end

% X: all variables excluding lambda
idxX     = [obj.idx.FreeState,obj.idx.FreeTime,obj.idx.FreeCtrl,obj.idx.Eps];
% x: all variables s.t. the optimization problem
idxx     = obj.idxVarInEqCon;

if isfield(output,'t')
    sol.rfpData.jacobianAug = [jacobian;output.t'];
end

[~,~,objective] = rootFunctionTDsingle(obj,init,decVar,2);

if ~isempty(objective)
    if isfield(objective,'cost')
        cost = objective.cost;
    else
        cost = [];
    end
    if isfield(objective,'cost')
        costGrad1 = objective.grad1(idxx)';
    else
        costGrad1 = [];
    end
else
    cost      = [];
    costGrad1 = [];
end
sol.rfpData.cost = cost;
output.cost      = cost;
sol.rfpData.OptimalCtrl.gradientCost = costGrad1;

if obj.optimalControl
    % get size of x = [x0,tDomain,xi]
    nx = length(idxx);
    nX = length(idxX);
    if obj.useLagrangeMultipliers
        multiplier  = decVar(obj.idx.Multiplier);
        if length(x)==length(decVar)
            jacobianCon = jacobian(1:nx,nX+1:end)';
            if length(obj.idx.Multiplier) ~= size(jacobianCon,1)
                % remove last row, since it is not part of the optimality
                % constraints
                jacobianCon(end,:) = [];
            end
            L_xx        = jacobian(1:nx,idxx);
        else
            jacobianCon = jacobian(1:nx,nx+1:end)';
            L_xx        = jacobian(1:nx,1:nx);
        end
        output.tCon = null(jacobianCon); %null space
        output.secondDerivative = output.tCon'*L_xx*output.tCon;
    else
        nCon        = obj.nCon;
        if length(x)==length(decVar)
            jacobianCon = jacobian(end-nCon+1:end,idxx);
            firstOpt_X  = jacobian(1:end-nCon,idxx);
            multiplier  = -jacobianCon'\costGrad1;
        else
            jacobianCon = jacobian(end-nCon+1:end,1:nx);
            firstOpt_X  = jacobian(1:end-nCon,1:nx);
            multiplier  = -jacobianCon'\costGrad1;
        end
        
        output.tCon = getTangent(jacobianCon); %null space
        output.secondDerivative = firstOpt_X*output.tCon;
    end
    output.jacobianCon      = jacobianCon;
    output.multiplier       = multiplier;
  
    sol.rfpData.OptimalCtrl.secondDerivative = output.secondDerivative;
    sol.rfpData.OptimalCtrl.jacobianCon      = output.jacobianCon;
    sol.rfpData.OptimalCtrl.tCon             = output.tCon;
    sol.rfpData.multiplier                   = multiplier;
    
else
    output.jacobianCon = jacobian;
    sol.rfpData.multiplier  = decVar(obj.idx.Multiplier); 
end

if exitflag~=1 || ~obj.options.Grad1
    return
end
%% construct Monodromy matrix and get cost
if strcmp(obj.integrationType,'TD') && strcmp(obj.shootingType,'single')
    fun  = @(x0,tDomain) obj.simulation.singleShootingTD(x0,tDomain,true,false);
else
    error('Has not been implemented yet!')
end

[~,simData] = fun(x0,tDomain);
% get eventGradients
nEvents         = length(obj.sequence.events);
eventGradients  = zeros(nEvents,obj.nState);
for iEvent = 1:nEvents
    % get x at event
    xe_i = simData.flowMapData.x_E(:,iEvent);
    te_i = sum(tDomain(1:iEvent));
    [~,eventData] = obj.model.getEventFunctional(te_i,xe_i,[],T,obj.sequence.events(iEvent),obj.controller,true,false);
    eventGradients(iEvent,:) = eventData.grad_x;
end

sol.dynamics.MonodromyMatrix = obj.simulation.computeMonodromyMatrix(simData,eventGradients);
sol.dynamics.f0              = simData.flowMapData.f_S(:,1);
sol.dynamics.E_x0            = obj.model.Grad1.E_x(epsilon,x0,obj.sequence.FlowMapOrder(:,1));

%% get eventVelocities
sol.rfpData.eventVelocities = diag(eventGradients*simData.flowMapData.f_E);

end
