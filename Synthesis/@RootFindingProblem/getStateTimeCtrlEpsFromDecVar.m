function [x0,tDomain,xi,epsilon,decVar] = getStateTimeCtrlEpsFromDecVar(obj,decVar,sol)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% if length(decVar) < obj.nDecVar
%     % add missing information from sol
%     varIdx = [obj.idx.FreeState,obj.idx.Time,obj.idx.Ctrl,obj.idx.Eps];
%     missingIndices = varIdx(~ismember(varIdx,obj.idxVarInEqCon));
%     if ~isempty(missingIndices)
%         if obj.idx.Time==missingIndices
%             decVar = [decVar(1:missingIndices-1);sol.tDomain;decVar(missingIndices:end)];
%         elseif obj.idx.Eps==missingIndices
%             decVar = [decVar(1:missingIndices-1);sol.epsilon;decVar(missingIndices:end)];
%         else
%             error('Please debug me!') % needs implementation
%         end
%     end
% end

x0                     = sol.x0;
x0(obj.freeStates)     = decVar(obj.idx.FreeState);
tDomain                = sol.tDomain;
tDomain(obj.freeTimes) = decVar(obj.idx.FreeTime);
xi                     = sol.xi;
xi(obj.freeCtrls)      = decVar(obj.idx.FreeCtrl);

% set new control parameters in Controller
obj.controller.setFreeParameters(xi);

% set/get new homotopy parameter in MechanicalModel
if isempty(obj.fixedEpsilon)
    epsilon = decVar(obj.idx.Eps);
else
    epsilon = sol.epsilon;
end
obj.model.setEpsilon(epsilon);

end

