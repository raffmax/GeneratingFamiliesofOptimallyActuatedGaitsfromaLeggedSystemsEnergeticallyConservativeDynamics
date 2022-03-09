function nDecVar = getNumDecVar(obj)
%getNumDecVar counts decision variables
%   returns number of decision variables of RootFindingProblem

nDecVar = 0;
% get number of free states
if strcmp(obj.shootingType,'single')
    nDecVar = nDecVar + obj.nFreeState;
elseif strcmp(obj.shootingType,'multiple')
    error('multiple shooting is not implemented yet')
end
% get number of domain times
nDecVar = nDecVar + obj.nFreeTime;
% get number of control parameters
nDecVar = nDecVar + obj.nFreeCtrl;
% get number of variable homotopy/model parameters
nDecVar = nDecVar + obj.nFreeEps;

if isempty(obj.nDecVar)
    % default when RootFindingProblem is initiated
    obj.idxVarInEqCon = 1:nDecVar;
end
% get number of lagrange multipliers
if obj.optimalControl && obj.useLagrangeMultipliers
	nDecVar = nDecVar + length(obj.idxEqConInOptim);
end

% set decVar index
obj.idx.FreeState = (1:obj.nFreeState);
obj.idx.FreeTime  = obj.nFreeState+(1:obj.nFreeTime);
obj.idx.FreeCtrl  = obj.nFreeState+obj.nFreeTime+(1:obj.nFreeCtrl);
if obj.nFreeEps==0
    obj.idx.Eps   = [];
else
    obj.idx.Eps   = obj.nFreeState+obj.nFreeTime+obj.nFreeCtrl + obj.nFreeEps;
end
if obj.optimalControl && obj.useLagrangeMultipliers
    obj.idx.Multiplier = obj.nFreeState+obj.nFreeTime+obj.nFreeCtrl+obj.nFreeEps+(1:length(obj.idxEqConInOptim));
else
    obj.idx.Multiplier = [];
end
end

