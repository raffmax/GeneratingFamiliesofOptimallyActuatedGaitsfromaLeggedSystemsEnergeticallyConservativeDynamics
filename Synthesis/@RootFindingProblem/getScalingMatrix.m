function scalingMatrix = getScalingMatrix(obj,useLagrangeVar,fixEpsilon)
%getScalingMatrix gets the Scaling matrix for decVar in RootFindingProblem
if nargin < 3
    fixEpsilon = false;
end

if fixEpsilon
    nEps = 0;
else
    nEps = length(obj.idx.Eps);
end
if useLagrangeVar
    scalingMatrix = diag([ones(1,obj.nFreeState),ones(1,obj.nTime),ones(1,obj.nCtrl)*obj.scaleControl,ones(1,length(obj.idx.Multiplier))*obj.scaleLagrange,ones(1,nEps)]);
else
    scalingMatrix = diag([ones(1,obj.nFreeState),ones(1,obj.nTime),ones(1,obj.nCtrl)*obj.scaleControl,ones(1,nEps)]);
end
end
