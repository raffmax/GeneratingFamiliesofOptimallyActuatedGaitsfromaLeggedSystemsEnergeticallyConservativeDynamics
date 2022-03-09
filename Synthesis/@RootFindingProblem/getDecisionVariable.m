function decVar = getDecisionVariable(obj,Sol,varargin)
%getDecisionVariable returns the decision variable.
if strcmp(obj.integrationType,'TD') && strcmp(obj.shootingType,'single')
    decVar = [Sol.x0(obj.freeStates);Sol.tDomain(obj.freeTimes);Sol.xi(obj.freeCtrls);Sol.epsilon(obj.freeEps)];
    if obj.optimalControl
        if obj.useLagrangeMultipliers
            if nargin > 3
                % are Lagrange multipliers provided in varargin?
                lambda = varargin{:};
                if length(lambda)~=obj.nCon
                   error('Size of Lagrange multipliers must be equal to number of constraints!') 
                end
            elseif ~isempty(Sol.rfpData.multiplier)
                lambda = Sol.rfpData.multiplier;
            else
                fun = @(x)rootFunctionTDsingle(obj,Sol,x,0);
                % assume optimal point and reconstruct lambda
                [F,jacobian,objective] = fun([decVar;nan*eye(size(obj.idx.Multiplier'))]);
                % lambda is unique if decVar is optimal point and jacobian has 
                % full rank 
                idxCon = obj.idxEqConInOptim-obj.idxEqConInOptim(1)+1;
                lambda = -jacobian(idxCon,obj.idxVarInEqCon)'\objective.grad1(obj.idxVarInEqCon)'; 
                if max(abs(F))>obj.options.FunctionTolerance || rank(jacobian)<obj.nCon
                   warning('Reconstructing Lagrange multipliers may not be accurate!') 
                end
            end   
        else
            lambda = [];
        end
        decVar = [decVar;lambda];
    end
elseif strcmp(obj.shootingType,'multiple')
    error('multiple shooting is not implemented yet')
end
end

