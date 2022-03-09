function [funOpt,jacobianOpt] = getOptimJacobian2rootFunction(obj,decVar,fun,jacobian,objective,init)
%getOptimJacobian2rootFunction computes rootFun of constrained optimization
%   uses Lagrange-Multipliers or Null-Space projections
%   is called in rootFunctionTDsingle()
h = fun; % equality constraints
% X: all variables excluding lambda
idxX     = [obj.idx.FreeState,obj.idx.FreeTime,obj.idx.FreeCtrl,obj.idx.Eps];
% x: all variables s.t. the optimization problem
idxx        = obj.idxVarInEqCon;
X           = decVar(idxX);
h_X         = jacobian(:,idxX);
h_x         = jacobian(:,idxx);
c_X         = objective.grad1(idxX);
c_x         = objective.grad1(idxx);
[h_XX,c_XX] = obj.get2ndOrderDerivatives(@(X)rootFunctionTDsingle(obj,init,X,2),X,h_X,c_X);
if obj.useLagrangeMultipliers
    % L   = c(X) + lambda'*h(X)
    % fun = [L_X,L_lambda]' = [L_X';h]
    
    idxCon   = obj.idxEqConInOptim-obj.idxEqConInOptim(1)+1;
    Lambda   = decVar(obj.idx.Multiplier);
    
    L_X         = c_X'+h_X(idxCon,:)'*Lambda;
    h_XX_Transp = permute(h_XX(idxCon,:,:), [2 1 3]);
    L_XX        = c_XX + t3m(h_XX_Transp,Lambda);

    L_x = L_X(idxx);

    funOpt      = [L_x;h];
    jacobianOpt = [L_XX(idxx,:), h_x(idxCon,:)'; ...
                   h_X ,zeros(obj.nCon,length(idxCon))];
    %L_XX2 = obj.getHessian(@(X)rootFunctionTDsingle(obj,init,X,2),X,Lambda);            
else
    % fun = [null(h'(x))*c'(x),h(x)]'
    D           = getTangent(h_x); %null space of h_x
    firstOpt    = D'*c_x';
    funOpt      = [firstOpt;h];
    
    
% %%    
%     % deactivate optimal control problem for gradient computation
%     obj.setOptimalControl(false);
%     obj.updateNumDecVar();
%     [nF,nX] = size(D);
%     % compute finite difference w.r.t. x_i
%     s     = obj.options.FiniteDifferenceStepSize;
%     D_X   = zeros(nF,nX,obj.nDecVar);
% 
%     for ix = 1:obj.nDecVar  
%            iX       = zeros(obj.nDecVar  ,1);
%            iX(ix)   = 1;
%            [~,h_X_s] = rootFunctionTDsingle(obj,init,X+s*iX,2);
%            D_s       = getTangent(h_X_s(:,idxx));
%            D_X(:,:,ix)  = (D_s-D)/s;
%     end
%     
%     %% reset
%     % reset controller
%     obj.controller.setFreeParameters(decVar(obj.idx.Ctrl));
%     % reset epsilon
%     if ~isempty(obj.idx.Eps)
%         obj.model.setEpsilon(decVar(obj.idx.Eps));
%     end
%     obj.setOptimalControl(true);
%     obj.updateNumDecVar();
% %%
% 
%     D_X_Transp  = permute(D_X, [2 1 3]);
%     firstOpt_X1  = t3m(D_X_Transp,c_x')+D'*c_XX(idxx,:);   
%     
%     firstOpt_X   = obj.getHessian(@(X)rootFunctionTDsingle(obj,init,X,2),X,[]); 
    h_xX_D      = t3m(h_XX(:,idxx,:),D);  
    D_X_Transp  = -t3m(permute(h_xX_D,[2 1 3]),pinv(h_x)');
    firstOpt_X  = t3m(D_X_Transp,c_x')+D'*c_XX(idxx,:); 
    

    jacobianOpt = [firstOpt_X; h_X];
end
end

