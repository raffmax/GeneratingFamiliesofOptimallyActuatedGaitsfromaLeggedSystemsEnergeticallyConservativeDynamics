function [h_xx,c_xx] = get2ndOrderDerivatives(obj,rootfun,x,h_x,c_x)
%getHessian returns second derivative of Lagrangian if
%useLagrangeMultiplier, else compute nullspace derivative

if obj.options.Grad2
    error('Analytic Hessian must be implemented first!')
else    
    % deactivate optimal control problem for gradient computation
    obj.setOptimalControl(false);
    obj.updateNumDecVar();
    if ~obj.options.Grad1
           [h_x,c_x] = get1stOrderDerivatives(obj,rootfun,x);    
    end
    
    [nF,nX] = size(h_x);
    % compute finite difference w.r.t. x_i
    s     = obj.options.FiniteDifferenceStepSize;
    h_xx = zeros(nF,nX,nX);
    c_xx = zeros(nX,nX);

    for ix = 1:nX
           iX       = zeros(nX,1);
           iX(ix)   = 1;
           [h_x_s,c_x_s] = get1stOrderDerivatives(obj,rootfun,x+s*iX);
           h_xx(:,:,ix)  = (h_x_s-h_x)/s;
           c_xx(ix,:)    = (c_x_s-c_x)/s;
    end
           
    %% reset
    % reset controller
    obj.controller.setFreeParameters(x(obj.idx.FreeCtrl));
    % reset epsilon
    if ~isempty(obj.idx.Eps)
        obj.model.setEpsilon(x(obj.idx.Eps));
    end
    obj.setOptimalControl(true);
    obj.updateNumDecVar();
end
end

function [h_x,c_x] = get1stOrderDerivatives(obj,rootfun,x)
if obj.options.Grad1
    [~,h_x,objective] = rootfun(x);
    c_x = objective.grad1;
else
    stepSize = obj.options.FiniteDifferenceStepSize;
    [~,h_x] = getJacobianFD(rootfun,x,stepSize);
    [~,c_x] = getJacobianFD(@(x) getCostValue(rootfun,x),x,stepSize);  
end
    function cost = getCostValue(rootfun,x)
            [~,~,objec] = rootfun(x);
            cost        = objec.cost;
    end
end
