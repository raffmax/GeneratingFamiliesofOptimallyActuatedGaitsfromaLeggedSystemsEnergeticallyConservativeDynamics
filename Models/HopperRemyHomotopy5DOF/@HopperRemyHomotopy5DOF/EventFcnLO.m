function [fun, grad_t, grad_x, grad_xi, grad_eps, grad_T] = EventFcnLO(obj, t, x, T, controller, Grad1, Grad2)
%EventFcnLO lift-off event of hopper 5DOF
%   Event Lift-Off is defined kinetically by the scaled normal force
%   The gradients of fun w.r.t. x, xi and eps are also provided if
%   Grad1=true

z       = [1;0]; % in stance before lift-off
epsilon = obj.getEpsilon();
k_l     = obj.parameters.param_k_l(epsilon);
b_l     = obj.parameters.param_b_l(epsilon);
l_0     = obj.parameters.l_0;

[~,lambda,~, lambdaGrad] = obj.FlowMap(t, x, z, T, controller, Grad1, Grad2);
[tau, tau_t, tau_x, tau_xi, tau_T] = controller.inputTau( t, x, z, obj, T);

F_l = k_l*(l_0-x(5)) - b_l*x(10) + tau(2); % joint force

fun = epsilon*lambda(2) + (1-epsilon)*F_l*cos(x(3)+x(4));

if ~Grad1 
    grad_x   = 0;
    grad_xi  = 0;
    grad_eps = 0;
    grad_t   = 0;
    grad_T   = 0;
else
    F_l_x   = [0,0,0,0,-k_l, 0,0,0,0,-b_l]+tau_x(2,:);
    F_l_xi  = tau_xi(2,:);
    
    grad_x  = epsilon*lambdaGrad.lambda_x(2,:)+(1-epsilon)*F_l_x*cos(x(3)+x(4))...
              -(1-epsilon)*F_l*sin(x(3)+x(4))*[0 0 1 1 zeros(1,6)];
    grad_xi = epsilon*lambdaGrad.lambda_xi(2,:)+(1-epsilon)*F_l_xi*cos(x(3)+x(4));

    F_l_eps  = obj.parameters.grad_k_l(epsilon)*(l_0-x(5))-obj.parameters.grad_b_l(epsilon)*x(10);
    grad_eps = epsilon*lambdaGrad.lambda_eps(2,:)+lambda(2)...
               +(1-epsilon)*F_l_eps*cos(x(3)+x(4)) -F_l*cos(x(3)+x(4));

    if controller.timeBased
        grad_t   = epsilon*lambdaGrad.lambda_t(2,:)...
                   +(1-epsilon)*tau_t(2,:)*cos(x(3)+x(4));
        grad_T   = epsilon*lambdaGrad.lambda_T(2,:)...
                   +(1-epsilon)*tau_T(2,:)*cos(x(3)+x(4));
    else
        grad_t   = 0;
        grad_T   = 0;
    end 
end
end

