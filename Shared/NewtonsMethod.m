function[x,fval,exitflag,output,jacobian] = NewtonsMethod(fun,x,options,varargin)
%NewtonsMethod implements Newton's method.
%   NewtonsMethod can be called by the corrector-step of a continuation
%   method or an arbitrary root-finding-problem
%   see also getJacobianFD

if nargin>3
    f        = varargin{1};
    jacobian = varargin{2};
    t        = varargin{3}; % nullspace of jacobian / tangent space
    
    funCounter = 0; % counter for function evaluations
    
    if options.Grad1
        useFD = false; % use finite differences
    else
        useFD    =  true;
        stepSize = options.FiniteDifferenceStepSize;
        % get dimensions of jacobian
        nF    = length(f);
        nX    = length(x);
    end
    
else
    t = []; % tangent space is assumed to be empty
    
    % this option is only important in Continuation-Corrector
    options.aimOnTarget = false;
    
   if options.Grad1
       useFD = false;
       [f,jacobian] = fun(x);
       funCounter   = 1; % counter for function evaluations
   else
       useFD        = true;
       stepSize     = options.FiniteDifferenceStepSize;
       [f,jacobian] = getJacobianFD(fun,x,stepSize);
       % get dimensions of jacobian
       nF = length(f);
       nX = length(x);
       funCounter   = 1 + nX*nF; % counter for function evaluations
   end   
end

err  = max(abs(f));
iter = 0;
while err > options.FunctionTolerance && iter < options.MaxIterations  
   % newton step
   if options.aimOnTarget
        idx = options.idxConPar;
        delta_x = -[jacobian;[zeros(1,idx-1),1,zeros(1,length(x)-idx)]]\[f;0]; % fix lambda
   else
        delta_x = -[jacobian;t']\[f;zeros(size(t,2))]; % minimze distance to curve if size(t,2)=1 / equivalent to Moore-Penrose Inverse
   end

   if norm(delta_x,1) < options.StepTolerance 
       warning(['The Newton step is smaller than StepTolerance ',num2str(options.StepTolerance),'!'])
       break
   end

   if options.LineSearch
       % update with line search (Globalization of Newton)
       alpha     = 1; % initial step length
       alphaFlag = false; % is step length alpha excepted

       % compute Merit function of problem
       [M,M_x] = computeMeritFunction(f,jacobian);
       while ~alphaFlag
           x_alpha = alpha*delta_x + x;
           if useFD
               [f_alpha,jacobian_alpha] = getJacobianFD(fun,x_alpha,stepSize);
               funCounter               = funCounter + 1 + nX*nF; 
           else
               [f_alpha,jacobian_alpha] = fun(x_alpha);
               funCounter               = funCounter + 1;
           end
           M_alpha = computeMeritFunction(f_alpha,jacobian_alpha);

           if alpha == 1
               x_alpha1        = x_alpha;
               f_alpha1        = f_alpha;
               jacobian_alpha1 = jacobian_alpha;
           end

           % check Goldstein-Armijo Condition (Betts 2010 - Chapter 1.11.2)
           kappa1  = 1e-4;
           kappa2  = 0.9;
           c_alpha = alpha*M_x*delta_x;
           c_1     = -kappa1*c_alpha;
           c_2     = M-M_alpha;
           c_3     = -kappa2*c_alpha;
           if ((0<c_1) && (c_1<=c_2) && (c_2<= c_3)) || norm(alpha*delta_x)^2 < options.StepTolerance 
               % except step length
               alphaFlag = true;
               % update x
               if norm(alpha*delta_x)^2 < options.StepTolerance 
                   x         = x_alpha1;
                   f         = f_alpha1;
                   jacobian  = jacobian_alpha1;
               else
                   x         = x_alpha;
                   f         = f_alpha;
                   jacobian  = jacobian_alpha;
               end
           else
               % update alpha
               alpha = alpha/2;
           end
       end
   else
       % Newton update
       x = delta_x + x;
       if useFD
           [f,jacobian] = getJacobianFD(fun,x,stepSize);
           funCounter               = funCounter + 1 + nX*nF; 
       else
           [f,jacobian] = fun(x);
           funCounter               = funCounter + 1;
       end
   end
   
   % compute new tangent vector
   t = getTangent(jacobian);

   iter = iter + 1;
   err  = max(abs(f));       
end  
fval   = f;
if err < options.FunctionTolerance
   exitflag = 1;
else
   exitflag = -1;
end
output.t          = t;
output.iterations = iter;
output.funcCount  = funCounter;

end

function [M,M_x] = computeMeritFunction(f,f_x)
M   = 0.5*(f')*f;
if nargout>1
    M_x = f'*f_x; 
else
    M_x = [];
end
end