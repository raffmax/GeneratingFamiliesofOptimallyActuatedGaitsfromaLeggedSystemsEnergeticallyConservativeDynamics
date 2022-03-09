classdef BezierController < Controller
    % BezierController is a general controller using Bezier curves
    %   This controller uses time based actuation, i.e. 
    %   tau = xi(1) + sum_i ( xi_i*sin(i*2*pi*t/T) + xi_i*cos(i*2*pi*t/T) )
    %
    %   Properties:
    %       nXi          number of control parameters
    %       model        reference to mechanical model
    %       xi           all control parameters
    %       instantData  instantaneous data ( M(q(t)), B(t), ...) at time t     
    %
    %   Methods:
    %       inputU       computes general input u in xDot = f(x) + u
    %       inputTau     computes torque input tau in M*ddq = B*tau + ...
    %       update       updates M and B at time instance t
    %
    %       getFreeParameters   get free ctrl parameter (decision variable)
    %       setFreeParameters   set free ctrl parameter (decision variable)
    %
    %   see also: Controller
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 15-September-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        nXi       % number of free control parameters
        timeBased = true
        degree
    end

    methods (Access = public)
        %% - constructor -
        function obj = BezierController(model, xi, degree)
            % Constructor
            obj.instantData.M = zeros(model.nQ);
            obj.instantData.B = zeros(model.nQ,model.nTau);
            obj.xi            = xi;
            obj.nXi           = degree*model.nTau;
            obj.degree        = degree;
            
            if length(obj.xi) ~= obj.nXi
                error('The size of xi must correspond to the degree of the closed Bezier curve!')
            end
        end
        %% - controller function, generates the input for the system -
        function [u,u_t,u_x,u_xi,u_T] = inputU(obj, ~, ~, ~, model, varargin)
            % u -> xDot = f(x) + u
            u       = zeros(2*model.nQ,1);
            u_t     = zeros(2*model.nQ,1);
            u_T     = zeros(2*model.nQ,1);
            u_x     = zeros(2*model.nQ);
            u_xi    = zeros(2*model.nQ,obj.nXi);
        end
        
        function [tau, tau_t, tau_x, tau_xi, tau_T] = inputTau(obj, t, ~, ~, model, varargin)
            % tau -> M*ddq = B*tau + ...
            nTau   = model.nTau;
            nDeg   = obj.degree;
            T      = varargin{:};
            w      = t/T;
            
            tau    = zeros(nTau,1);
            tau_w  = zeros(nTau,1);
            tau_x  = zeros(nTau,2*model.nQ);
            tau_xi = zeros(nTau,obj.nXi);
            
            for iDeg = 0:nDeg  
                B_i  = getBernsteinBasis(w,iDeg,nDeg); % Bernstein basis 
                % derivative
                if iDeg == 0
                    BPrime = -getBernsteinBasis(w,iDeg,nDeg-1);
                elseif iDeg == nDeg
                    BPrime = getBernsteinBasis(w,iDeg-1,nDeg-1);
                else
                    BPrime = getBernsteinBasis(w,iDeg-1,nDeg-1)-getBernsteinBasis(w,iDeg,nDeg-1);
                end            

                for iTau = 1:nTau
                    if iDeg < nDeg
                        idx_xi = (iTau-1)*nDeg+(iDeg+1);
                    else
                        idx_xi = (iTau-1)*nDeg+1;
                    end
                    xi_i = obj.xi(idx_xi);

                    tau(iTau)           = tau(iTau) + B_i*xi_i;
                    tau_w(iTau)         = tau_w(iTau) + nDeg*BPrime*xi_i;
                    tau_xi(iTau,idx_xi) = tau_xi(iTau,idx_xi)+B_i;
                end
            end

            tau_t = tau_w/T;
            tau_T = tau_w*(-t/(T^2));
            
            function B_i = getBernsteinBasis(t,i,n)
                B_i  = nchoosek(n,i)*(t^(i))*((1-t)^(n-i));
            end
        end
 
        
        function param = getFreeParameters(obj)
            % set xi
            param = obj.xi;
        end

        function setFreeParameters(obj, param)
            % get xi
            obj.xi = param;
        end

    end
    
end