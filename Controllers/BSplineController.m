classdef BSplineController < Controller
    % BSplineController is a general controller using periodic B-Splines
    %   This controller uses time based actuation, i.e. 
    %   tau = sum_i x_i*B_(i,d)(t)
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
    % Last revision: 20-September-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        nXi       % number of free control parameters
        timeBased = true
        degree
        segments
        symmetric = false
        periodic = true
    end

    methods (Access = public)
        %% - constructor -
        function obj = BSplineController(model, xi, degree, segments, symmetric, periodic)
            % Constructor
            obj.instantData.M = zeros(model.nQ);
            obj.instantData.B = zeros(model.nQ,model.nTau);
            obj.xi            = xi;
            obj.segments      = segments;
            obj.degree        = degree;
            
            if nargin > 4 % check if controller is for symmetric gait
                obj.symmetric = symmetric;
                if nargin == 6 % check if tau needs to be periodic
                    obj.periodic  = periodic;
                end
            end
            
            if obj.periodic
                obj.nXi = model.nTau*segments;
                if length(obj.xi) ~= obj.nXi
                    error('The size of xi must be equal to the number of segments!')
                end
            else
                obj.nXi = model.nTau*segments+obj.degree;
                if length(obj.xi) ~= obj.nXi
                    error('The size of xi must be equal to the number of segments + degree of splines!')
                end                
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
            nSeg   = obj.segments; %number of segments
            nCtrl  = nSeg+nDeg; % number of control points
            dKnot  = 1/nSeg;
            % construct knot vector
            knotVec = -nDeg*dKnot:dKnot:1+nDeg*dKnot;

            tau    = zeros(nTau,1);
            tau_w  = zeros(nTau,1);
            % tau_ww = zeros(nTau,1);
            tau_x  = zeros(nTau,2*model.nQ);
            tau_xi = zeros(nTau,obj.nXi);

            xi_idx = 0; %index xi
            useSym = false; % use symmetry in periodicity
            for i = 0:nCtrl-1 
                % compute base spline
                B = baseSplineB(w,knotVec,i,nDeg);
                % comput deriv. of base spline
                if nDeg>0
                    BPrime = diffB(1,w,knotVec,i,nDeg);
                end
                if i == nSeg && obj.periodic% reset for periodicity
                    xi_idx = 0;
                    if obj.symmetric
                        useSym = true;
                    end
                end
                xi_idx  = xi_idx+1;

                for iTau = 1:nTau
                    if useSym
                    	xi = -obj.xi(xi_idx+(iTau-1)*nSeg);
                    else
                        xi = obj.xi(xi_idx+(iTau-1)*nSeg);
                    end
                    tau(iTau) = tau(iTau) + xi*B;

                    if nDeg>0
                        tau_w(iTau) = tau_w(iTau) + xi*BPrime;
                    end
                    if useSym
                        tau_xi(iTau,xi_idx+(iTau-1)*nSeg) = tau_xi(iTau,xi_idx+(iTau-1)*nSeg) - B;
                    else
                        tau_xi(iTau,xi_idx+(iTau-1)*nSeg) = tau_xi(iTau,xi_idx+(iTau-1)*nSeg) + B;
                    end
                end
            end
          
            tau_t  = tau_w/T;
            %tau_tt = tau_ww/(T^2);
            tau_T  = tau_w*(-t/(T^2));

            function BPrime = diffB(order,t,knotVec,i,k)
                denom1 = knotVec(i+1+k)-knotVec(i+1);
                denom2 = knotVec(i+2+k)-knotVec(i+2);
                if denom1==0
                    B_term1 = 0;
                else
                    if order == 1
                        B_term1 = k*baseSplineB(t,knotVec,i,k-1)/denom1;
                    else
                        B_term1 = k*diffB(order-1,t,knotVec,i,k-1)/denom1;
                    end
                end
                if denom2==0
                    B_term2 = 0;
                else
                    if order == 1
                        B_term2 = -k*baseSplineB(t,knotVec,i+1,k-1)/denom2;
                    else
                        B_term2 = -k*diffB(order-1,t,knotVec,i+1,k-1)/denom2;
                    end
                end
               BPrime = B_term1+B_term2;
            end

            function B = baseSplineB(t,knotVec,i,k)
                % Cox–de Boor recursion 
                if k==0
                    if t>=knotVec(i+1) && t<knotVec(i+2)
                        B=1;
                    else
                        B=0;
                    end
                else
                    denom1 = knotVec(i+1+k)-knotVec(i+1);
                    if denom1==0
                        B_term1 = 0;
                    else
                        B_term1 = (t-knotVec(i+1))/denom1*baseSplineB(t,knotVec,i,k-1);
                    end
                    denom2 = knotVec(i+2+k)-knotVec(i+2);
                    if denom2 == 0
                        B_term2 = 0;
                    else
                        B_term2 = (knotVec(i+2+k)-t)/denom2*baseSplineB(t,knotVec,i+1,k-1);
                    end
                    B = B_term1+B_term2;
                end
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