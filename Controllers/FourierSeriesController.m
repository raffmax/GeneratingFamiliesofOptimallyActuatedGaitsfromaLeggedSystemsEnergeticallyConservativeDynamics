classdef FourierSeriesController < Controller
    % FourierSeriesController is a general fourier-series controller
    %   This controller uses time based actuation, i.e. 
    %   tau = xi(1) + sum_i ( xi_i*sin(i*2*pi*t/T) + xi_i*cos(i*2*pi*t/T) )
    %
    %   Properties:
    %       nXi          number of control parameters
    %       model        reference to mechanical model
    %       xi           all control parameters
    %       instantData  instantaneous data ( M(q(t)), B(t), ...) at time t
    %       timeBased    this controller is time-based
    %       degree       degree of series sumation
    %       symmetric    use only odd terms in series (excluding the bias)
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
    % Last revision: 28-Oct-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        nXi       % number of free control parameters
        timeBased = true
        degree
        symmetric = false
    end

    methods (Access = public)
        %% - constructor -
        function obj = FourierSeriesController(model, xi, degree, symmetric)
            % Constructor
            obj.instantData.M = zeros(model.nQ);
            obj.instantData.B = zeros(model.nQ,model.nTau);
            obj.xi            = xi;
            obj.degree        = degree;
            
            if nargin > 3 % check if controller is for symmetric gait
                obj.symmetric = symmetric;
            end
            
            if obj.symmetric
                obj.nXi       = 2*degree*model.nTau;
            else
                obj.nXi       = (1+2*degree)*model.nTau;
            end
            
            if length(obj.xi) ~= obj.nXi
                error('The size of xi must be correspond to the degree of the Fourier series!')
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
            if obj.symmetric
                T      = 2*varargin{:};
            else
                T      = varargin{:};
            end
            w      = t/T;
            
            tau    = zeros(nTau,1);
            tau_w  = zeros(nTau,1);
            tau_x  = zeros(nTau,2*model.nQ);
            tau_xi = zeros(nTau,obj.nXi);
            
            for iTau = 1:nTau
                if ~obj.symmetric
                    % add offset
                    tau(iTau) = obj.xi((iTau-1)*(2*nDeg+1)+1);
                    tau_xi(iTau,(iTau-1)*(2*nDeg+1)+1) = 1;
                end
                
                for iDeg = 1:nDeg 
                    
                    if obj.symmetric
                        aIDX = (iTau-1)*(2*nDeg+1)+(iDeg-1)*2+1;
                        bIDX = (iTau-1)*(2*nDeg+1)+(iDeg-1)*2+2;
                        i    = iDeg*2-1;
                    else
                        aIDX = (iTau-1)*(2*nDeg+1)+1+(iDeg-1)*2+1;
                        bIDX = (iTau-1)*(2*nDeg+1)+1+(iDeg-1)*2+2;
                        i    = iDeg;
                    end
                    a    = obj.xi(aIDX);
                    b    = obj.xi(bIDX);
                    s    = sin(i*2*pi*w);
                    c    = cos(i*2*pi*w);
                    
                    tau(iTau)         = tau(iTau)   + a*c+b*s;
                    tau_w(iTau)       = tau_w(iTau) -(a*s-b*c)*i*2*pi;
                    tau_xi(iTau,aIDX) = c;
                    tau_xi(iTau,bIDX) = s;
                end
            end
            tau_t  = tau_w/T;
            %tau_tt = tau_ww/(T^2);
            if obj.symmetric
                tau_T  = tau_w*(-t/(T^2))*2;
            else
                tau_T  = tau_w*(-t/(T^2));
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