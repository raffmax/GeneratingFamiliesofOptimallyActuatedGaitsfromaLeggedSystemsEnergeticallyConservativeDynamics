classdef ZeroContoller < Controller
    % ZeroContoller is a pseudo-controller that applies zero actuation
    %   ZeroContoller applies constant zero actuation to the mechanical
    %   model. Note, the default controller for conservative systems is 
    %   ConservativeController.
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
    %  see also: Controller, ConservativeController
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 13-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        nXi       = 0       % number of free control parameters
        timeBased = false
    end
    methods (Access = public)
        %% - constructor -
        function obj = ZeroContoller(model, ~)
            obj.instantData.M = zeros(model.nQ);
            obj.instantData.B = zeros(model.nQ,1);
            obj.xi            = [];
        end
        %% - controller function, generates the input for the system -
        function [u,u_t,u_x,u_xi,u_T] = inputU(~, ~, ~, ~, model, varargin)
            % u -> xDot = f(x) + u
            u       = zeros(2*model.nQ,1);
            u_x     = zeros(2*model.nQ);
            u_xi    = zeros(2*model.nQ,1);

            u_t     = [];
            u_T     = [];
        end
        
        function [tau, tau_t, tau_x, tau_xi, tau_T] = inputTau(obj, ~, ~, ~, model, varargin)
            % tau -> M*ddq = B*tau + ...
            nTau     = model.nTau;
            tau      = zeros(nTau,1);
            tau_x    = zeros(nTau,2*model.nQ);
            tau_xi   = zeros(nTau,obj.nXi);
            
            tau_t = [];
            tau_T = [];
        end
 
        
        function param = getFreeParameters(obj)
            % set xi
            param = obj.xi;
        end

        function setFreeParameters(obj, ~)
            % get xi
            obj.xi = [];
        end
    end
end