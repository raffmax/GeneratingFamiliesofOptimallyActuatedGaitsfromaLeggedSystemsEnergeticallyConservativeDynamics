classdef HopperRemyHomotopy5DOF < MechanicalModel
    % HopperRemyHomotopy5DOF represents a mechanical model with a model homotopy.
    %   This model represents a pogo-stick constrained to planar motions.
    %   Thus, the main mass is described by vertical position y, 
    %   orientation phi and its relative foot position position by l. 
    %   The model homotopy connects a 2D SlIP model at
    %   epsilon = 0 and a 2D pogo-stick at epsilon = 1.
    %
    %   Properties:  
    %      nZ            number of (physical) constraints
    %      nQ            number of generalized coordinates
    %
    %      more properies are defined in super-class
    %
    %   Methods:
    %      Hopper2Dpitch       The constructor for this class
    %      EventMap1           computes touch-down event
    %      EventMap2           computes lift-off event
    %      getEventFunctional  decideds which EventMap to call
    %      CreateEOM           autogenerate EOM
    %
    %   Example:
    %       hopper = Hopper2Dpitch('load', false);
    %
    %   See also MechanicalModel

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 12-July-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Constant)
        nZ   = 2       % number of (physical) constraints
        nQ   = 5       % number of generalized coordinates
        nTau = 2       % number of actuators
        
        nLambda     = [2 1];  % size of constraint forces
        homotopyCon = [0 1];  % homotopy in "physical" constraints
        % z = [z_S z_F]; stance, flight
    end
 
    methods
        % - constructor - 
        function obj = HopperRemyHomotopy5DOF(varargin)
            
            % default values for the configuration
            default_coefRes  = 0;
            default_load_    = false;
            default_id       = 'hopper5dof';
            default_epsilon  = 0;
            
            % default path to cache of model
            main_dir = which('install_project.m'); % main dir is where install_project is located
            [main_dir ,~,~] = fileparts(main_dir);
            default_path    = fullfile(main_dir, 'AppData');
            
            p = inputParser;
            % cfg 
            addParameter(p, 'coefRes',  default_coefRes,  @isnumeric);
            addParameter(p, 'epsilon',  default_epsilon,  @isnumeric);
            addParameter(p, 'load',     default_load_,    @islogical);
            addParameter(p, 'path',     default_path,     @ischar);
            addParameter(p, 'ID',       default_id,       @ischar);
            
            
            parse(p, varargin{:});
            
            load_        = p.Results.load;
            coefRes      = p.Results.coefRes;
            obj.epsilon = p.Results.epsilon;

            obj.nEps    = 1; % number of homotopy parameters

            % place where to 'deserialize' the object
            obj.identifier   = p.Results.ID;
            obj.folder_path  = fullfile(p.Results.path, obj.identifier); % default folder is AppData/ID
            obj.load_path    = fullfile(obj.folder_path, [obj.identifier, '.mat']); % .mat file of obj
            
            if load_
                disp(strcat('Loading ', " ",  obj.identifier, " {2D Hopper}..."))
                tmp  = load(obj.load_path);
                obj = tmp.obj;
                addpath(obj.folder_path)
                disp('done!')
            else
                deleteAppData(obj.identifier); % ask if any cache should be deleted
                
                disp('Creating 2D Hopper ...')
                % parameters
                obj.parameters = HopperParametersRemyHomotopy5DOF(coefRes);
                
                % - compute EOM function handles - 
                obj.CreateHybridDynamics();
                
                % save obj as .mat for loading later
                save(obj.load_path, 'obj');  
                disp('done!')
            end
        end
        
        
        
        %% --- Flow Map ---
        %[f, lambda, fGrad]               = FlowMap(obj, t, x, z, controller, Grad1, Grad2, epsGrad, varargin)
        %% --- Jump Map ---
        %[xPlus, Lambda, jumpGrad]        = JumpMap(obj, xMinus, z_next, Grad1, Grad2, epsGrad)
        %% --- Event Maps ---
        % touch down event
        [fun, grad_t, grad_x, grad_xi, grad_eps, grad_T] = EventFcnTD(obj, t, x, T, controller, Grad1, Grad2) % touch down event
        % lift-off event
        [fun, grad_t, grad_x, grad_xi, grad_eps, grad_T] = EventFcnLO(obj, t, x, T, controller, Grad1, Grad2) % lift-off event
        %apax
        [fun, grad_t, grad_x, grad_xi, grad_eps, grad_T] = EventFcnApax(obj, t, x, T, controller, Grad1, Grad2) % apax event
        
    end

    methods %% --- Misc ----
        function [conEvent,eventData] = getEventFunctional(obj,t,x,z,T,iEvent,controller, Grad1, Grad2)
            % select event functionals 
            switch iEvent
                case 1 % touch-down
                        [conEvent, grad_t, grad_x, grad_xi, grad_eps, grad_T] = obj.EventFcnTD(t, x, T, controller, Grad1, Grad2);
                        eventData.grad_t    = grad_t;
                        eventData.grad_x    = grad_x;
                        eventData.grad_xi   = grad_xi;
                        eventData.grad_eps  = grad_eps;
                        eventData.grad_T    = grad_T;
                case 2 % lift-off
                        [conEvent, grad_t, grad_x, grad_xi, grad_eps, grad_T] = obj.EventFcnLO(t, x, T, controller, Grad1, Grad2);
                        eventData.grad_t    = grad_t;
                        eventData.grad_x    = grad_x;
                        eventData.grad_xi   = grad_xi;
                        eventData.grad_eps  = grad_eps;
                        eventData.grad_T    = grad_T;
                case 3 % apax
                        [conEvent, grad_t, grad_x, grad_xi, grad_eps, grad_T] = obj.EventFcnApax(t, x, T, controller, Grad1, Grad2);
                        eventData.grad_t    = grad_t;
                        eventData.grad_x    = grad_x;
                        eventData.grad_xi   = grad_xi;
                        eventData.grad_eps  = grad_eps;
                        eventData.grad_T    = grad_T;
            end
        end
        
        function obj = CreateEOM(obj)  
            % call script to autogenerate dynamics and function handles
            obj = generateEOMhopperRemyHomotopy5DOF(obj);
        end
    end
end

