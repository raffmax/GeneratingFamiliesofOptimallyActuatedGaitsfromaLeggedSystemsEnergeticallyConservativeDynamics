classdef RootFindingProblem < matlab.mixin.Copyable
    % RootFindingProblem defines the general structure of the root finding problem
    %   This 
    %
    %   Properties:  
    %
    %
    %   Methods:
    %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 20-Oct-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (SetAccess = private) % --- Class-References ---
        model          % reference to mechanical model                                     
        controller     % reference to used controller
        simulation     % reference to simulation object                                                      
    end
    
    properties (SetAccess = private) % --- Model dependent props ---
        freeStates    % states that are part of decision variables
        freeTimes     % time durations that are part of decision variables
        freeCtrls     % control parameters that are part of decision variables
        freeEps       % epsilon parameter that is part of decision variables
        
        periodicStates % states on which we enforce periodicity with period T
        periodicStatesHomotopy % states that start being fixed (epsilon=0) and transition to periodicity (epsilon=1)
        
        fixedStates    % states that are not part of decision variables
        fixedTimes     % time durations that are not part of decision variables
        fixedCtrls     % control parameters that are not part of decision variables
        fixedEpsilon   % epsilon is not part of decision variables

        sequence %Predefined Footfall pattern                                                         

        flipOperator % exploit symmetry or use coordinate transformation
        % - active implicit equality constraint (e.g. fix initial system energy)
        fixedParameter
        % - functionals for implicit equality constraints or costs
        functionals
        % - define optimization Problem in Control-Parameters
        optimalControl
        % - use Lagrange Mulipliers in Optimization/Continuation
        useLagrangeMultipliers
%         % -  scale variables
%         scaleControl
%         scaleLagrange
    end

    properties % --- Solver/Implementation dependent props ---
        integrationType % event-driven (ED), time-driven (TD), 
        shootingType    % single shooting (single), multiple shooting (multiple)
        nDecVar         % number of Decision Variables
        nCon            % number of Constraints
        nInt            % number of additional integrals
        idx             % struct indexing decVar for free-state, time and param 
        idxCon          % struct indexing equality constraints
        idxVarInEqCon   % variables that are part of the optimization problem
        idxEqConInOptim % equality constraints that are part of the optimization problem
        nState          % number of total states
        nFreeState      % number of free states
        nFreeTime       % number of free time durations
        nFreeCtrl       % number of free control parameters
        nFreeEps        % number of free epsilon parameter
        nTime           % number of domain time
        nCtrl           % number of control parameters (not including epsilon)
        options         % tolerances, gradients, ...
    end
    


    methods 
        %  - constructor - 
        function obj = RootFindingProblem(model,sequence,varargin) 
            if nargin < 2
                error('To instantiate RootFinding() the minimum amount of arguments must be model and sequence!')
            end
            
            % check if model is of type MechanicalModel
            if ~strcmp(model.type,'MechanicalModel')
                disp(pogoStick)
                error('model must be of type MechanicalModel!')
            else
                obj.model    = model;
            end
            
            if size(sequence.FlowMapOrder,1) ~= model.nZ && model.nZ>0
                error(['dimension of discrete states in sequence must be ',num2str(model.nZ),'!'])
            else
                obj.sequence = sequence;
            end
                
            obj.generalSetup(varargin{:});
        end
        
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    methods (Access = private)
        generalSetup(obj,varargin) % setup in constructor
        nDecVar = getNumDecVar(obj)
        function updateNumDecVar(obj)
            obj.nDecVar = obj.getNumDecVar();
        end 
        nCon    = getNumConstraints(obj)
        nInt    = getNumIntegrands(obj)
        nTime   = getNumDomainTime(obj)
    end
    
    methods (Access = public) % get and set functions
        [x0,tDomain,xi,epsilon,decVar] = getStateTimeCtrlEpsFromDecVar(obj,decVar,sol)
        Lxx = getHessian(obj,rootFun,x,lambda)
        [hxx,cxx] = get2ndOrderDerivatives(obj,rootfun,x,hx,cx)
        [funOpt,jacobianOpt] = getOptimJacobian2rootFunction(obj,decVar,fun,jacobian,objective,init)
        setOptimalControl(obj,optLogical)
        setFunctional(obj,idx,type,varargin)
        removeFunctional(obj,idx)
        removeVarFromOptimization(obj,type,idx)
        setSolverOptions(obj,options)
        fixParameter(obj,paramLogical)
        decVar  = getDecisionVariable(obj,Sol,varargin)
        [sol,output] = getOutputSolvers(obj,init,decVar,output,jacobian,exitflag)
        function trajectory = getTrajectory(obj,x0,tDomain)
            trajectory = obj.simulation.getTrajectory(x0,tDomain);
        end
    end
    
    methods (Access = public) % root functions
        % compute rootFunction by time-driven single shooting
        [fun,jacobian,objective] = rootFunctionTDsingle(obj,init,decVar,computeLxx)
    end

    methods (Access = public) % - Root-Finding-Methods
        % executable solvers and functions
        [sol,fval,exitflag,output,jacobian] = fsolveRFP(obj,init)
        [sol,fval,exitflag,output,jacobian] = solveNewtonRFP(obj,init)
    end

    methods (Access = private)
        % additional methods for executable solvers
        scalingMatrix = getScalingMatrix(obj,useLagrangeVar,varargin)
        options   = setOptionsfsolve(obj)
        [xE,intE] = getGradientsFromSimData(obj,simData)
    end

end

