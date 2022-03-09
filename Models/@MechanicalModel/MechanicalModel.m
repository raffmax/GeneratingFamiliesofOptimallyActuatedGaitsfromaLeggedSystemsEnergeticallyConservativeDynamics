classdef (Abstract) MechanicalModel < matlab.mixin.Copyable
    % MechanicalModel is a blueprint for a mechanical model with
    % a model homotopy.
    %   This abstract model provides the important attributes and methods
    %   for a mechanical model homotopy.
    %
    %   Properties:
    %      epsilon           model-homotopy parameter
    %      nEps              number of homotopy parameters
    %      parameters        model parameters
    %      identifier        model ID
    %      folder_path       directory for autogenerated data
    %      load_path         name/directory of object mat-file
    %
    %      M                 mass matrix
    %      h                 force collection (potential, centrifugal, coriolis, ...)
    %      B                 actutation matrix (projects active torques)
    %      Wdyn              contact jacobian in dynamics equation
    %      W                 contact jacobian in constraints equation
    %      Wdot              time derivative of contact jacobian
    %      T                 kinetic energy
    %      V                 potential energy
    %
    %      Grad1             struct with first order gradients
    %      type              identifier for sub-classes of this type
    %
    %
    %   Methods:
    %      PogoStick1D           The constructor for this class
    %      CreateHybridDynamics  create and save autogenerated functions (calls CreateEOM)
    %      setEpsilon            set new homotopy parameter
    %      getEpsilon            get current homotopy parameter
    %      SystemEnergy          get energy at current state x and active constraints z
    %
    %      CreateEOM             (abstract) autogenerate EOM
    %      getEventFunctional    (abstract) % select pre-defined event functionals
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 12-Apr-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        epsilon    % model-homotopy parameter
        parameters % general model parameters (mass, spring stiffnesses ...)
        
        %% management
        %
        identifier    % model ID
        folder_path   % directory for autogenerated data
        load_path     % name/directory of object mat-file
        
        %% EOM
        %
        M             % mass matrix
        h             % force collection (potential, centrifugal, coriolis, ...)
        B             % actutation matrix (projects active torques)
        Wdyn          % contact jacobian in dynamics equation
        
        %% constraints
        %
        W             % contact jacobian in constraints equation
        Wdot          % time derivative of contact jacobian
        
        %% energy
        %
        T             % kinetic energy
        V             % potential energy
        
        nEps          % number of homotopy parameters
    end
    
    properties % - Gradients -
        % struct with first order gradients
        Grad1
    end
    
    properties (Constant)
        % identifier for sub-classes of this type
        type = 'MechanicalModel'
    end
    
    properties (Abstract, Constant)
        nZ            % number of (physical) constraints
        % 0 or 1 if there are no active constraints
        nQ            % number of generalized coordinates
        nTau          % number of actuators
    end
    
    
    methods (Abstract) %% - must be implemented by the subclasses -
        % --- EOM creation ---
        CreateEOM(obj)  % computes function handles for the FlowMap and JumpMap
        [conEvent,eventData] = getEventFunctional(obj,t,x,z,T,iEvent,controller, Grad1, Grad2)
    end
    
    methods
        %% --- Energy ---
        function E = SystemEnergy(obj, x, z)
            % get energy at current state x and active constraints z
            q   = x(1:obj.nQ);
            dq  = x(obj.nQ+(1:obj.nQ));
            E   = obj.T(obj.epsilon, q, dq, z) + obj.V(obj.epsilon, q, dq, z);
        end
        
        %% --- Flow Map ---
        [f, lambda, fGrad, lambdaGrad, tau] = FlowMap(obj, t, x, z, T, controller, Grad1, Grad2, epsFix)
        %% --- Jump Map ---
        [xPlus, Lambda, jumpGrad]           = JumpMap(obj, xMinus, z_prev, z_next, Grad1, Grad2, epsFix)
    end
    
    methods %% --- Management ---
        function CreateHybridDynamics(obj)
            % create and save autogenerated functions
            if ~exist(obj.folder_path, 'dir')
                mkdir(obj.folder_path)
            end
            addpath(obj.folder_path);
            
            % - compute function handels for the EOM (Flow Map and Jump Map) -
            CreateEOM(obj);
        end
        
        function setEpsilon(obj,epsilon)
            % set new homotopy parameter
            obj.epsilon = epsilon;
        end
        function epsilon = getEpsilon(obj)
            % get current homotopy parameter
            epsilon = obj.epsilon;
        end
    end
    
    methods %% --- generation of EoM ---
        function [M, h] = EulerLagrange(~, T, V, q, dq)
            % Lagrangian
            L = simplify(T-V);
            
            dL_dq   = jacobian(L, q).';
            dL_dqdt = jacobian(L, dq).';
            
            % d/dt (dL_dqdt)
            dd_L_dqdt2  = jacobian(dL_dqdt, dq);
            d_dLdqdt_dq = jacobian(dL_dqdt, q);
            
            % assign to notation
            M = simplify(dd_L_dqdt2);
            h = dL_dq - d_dLdqdt_dq*dq.';
        end
        
        function [W, W_dot] = computeContactMatrix(~, con, q, dq)
            % compute projection matrix W and its derivative W_dot of a
            % constraint con
            
            W     = jacobian(con, q)';
            W_dot = W;
            
            % compute for each column of W the jacobian
            for j=1:size(W,2)
                W_dot(:,j) = jacobian(W(:,j), q)*dq'; % + 0 from d/dqdt
            end
        end
        
        function S = computeDerivativesEOM(~, S, var)
            % S - cell matrix with symbolic functions which will be derived
            % wrt to var
            
            % computation
            for i=1:length(S)
                if size(S{i,2},2)>1 % e.g. W
                    S{i,2} = cellfun(@(sfun) jacobian(sfun, var), S{i,2}, 'UniformOutput', false);
                    S{i,2} = cellfun(@simplify,                   S{i,2}, 'UniformOutput', false);
                else
                    S{i,2} = jacobian(S{i,2}, var); % differentiate wrt var
                    S{i,2} = simplify(S{i,2});      % simplify derivative
                end
            end
        end
        
        
        function generateMatlabHandels(obj, S)
            % S: cell matrix with symbolic terms
            
            % t: timestemp for unique matlab file names
            t = string(datetime('now', 'Format', 'yyMMdd_HHmm_')); % add timestamp to name of fcn

            
            file = fullfile(obj.folder_path, strcat('F_', t)); % 'F_date_time_
            epsilon = S{1}{1};
            q       = S{1}{2};
            dq      = S{1}{3};
            x       = S{1}{4};
            z       = S{1}{5};
            
            obj.T  = matlabFunction(S{2}, 'Optimize', true, 'file', strcat(file,'_T'), 'vars', {epsilon, q, dq, z});
            obj.V  = matlabFunction(S{3}, 'Optimize', true, 'file', strcat(file,'_V'), 'vars', {epsilon, q, dq, z});
            % - EOM -
            obj.M  = matlabFunction(S{4}{1}, 'Optimize', true, 'file', strcat(file,'_M'), 'vars', {epsilon, q});
            obj.h  = matlabFunction(S{5},    'Optimize', true, 'file', strcat(file,'_h'), 'vars', {epsilon, q, dq, z});
            obj.B  = matlabFunction(S{6}{1}, 'Optimize', true, 'file', strcat(file,'_B'), 'vars', {epsilon, q, dq, z});    % not every input is needed
            
            obj.Wdyn = [];
            obj.W    = [];
            obj.Wdot = [];
            for i=1:length(S{7})
                obj.Wdyn{i} = matlabFunction(S{7}{i}{1}, 'Optimize', true, 'file', strcat(file,'_Wdyn',num2str(i)), 'vars', {epsilon, q, z});
                obj.W{i}    = matlabFunction(S{8}{i}{1}, 'Optimize', true, 'file', strcat(file,'_W',   num2str(i)), 'vars', {epsilon, q, z});
                obj.Wdot{i} = matlabFunction(S{9}{i}{1}, 'Optimize', true, 'file', strcat(file,'_Wdot',num2str(i)), 'vars', {epsilon, q, dq, z});
            end

            % - Derivatives -
            E_x             = jacobian(S{2}+S{3},x);
            obj.Grad1.E_x   = matlabFunction(E_x,    'Optimize', true, 'file', strcat(file,'_E_x'),  'vars', {epsilon, x, z});
            E_xx            = jacobian(E_x,x);
            obj.Grad1.E_xx  = matlabFunction(E_xx,   'Optimize', true, 'file', strcat(file,'_E_xx'), 'vars', {epsilon, x, z});
            E_eps           = jacobian(S{2}+S{3},epsilon);
            obj.Grad1.E_eps = matlabFunction(E_eps,  'Optimize', true, 'file', strcat(file,'_E_eps'),'vars', {epsilon, x, z});    
            M_x             = simplify(jacobian(S{4}{1}*S{4}{2}, x));
            obj.Grad1.M_x   = matlabFunction(M_x,    'Optimize', true, 'file', strcat(file,'_M_x'),  'vars', {epsilon, x, S{4}{2}});
            M_eps           = simplify(jacobian(S{4}{1}*S{4}{2}, epsilon));
            obj.Grad1.M_eps = matlabFunction(M_eps,  'Optimize', true, 'file', strcat(file,'_M_eps'),'vars', {epsilon, x, S{4}{2}});
            h_x             = simplify(jacobian(S{5},x));
            obj.Grad1.h_x   = matlabFunction(h_x,    'Optimize', true, 'file', strcat(file,'_h_x'),  'vars', {epsilon, x, z});
            h_eps           = simplify(jacobian(S{5,1}, epsilon));
            obj.Grad1.h_eps = matlabFunction(h_eps,  'Optimize', true, 'file', strcat(file,'_h_eps'),'vars', {epsilon, x, z});
            B_xtau          = simplify(jacobian(S{6}{1}*S{6}{2}, x));
            obj.Grad1.B_x   = matlabFunction(B_xtau, 'Optimize', true, 'file', strcat(file,'_B_xtau'), 'vars', {epsilon, x, z, S{6}{2}});
            B_epstau        = simplify(jacobian(S{6}{1}*S{6}{2}, epsilon));
            obj.Grad1.B_eps = matlabFunction(B_epstau,'Optimize', true, 'file', strcat(file,'_B_epstau'), 'vars', {epsilon, x, z, S{6}{2}});
            
            obj.Grad1.Wdyn_x     = [];
            obj.Grad1.Wdyn_eps   = [];
            obj.Grad1.W_x_T      = [];
            obj.Grad1.W_eps_T    = [];
            obj.Grad1.Wdot_x_T   = [];
            obj.Grad1.Wdot_eps_T = [];
            for i=1:length(S{7})
                Wdyn_xlambda            = simplify(jacobian(S{7}{i}{1}*S{7}{i}{2}, x));
                obj.Grad1.Wdyn_x{i}     = matlabFunction(Wdyn_xlambda, 'Optimize', true, 'file', strcat(file,'_Wdyn_xlambda',num2str(i)), 'vars', {epsilon, x, z, S{7}{i}{2}});
                Wdyn_epslambda          = simplify(jacobian(S{7}{i}{1}*S{7}{i}{2}, epsilon) );
                obj.Grad1.Wdyn_eps{i}   = matlabFunction(Wdyn_epslambda,'Optimize', true, 'file', strcat(file,'_Wdyn_epslambda',num2str(i)), 'vars', {epsilon, x, z, S{7}{i}{2}});
                W_x_T                   = simplify(jacobian(S{8}{i}{1}.'*S{8}{i}{2}, x));
                obj.Grad1.W_x_T{i}      = matlabFunction(W_x_T       , 'Optimize', true, 'file', strcat(file,'_W_x_T',num2str(i)), 'vars', {epsilon, x, z, S{8}{i}{2}});
                W_eps_T                 = simplify(jacobian(S{8}{i}{1}.'*S{8}{i}{2}, epsilon));
                obj.Grad1.W_eps_T{i}    = matlabFunction(W_eps_T,     'Optimize', true, 'file', strcat(file,'_W_eps_T',num2str(i)),    'vars', {epsilon, x, z, S{8}{i}{2}});
                Wdot_x_T                = simplify(jacobian(S{9}{i}{1}'*S{9}{i}{2}, x));
                obj.Grad1.Wdot_x_T{i}   = matlabFunction(Wdot_x_T, 'Optimize', true, 'file', strcat(file,'_Wdot_x_T',num2str(i)), 'vars', {epsilon, x, z, S{9}{i}{2}});
                Wdot_eps_T              = simplify(jacobian(S{9}{i}{1}.'*S{9}{i}{2}, epsilon));
                obj.Grad1.Wdot_eps_T{i} = matlabFunction(Wdot_eps_T,  'Optimize', true, 'file', strcat(file,'_Wdot_eps_T',num2str(i)), 'vars', {epsilon, x, z, S{9}{i}{2}});
            end
        end
        
    end
end

