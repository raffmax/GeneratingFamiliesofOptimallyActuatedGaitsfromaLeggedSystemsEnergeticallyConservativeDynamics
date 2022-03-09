classdef HopperParametersRemyHomotopy5DOF < handle
    % HopperParametersRemyHomotopy5DOF contains all the mechanical parameters for the 1D pogo stick.
    %   This class exists in order to not clutter the functionality in
    %   actual pogo stick class.
    %
    %   Input:  
    %      coefRes       coefficient of restitution
    %
    %   Methods:
    %      PogoStickParameters   The constructor for this class
    %
    %   Example:
    %      parameters = PogoStickParameters2D(coefRes);
    %
    %   See also PogoStick1D, MechanicalModel
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        University of Stuttgart, Institute for Nonlinear Mechanics
    % Author:        Maximilian Raff
    % email address: raff@inm.uni-stuttgart.de
    % Website:       https://www.inm.uni-stuttgart.de/en
    % Last revision: 12-July-2021
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        coefRes            % coefficient of restitution                                                  
         
        % - pogo stick parameters -
        m       % total mass
        m_t     % mass of torso / upper body relative to m
        m_l     % mass of leg
        m_f     % mass of foot  / lower body relative to m
        theta_l % leg inertia
        theta_f % foot inertia
        theta_t % torso inertia 
        r_f     % foot radius
        g       % gravity
        l_0     % natural length of spring
        d_l     % distance hip to COG leg
        d_f     % distance foot center to COG foot
        k_l     % leg spring stiffness
        b_l     % leg spring damping   
        k_h     % hip spring stiffness
        b_h     % hip spring damping        
    end
    
%     properties (Constant)
%         g                           = 1                                         % normalized gravity
%         l_0                         = 1                                         % normalized leglenght
%     end
    
    methods
        % - constructor - 
        function this = HopperParametersRemyHomotopy5DOF(coefRes,epsilon)
            % set mass and actuation scaling values
            this.coefRes = coefRes;
            if nargin <2
                epsilon = sym('epsilon','real');
            end
            
            % physical values
            this.m          = this.param_m(epsilon);  % total mass
            this.m_f        = this.param_m_f(epsilon);
            this.m_l        = this.param_m_l(epsilon);
            this.m_t        = this.param_m_t(epsilon);
            this.g          = this.param_g(epsilon);
            this.l_0        = this.param_l_0(epsilon);
            this.r_f        = this.param_r_f(epsilon);
            this.theta_l    = this.param_theta_l(epsilon);
            this.theta_f    = this.param_theta_f(epsilon);
            this.theta_t    = this.param_theta_t(epsilon);
            this.d_l        = this.param_d_l(epsilon);
            this.d_f        = this.param_d_f(epsilon);
            this.k_l        = this.param_k_l(epsilon);
            this.b_l        = this.param_b_l(epsilon);
            this.k_h        = this.param_k_h(epsilon);
            this.b_h        = this.param_b_h(epsilon);
        end
        function m = param_m(~,~)
            m = 1;
        end
        function m_t = param_m_t(this,~)
            m_t = .7*this.m;
        end
        function m_l = param_m_l(this,epsilon)
            m_l = .2*this.m*epsilon;
        end
        function m_f = param_m_f(this,epsilon)
            m_f = .1*this.m*epsilon;
        end
        function g = param_g(~,~)
            g = 1;
        end
        function l_0 = param_l_0(~,~)
            l_0 = 1;
        end
        function theta_l = param_theta_l(this,epsilon)
            l_0 = this.l_0;
            r   = sqrt(2)/10*l_0; % radius of gyration
            m_l = this.param_m_l(epsilon);
            theta_l = m_l*r^2;
        end
        function theta_f = param_theta_f(this,epsilon)
            l_0 = this.l_0;
            r   = 0.2*l_0; % radius of gyration
            m_f = this.param_m_f(epsilon);
            theta_f = m_f*r^2;
        end
        function theta_t = param_theta_t(this,epsilon)
            l_0 = this.l_0;
            r   = sqrt(40/70)*l_0; % radius of gyration
            m_t = this.param_m_t(epsilon);
            theta_t = epsilon*m_t*r^2;
        end
        function r_f = param_r_f(this,~)
            r_f = 0.05*this.l_0;
        end
        function d_l = param_d_l(this,~)
            d_l = 0.25*this.l_0;
        end
        function d_f = param_d_f(this,~)
            d_f = 0.25*this.l_0;
        end
        function k_l = param_k_l(this,~)
            m   = this.m;
            l_0 = this.l_0;
            g   = this.g;
            k_l = 20*m*g/l_0;
        end
        function dk_l = grad_k_l(~, ~)
            dk_l = 0;
        end
        function b_l = param_b_l(this,epsilon)
            k_l = this.param_k_l(1);
            m_f = this.param_m_f(1);
            b_l = 0.2*2*sqrt(k_l*m_f)*epsilon; % 0.4*sqrt(2)*epsilon
        end
        function db_l = grad_b_l(this,~)
            k_l  = this.param_k_l(1);
            m_f = this.param_m_f(1);
            db_l = 0.2*2*sqrt(k_l*m_f);
        end
        function k_h = param_k_h(this,epsilon)
            m_f  = this.param_m_f(epsilon);
            g    = this.g;
            l_0  = this.l_0;
            k_h  = 100*m_f*g*l_0; % k_h = 5*m_f*g*l_0
        end
        function b_h = param_b_h(this,epsilon)
%             l_0     = this.l_0;
%             k_h     = this.param_k_h(1);
%             theta_f = this.param_theta_f(1);
%             theta_l = this.param_theta_l(1);
%             m_f     = this.param_m_f(1);
%             m_l     = this.param_m_l(1);
%             d_f     = this.param_d_f(1);
%             d_l     = this.param_d_l(1);
%             
%             theta_legHip = theta_f + (l_0 - d_f)^2*m_f + theta_l + d_l^2*m_l;
%             b_h = 0.2*2*sqrt(k_h*theta_legHip)*epsilon^2;
            b_h = 0.35*epsilon^2;
        end
    end
end

