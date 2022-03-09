% DEMO_MONOPED A demo on how to use the ModelHomotopy framework.
%   This script demonstrates how to systematically explore passive gaits
%   of an energetically conservative system (epsilon=0). 
%   It finds a forward hopping motion at xDot_avg = 1 and utilizes the
%   first order optimality to trace the gait into the desired dissipative
%   system at epsilon=1
clear, clc, close all;
install_project()

%% load model
% - sub-class of MechanicalModel -
hopper = HopperRemyHomotopy5DOF('load', false);

%% fixed contact sequence

% configure sequence:
% for this model, we have a constraint during flight and a constraint
% during stance, thus the phase variable is dim(z) = 2
% z = [stance constraint active;....
%      flight constraint active]
sequence.FlowMapOrder   = [1, 0;... % stance constraint
                           0, 1];   % flight constraint
% the jump map follows after the flowMap and needs the active
% constraint configuration of the next flowMap.
% Thus, JumpMapOrder is usually just a shifted FlowMapOrder.
sequence.JumpMapOrder   = [0, 1;...
                           1, 0];
% define sequence of events, that should be triggered at domain times t_i
% the specification (here 1,2) is user-defined and corresponds to the event
% implementation in their MechanicalModel
% 1: touch-down
% 2: lift-off
% 3: apax
sequence.events         = [2, 1];

% coordinates
% q = [x,y,phi,alpha,l]

% fixing the state to its inital value is done here.
fixedStates             = 1; % fixed at initial state -> A_np
% define the states which need to be periodic
periodicStates          = (2:10); % -> A_p

%% hopping-in-place

% get model parameters to compute mode durations
m   = hopper.parameters.param_m_t(0);
g   = hopper.parameters.g;
k   = hopper.parameters.param_k_l(0);
l0  = hopper.parameters.l_0;
r_f = hopper.parameters.param_r_f(0);

% set initial flight phase
tF  = 0.1;


y0  = l0+r_f;
dy0 = -g*tF/2;
% stance phase: x = a*sin( w*t + phi )
% m ddx = -k x
% x = y - l0 - r_f + mg/k
x0  = y0-l0-r_f+m*g/k;
dx0 = dy0;

w   = sqrt(k/m);
phi = atan(w*(x0/dx0));

tS  = 2*(pi/2-phi)/w;

x_init = [0; l0+r_f; 0;0; l0; 0; dy0; 0;0; dy0];
% domain times: stance, flight
t_init  = [tS, tF]';
% control parameter (numerical damping)
xi_init = 0;
% model homotopy parameter
epsilon = 0;
% store initial results in Solution
sol = Solution(x_init,t_init,xi_init,epsilon);
sol.param.E = hopper.SystemEnergy(x_init,sequence.FlowMapOrder(:,1));

% set anallytic first order gradients
opts.Grad1 = true; 

% set up RootSearch to define the search for a limit cycle
rfp = RootFindingProblem(hopper,sequence,'fixedStates',fixedStates,...
                                         'fixedEpsilon',1,...
                                         'periodicStates',periodicStates,...
                                         'simulation',Simulation('odeSolver',@ode113,'odeOptions',odeset('RelTol',1e-8,'AbsTol',1e-9)),...
                                         'options',opts);
% rfp.removeFunctional(1);
% rfp.setFunctional(1,'xDot_avg');                                    
                                    
% create limit cycle (linked list object)
solAn  = rfp.solveNewtonRFP(sol);
LC     = LimitCycle(['E=',num2str(sol.param.E)],solAn,rfp);

%% run Continuation in E
arcPC  = PseudoArclengthPC(rfp,LC,'fixedStepSize',true,'StepSize',1e-4);

% run continuation in energy level and terminate after 1000 interations
[exitflag,PC_steps] = arcPC.runContinuation('E');

% approximate all simple bifurcations
[LCBefore,LCAfter] = LC.findNextBifurcation();
[LCBif1, flag1]    = arcPC.approximateRoot(LCBefore,LCAfter, 'detJacAug'); 
[LCBefore,LCAfter] = LCBif1.findNextBifurcation();
[LCBif2, flag2]    = arcPC.approximateRoot(LCBefore,LCAfter, 'detJacAug'); 

% get new gaits with forward/backward motion
[LC1,LC2] = arcPC.branchingOff(LCBif1);
[LC3,LC4] = arcPC.branchingOff(LCBif2);

%% run Continuation in xDot_avg
% change parameterization of energetically conservative gaits to xDot_avg
rfp.removeFunctional(1);
rfp.setFunctional(1,'xDot_avg');
% select forward motion
sol=LC2.Sol;
sol.param=[];
trajectory = LC2.getTrajectory;
sol.param.xDot_avg=trajectory(end).x(end,1)/trajectory(end).t(end);
% create linked list object
M1 = LimitCycle(['xDot_avg=',num2str(sol.param.xDot_avg)],sol,rfp);
arcPC  = PseudoArclengthPC(rfp,M1,'fixedStepSize',true,'Stepsize',1e-2);

tic
% find forward hopping at xDot_avg=1
arcPC.runContinuation('xDot_avg',1);
toc

%% set up constrained optimization problems for xDot_avg=1
% controller Fourier series
degree     = 3;
xiOpt      = zeros(2*(2*degree+1),1);
controller = FourierSeriesController(hopper, xiOpt, degree);
                                                              

% set up RootSearch with First-Order Optimality
rfpOpt = RootFindingProblem(hopper,sequence,'fixedStates',fixedStates,...
                                         'periodicStates',periodicStates,...
                                         'simulation',Simulation('odeSolver',@ode113,'odeOptions',odeset('RelTol',1e-8,'AbsTol',1e-9)),...
                                         'controller',controller,...
                                         'useLagrangeMultiplier',true,...
                                         'optimalControl',true,...
                                         'options',opts);

% epsilon is not part of optimization
rfpOpt.removeVarFromOptimization('epsilon',1);
rfpOpt.removeFunctional(1);
rfpOpt.setFunctional(1,'xDot_avg');                                  
% add cost functional tau^2
rfpOpt.setFunctional(2,'tau^2');

% load solution - xDot_avg = 1
sol = M1.getLCend.Sol;
% store initial results (optimal solution) in Solution
solOpt = Solution(sol.x0,sol.tDomain,xiOpt,sol.epsilon,sol);
solOpt.rfpData = [];
solOpt.rfpData.multiplier = []; % is computed in solveNewtonRFP()
% compute Jacobian and Lagrange Multipliers
[solOpt,fvalOpt,exitflagOpt,optionsOpt,jacobianOpt] = rfpOpt.solveNewtonRFP(solOpt);

% create linked list object
LCopt  = LimitCycle('eps=0',solOpt,rfpOpt);
arcPC  = PseudoArclengthPC(rfpOpt,LCopt,'MaxIterPred',10,'Stepsize',0.1);

tic
% continuation on 1st Order Optimality from epsilon=0 to epsilon=1
[exitflag,PC_steps] = arcPC.runContinuation('epsilon',1);
toc

%% Plots
continuationEpsilon = LCopt.getLCstart.getContinuationData;
continuationPassive = M1.getLCstart.getContinuationData;

plotTestFunctions(continuationEpsilon,continuationEpsilon.epsilon)

figure
grid on
hold on
x1 = continuationPassive.xDot_avg';
x2 = continuationPassive.x0(:,2)';
x3 = continuationPassive.x0(:,4)';
plot3(x1,x2,x3,'r','LineWidth',2.0)
scatter3(x1(end),x2(end),x3(end),'filled','ro')

x1 = continuationEpsilon.xDot_avg';
x2 = continuationEpsilon.x0(:,2)';
x3 = continuationEpsilon.x0(:,4)';
plot3(x1,x2,x3,'k')
scatter3(x1(end),x2(end),x3(end),'filled','bo')

xlabel('$\dot{x}_\mathrm{avg}$ [$\sqrt{gl_o}$]','Interpreter','latex','FontSize',8)
ylabel('$y_0$ [$l_o$]','Interpreter','latex','FontSize',8)
zlabel('$\alpha_0$ [rad]','Interpreter','latex','FontSize',8)
view(-30,30)

%% Animation
% gifs of optimal forward motions (xDot_avg=1) for epsilon=0 and epsilon=1
dt = 3e-2;
frameEps0 = getAnimationHopper5DOF(LCopt.getLCstart.getTrajectory(dt),0,false);
frameEps1 = getAnimationHopper5DOF(LCopt.getLCend.getTrajectory(dt),1,false);

% export as gif
movie2gif(frameEps0, 'MonopedPassive.gif','DelayTime',.05,'LoopCount',Inf)
movie2gif(frameEps1, 'MonopedActive.gif','DelayTime',.05,'LoopCount',Inf)
