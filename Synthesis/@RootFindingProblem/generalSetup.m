function generalSetup(obj, varargin)
    % default values for the configuration
    default_fixedStates      = [];
    default_fixedTimes       = [];
    default_fixedCtrls       = [];
    default_fixedEpsilon     = [];
    default_periodicStates   = (1:2*obj.model.nQ);
    default_periodicStatesHomotopy = [];
    default_fixedParameter   = true;
    default_optimalControl   = false;
    default_useLagrangeMultipliers = true;
    default_flipOperator     = eye(2*obj.model.nQ);
    default_controller       = ConservativeController(obj.model, 0);
    default_simulation       = Simulation();
%     default_scaleControl     = 1;
%     default_scaleLagrange    = 1;

    default_integrationType  = 'TD';
    default_shootingType     = 'single';
    default_options          = []; % will be set after parser


    p = inputParser;
    % cfg 
    addParameter(p, 'fixedStates',    default_fixedStates);
    addParameter(p, 'fixedTimes',     default_fixedTimes);
    addParameter(p, 'fixedCtrls',     default_fixedCtrls);
    addParameter(p, 'fixedEpsilon',   default_fixedEpsilon);

    addParameter(p, 'periodicStates',  default_periodicStates);
    addParameter(p, 'periodicStatesHomotopy',  default_periodicStatesHomotopy);
    addParameter(p, 'flipOperator',    default_flipOperator);
    addParameter(p, 'fixedParameter',  default_fixedParameter);
    addParameter(p, 'optimalControl',  default_optimalControl);
    addParameter(p, 'useLagrangeMultipliers',  default_useLagrangeMultipliers);
    addParameter(p, 'controller',      default_controller);
    addParameter(p, 'simulation',      default_simulation);
%     addParameter(p, 'scaleControl',    default_scaleControl);
%     addParameter(p, 'scaleLagrange',   default_scaleLagrange);

    addParameter(p, 'integrationType', default_integrationType);
    addParameter(p, 'shootingType',    default_shootingType);
    addParameter(p, 'options',         default_options);

    parse(p, varargin{:}); 

    obj.fixedStates        = p.Results.fixedStates;
    obj.fixedTimes         = p.Results.fixedTimes;
    obj.fixedCtrls         = p.Results.fixedCtrls;
    obj.fixedEpsilon       = p.Results.fixedEpsilon;
    obj.periodicStates     = p.Results.periodicStates;
    obj.periodicStatesHomotopy = p.Results.periodicStatesHomotopy;
    obj.flipOperator       = p.Results.flipOperator;
    obj.fixedParameter     = p.Results.fixedParameter;
    obj.optimalControl     = p.Results.optimalControl;
    obj.useLagrangeMultipliers = p.Results.useLagrangeMultipliers;
    obj.controller         = p.Results.controller;
    obj.simulation         = p.Results.simulation;
%     obj.scaleControl       = p.Results.scaleControl;
%     obj.scaleLagrange      = p.Results.scaleLagrange;
    obj.simulation.setSystemRef(obj); % add reference
    obj.integrationType    = p.Results.integrationType;
    obj.shootingType       = p.Results.shootingType;
    
    obj.nState     = 2*obj.model.nQ;
    obj.nTime      = obj.getNumDomainTime();
    obj.nCtrl      = obj.controller.nXi;
    obj.freeStates = default_periodicStates(~ismember(default_periodicStates,obj.fixedStates));
    idxTimes       = 1:obj.nTime;
    obj.freeTimes  = idxTimes(~ismember(idxTimes,obj.fixedTimes));
    idxCtrls       = 1:obj.nCtrl;
    obj.freeCtrls  = idxCtrls(~ismember(idxCtrls,obj.fixedCtrls));
    if isempty(obj.fixedEpsilon)
        obj.freeEps = 1;
    end
    obj.nFreeState = length(obj.freeStates);
    obj.nFreeTime  = length(obj.freeTimes);
    obj.nFreeCtrl  = length(obj.freeCtrls);
    obj.nFreeEps   = length(obj.freeEps);
    obj.nCon       = obj.getNumConstraints();
    obj.nInt       = obj.getNumIntegrands();
    obj.nDecVar    = obj.getNumDecVar();
    
    % set solver options
    obj.setSolverOptions(p.Results.options)
    
    % set default structure for functionals
    if obj.fixedParameter
        obj.setFunctional(1,'E')
    end
end