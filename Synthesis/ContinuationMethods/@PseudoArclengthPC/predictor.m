function [uP,H_uP,Hprime_uP,tP,exitflagP,h,aimOnTarget] = predictor(obj,d,lambda_target,u,H_u,t_u)
    % Allgower[2003] chapter 6.2
    % returns u_Predictor

    % set default step size in Predictor
    h = obj.Direction*obj.StepSize;
    
    % try to get as close to the target homotopy parameter as possible
    aimOnTarget = false;
    
    f       = 2; % deceleration factor
    iter    = 0;
    
    h_trials    = zeros(obj.MaxIterPred,1); % trials for step lengths
    f_trials    = zeros(obj.MaxIterPred,1); % trials for step lengths
    
    v_trials        = cell(obj.MaxIterPred,1);
    H_v_trials      = cell(obj.MaxIterPred,1);
    Hprime_v_trials = cell(obj.MaxIterPred,1);
    t_v_trials      = cell(obj.MaxIterPred,1);
    

    while (f==0.5 || f==2) && iter < obj.MaxIterPred
        iter = iter + 1;
        if iter == 1
            f = 1;
        end
        
        if aimOnTarget
            h = (lambda_target-u(obj.idxConPar))/t_u(obj.idxConPar);
        else
            h = h/f; % step length adaptation
            h_trials(iter) = h;
        end
        % Newton step
        v   = u + h*t_u;
        
        lambdaP        = v(obj.idxConPar); 
        [H_v,Hprime_v] = obj.getH(v); % compute value and jacobian of the root function rootFunctionTDsingle();
        t_v            = getTangent(Hprime_v);
        
        if aimOnTarget
            % end corrector successfully
            f = 1;
            break
        end
        
        % determine variable step size
        if ~obj.fixedStepSize
            % pseudo-inverse of Hprime_v
            HprimePlus = pinv(Hprime_v);

            % first corrector steplength
            delta = norm(HprimePlus*H_v);
            % contraction rate
            kappa = norm(HprimePlus*H_u)/delta;
            if obj.RFP.optimalControl && obj.RFP.useLagrangeMultipliers  
                t_orient=1;
            else
                % angle between two consecutive steps
                if t_u'*t_v<0 % simple bifuraction
                    t_orient = -1;
                else
                    t_orient = 1;
                end
            end
            alpha = acos(t_orient*t_u'*t_v);
            
            % save trial step
            v_trials{iter}        = v;
            H_v_trials{iter}      = H_v;
            Hprime_v_trials{iter} = Hprime_v;
            t_v_trials{iter}      = t_v;

            % get new h
            f = max([sqrt(kappa/obj.kappaTilde),...
                     sqrt(delta/obj.deltaTilde),...
                     alpha/obj.alphaTilde]);
            f_trials(iter) = f;
            f = max(min(f,2),.5);
            
            if iter>2 && f == f_trials(iter-2)
                if f == 2
                    f = f*0.75;
                else
                    f = f*1.5;
                end
            end
            
            if max(abs(H_v)) <= obj.RFP.options.FunctionTolerance
                f=1;
            end
            
        end
        if ~isempty(lambda_target)
            if f < 2 && d*lambdaP >= d*lambda_target
                f = 2;
                aimOnTarget = true;
            end
        end
    end
    if (f>0.5 && f<2)
        exitflagP = 1;
        uP        = v;
        H_uP      = H_v;
        Hprime_uP = Hprime_v;
        tP        = t_v;
    else
        exitflagP = -1;
        [~,idx]   = min(abs(f_trials)-1.25);
        uP        = v_trials{idx};
        H_uP      = H_v_trials{idx};
        Hprime_uP = Hprime_v_trials{idx};
        tP        = t_v_trials{idx};
        warning('The predictor-step was not successful!')
    end
end