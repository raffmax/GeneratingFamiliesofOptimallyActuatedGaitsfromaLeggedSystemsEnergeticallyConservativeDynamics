function removeVarFromOptimization(obj,type,idx)
%removeVarFromOptimization remove variable as part of necessary condition
%   the 1st order condition excludes this variable
if nargin<2
    obj.idxVarInEqCon = [obj.idx.FreeState,obj.idx.Time,obj.idx.Ctrl,obj.idx.Eps];
else
    switch type
        case 'x0'
            idxInFreeStates = find(obj.freeStates==idx);
            if isempty(idxInFreeStates)
                warning(['The state x_0(',num2str(idx),') is already a fixed state!'])
            else
                obj.idxVarInEqCon(idxInFreeStates) = [];
            end           
        case 'tDomain'
            idxTdomain = find(obj.idxVarInEqCon == obj.idx.FreeTime(1))-1+idx;
            obj.idxVarInEqCon(idxTdomain) = [];
        case 'xi'
            idxXi = find(obj.idxVarInEqCon == obj.idx.FreeCtrl(1))-1+idx;
            obj.idxVarInEqCon(idxXi) = [];      
        case 'epsilon'
            idxEps = find(obj.idxVarInEqCon == obj.idx.Eps(1))-1+idx;
            obj.idxVarInEqCon(idxEps) = [];      
        otherwise
            warning(['The type -',type,'- is not supported!'])    
    end
end

end

