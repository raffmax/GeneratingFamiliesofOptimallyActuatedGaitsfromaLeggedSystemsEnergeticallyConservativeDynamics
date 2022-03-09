function createSimpleAnimationContinuationPlanarBiped(LC,fixedParam,dt,type,full3D)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

if full3D % change the 3D lighting effects
    edges         = false; % show edges for torso and legs
else
    edges         = true; % show edges for torso and legs
end

COG = true;
% Colours [R,G,B] for torso;legs;springs
colours = [230, 230, 230; 160, 160, 160; 255, 0, 0]/255; 


epsilon = 0;
Biped_graphic = PlanarBipedV1graphicX(full3D,COG,edges,colours,epsilon,epsilon,true);
fig = gcf; % get figure
fig.Position =  [1 1 1000 1000];

% check screen size
screensize = get( groot, 'Screensize' );
if min(screensize(3:4)) < 1000
    warning('Your screen size is too small to cover the required 1000x1000 pixels!') 
    fig.Position =  [1 1 700 700];
end

nDomains   = LC.nTime; % number of domains
if ~isequal(LC.flipOperator,eye(LC.nState))
    nDomains = 2*nDomains;
end
nCon       = LC.treeDepth+1; % number of continuation steps
F_complete = []; % figure data

for iLC = 1:nCon
    if iLC >LC.treeDepth/2
    
    data = LC.getTrajectory(dt);
    if mod(iLC,2) == 1
        rangeDomain = (1:nDomains/2);
    else
        rangeDomain = (nDomains/2+1:nDomains);
    end
    for iDomain = rangeDomain
        z = data(iDomain).z;
        t = data(iDomain).t;
        x = data(iDomain).x;
    
        switch type
            case  'epsilon'
                epsilonArray   = LC.Sol.epsilon*ones(1,length(t));
                amplitudeArray = zeros(1,length(t));
            otherwise
                amplitudeArray   = LC.Sol.param.(type)*ones(1,length(t));
                epsilonArray = fixedParam.value*ones(1,length(t));
        end
    
    
        F(length(t)) = struct('cdata',[],'colormap',[]); 
        frame_counter = 1;

        for i=1:length(t)
            state_  = x(i,:);
            z_      = z(i,:)';
            epsilon = epsilonArray(i);
            a       = amplitudeArray(i);

            % set actuation 
            u = [0;0;0;0];

            update(Biped_graphic, state_', u, z_,epsilon,epsilon);

            switch type
                case 'epsilon'
                    txtEps   = text(state_(1)-1.2,0,2.2,['$\varepsilon= ',num2str(round(epsilon,2)),'$'],'FontSize',36,'Color','blue','Interpreter','latex');
                    txtParam = text(state_(1)+.15,0,2.2,['$',fixedParam.symbol,' = ',num2str(round(fixedParam.value,2)),'~',fixedParam.unit,'$'],'FontSize',36,'Color','black','Interpreter','latex');
                    stretchRectanglePlot = false;
                otherwise
                    switch type
                        case 'E0'
                            txtParam1 = text(state_(1)-1.2,0,2.2,['$E_0 = ',num2str(round(a,2)),'$'],'FontSize',36,'Color','blue','Interpreter','latex');
                            txtParam2 = text(state_(1)-0.55,0,2.2,'$m g l_0$','FontSize',36,'Color','blue','Interpreter','latex');
                        case 'E_avg'
                            txtParam1 = text(state_(1)-1.2,0,2.2,['$E_\mathrm{avg} = ',num2str(round(a,2)),'$'],'FontSize',36,'Color','blue','Interpreter','latex');
                            txtParam2 = text(state_(1)-0.55,0,2.2,'$m g l_0$','FontSize',36,'Color','blue','Interpreter','latex');
                        case 'x6'
                            txtParam1 = text(state_(1)-1.2,0,2.2,['$\dot x_0 = ',num2str(round(a,2)),'$'],'FontSize',36,'Color','blue','Interpreter','latex');
                            txtParam2 = text(state_(1)-0.55,0,2.2,'$\sqrt{g l_0}$','FontSize',36,'Color','blue','Interpreter','latex');
                    end
                    txtEps    = text(state_(1)+.15,0,2.2,['$\varepsilon= ',num2str(round(epsilon,2)),'$'],'FontSize',36,'Color','black','Interpreter','latex');
                    stretchRectanglePlot = true;
            end
            txtT   = text(state_(1)+.15,0,2.0,['$T~= ',sprintf('%.2f',data(end).t(end)),'~\sqrt{l_o/g}$'],'FontSize',36,'Color','black','Interpreter','latex');
            % plot rectangle
            x1 = state_(1)-1.25;
            if ~stretchRectanglePlot
                x2 = state_(1)-0.65;
            else
                x2 = state_(1)-0.25;
            end
            x3 = x2;
            x4 = x1;
            y  = 0;
            z1 = 2.3;
            z2 = 2.3;
            z3 = 2.1;
            z4 = 2.1;
            hold on
            rectangle = plot3( [x1 x2 x3 x4 x1], [y y y y y], [z1 z2 z3 z4 z1],'b');

            % append
            F(frame_counter) = getframe(fig);
            if stretchRectanglePlot
                delete(txtEps)
                delete(txtT)
                delete(txtParam1)
                delete(txtParam2)
                delete(rectangle)
            else
                delete(txtEps)
                delete(txtT)
                delete(txtParam)
                delete(rectangle)
            end
            frame_counter = frame_counter + 1;
        end
        F_complete = [F_complete, F];
        clearvars F % s.t. F can be resized
    end
    end
    LC = LC.Next;
end

video_gif_ID = '\planarBiped';
v = VideoWriter(strcat(pwd, video_gif_ID, '.mp4'),'MPEG-4');
open(v)
writeVideo(v,F_complete)
close(v)
clearvars v 
close gcf
end