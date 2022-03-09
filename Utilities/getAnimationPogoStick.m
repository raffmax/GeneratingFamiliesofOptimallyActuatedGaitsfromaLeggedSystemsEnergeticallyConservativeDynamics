function F_complete = getAnimationPogoStick(dataIN,epsilon,full3D)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

if full3D % change the 3D lighting effects
    edges         = false; % show edges for torso and legs
else
    edges         = true; % show edges for torso and legs
end
    
% Colours [R,G,B] for torso;spring;foot
colours = [230, 230, 230; 255, 0, 0; 0, 0, 0]/255; 

% check if we simulate 1D or 2D
if size(dataIN(1).x,2) == 4
   % add states
   % stateIn  = [y,l]
   % stateOut = [x,y,phi,l]
   data = dataIN;
   for iDomain = 1:size(data,2)
       n    = length(data(iDomain).t);
       col0 = zeros(n,1);
       data(iDomain).x  = [col0,dataIN(iDomain).x(:,1),col0,dataIN(iDomain).x(:,2),...
                           col0,dataIN(iDomain).x(:,3),col0,dataIN(iDomain).x(:,4)];
   end
else
    data = dataIN;
end


pogoStick_graphic = PogoStickgraphicX(full3D,edges,colours,epsilon);

% --- create animation for this data case ---
F_complete = [];
for iDomain = 1:size(data,2)
    F(length(data(iDomain).t)) = struct('cdata',[],'colormap',[]); 
    frame_counter = 1;
    for i=1:length(data(iDomain).t)
        state_ = data(iDomain).x(i,:)';

        update(pogoStick_graphic, state_', epsilon);
        % append
        F(frame_counter) = getframe(gcf);
        frame_counter = frame_counter + 1;
    end
    F_complete = [F_complete, F];

    clearvars F % s.t. F can be resized
end
end

