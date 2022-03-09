function F_complete = getAnimationHopperPitch(data,epsilon,full3D)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

 if full3D % change the 3D lighting effects
    edges         = false; % show edges for torso and legs
else
    edges         = true; % show edges for torso and legs
end
    
% Colours [R,G,B] for torso;spring;foot
colours = [230, 230, 230; 255, 0, 0; 0, 0, 0]/255; 

% stateIn  = [x,y,phi,alpha,l]

hopper_graphic = HoppergraphicX(full3D,edges,colours,epsilon);

F_complete = [];
for iDomain = 1:size(data,2)
    F(length(data(iDomain).t)) = struct('cdata',[],'colormap',[]); 
    frame_counter = 1;
    for i=1:length(data(iDomain).t)
        state_ = data(iDomain).x(i,:)';

        update(hopper_graphic, state_', epsilon);
        % append
        F(frame_counter) = getframe(gcf);
        frame_counter = frame_counter + 1;
    end
    F_complete = [F_complete, F];

    clearvars F % s.t. F can be resized
end
end

