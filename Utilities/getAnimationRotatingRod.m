function F_complete = getAnimationRotatingRod(data,varargin)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

    
% Colours [R,G,B] for torso;spring;foot
colours = [230, 230, 230; 255, 0, 0; 0, 0, 0]/255; 


rotatingRodGraphic = RotatingRodGraphicX(colours,varargin{:});

% --- create animation for this data case ---
F_complete = [];
for iDomain = 1:size(data,2)
    F(length(data(iDomain).t)) = struct('cdata',[],'colormap',[]); 
    frame_counter = 1;
    for i=1:length(data(iDomain).t)
        state_ = data(iDomain).x(i,:)';
        %z      = data(iDomain).z(i,:)';

        update(rotatingRodGraphic, state_');
        % append
        F(frame_counter) = getframe(gcf);
        frame_counter = frame_counter + 1;
    end
    F_complete = [F_complete, F];

    clearvars F % s.t. F can be resized
end
end

