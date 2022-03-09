function getAnimationCompassGait(dataIN,motion,epsilon)
%getAnimationCompassGait Summary of this function goes here
%   Detailed explanation goes here

t_ = [];
y_ = [];
for i = 1:size(dataIN,2)
    t_ = [t_;dataIN(i).t];
    y_ = [y_;dataIN(i).x];
end
q0 = y_(1:4)';
gamma = epsilon;

animate_FloatingBase(t_,y_,q0,gamma,motion(:,1))

end

