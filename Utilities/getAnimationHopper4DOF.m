function F_complete = getAnimationHopper4DOF(data,epsilon,full3D)
%UNTITLED21 Summary of this function goes here
%   Detailed explanation goes here

for i=1:size(data,2)
    data(i).x = [data(i).x(:,1:2),zeros(size(data(i).t)),data(i).x(:,3:6),zeros(size(data(i).t)),data(i).x(:,7:8)];
end

F_complete = getAnimationHopper5DOF(data,epsilon,full3D);
end