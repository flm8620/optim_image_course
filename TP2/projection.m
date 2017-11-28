function [ xp ] = projection( x )
%PROJECTION Summary of this function goes here
%   Detailed explanation goes here
r = norm(x(1:2));
scale = 0;
if r>1
    scale = 1/r;
else
    scale = 1.0;
end
x(1:2)=x(1:2)*scale;
xp=x;
end

