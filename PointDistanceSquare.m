function [distance] = PointDistanceSquare(x0,y0,x1,y1)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

dx=x1-x0;
dy=y1-y0;
distance=dx*dx+dy*dy;

end

