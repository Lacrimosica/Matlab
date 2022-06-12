% This function implements a Heaviside Lambda signal
%
% The input parameters are:
% t --> time, it can be an array
% T --> the time-width at half height of the triangle
%
% The output parameter is:
% y --> it has the same size as t
%
function [y] = HLambda(T,t)
    y= (1-t/T).*(t>0 & t<T) +...
       (1+t/T).*(t>-T & t<0)+1*(t==0); 
end