% This function implements a Heaviside Pi signal
%
% The input parameters are:
% t --> time, it can be an array
% T --> the duration of the signal (the width of the rectangle)
%
% The output parameter is:
% y --> it has the same size as t
%
function [y] = HPi(T,t)
    y=1*((t<T/2) & (t>-T/2))+1/2*((t==-T/2)+(t==T/2));
end