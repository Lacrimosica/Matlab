%small pi function

function [y] = smallpi (T,t)
  y = HPi (t .- T/2);
end
