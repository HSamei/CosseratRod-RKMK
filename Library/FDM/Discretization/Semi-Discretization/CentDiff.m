function [c,offset] = CentDiff(delta)
% Coefficients for Explicit 2nd Order Central Difference Formula
%
%    c     :  Coefficient for the all indices in sequence (next - curr - prev)
%  offset  :  Index in c corresponding to current step
% 
% Written by: BD Bhai
    
    c = 1/(2*delta) .* [ 1; 0; -1 ];    %[]     Weights for Evaluating Intermediate Steps
    offset = 2;                         %[]     Index in c Corresponding to Current Step
    
end