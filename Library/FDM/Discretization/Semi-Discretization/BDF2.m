function [c,offset] = BDF2(delta)
% Coefficients for 2nd Order Implicit Backward Difference Formula
%
%    c     :  Coefficient for the all indices in sequence (next - curr - prev)
%  offset  :  Index in c corresponding to current step
% 
% Written by: BD Bhai
    
    c = 1/delta .* [ 1.5; -2; .5 ];     %[]     Weights for Evaluating Intermediate Steps
    offset = 1;                         %[]     Index in c Corresponding to Current Step

end