function [c,offset] = BDF1(delta)
% Coefficients for 1st Order Implicit Backward Difference Formula
%
%    c     :  Coefficient for the all indices in sequence (next - curr - prev)
%  offset  :  Index in c corresponding to current step
% 
% Written by: BD Bhai
    
    c = 1/delta .* [ 1; -1 ];   %[]     Weights for Evaluating Intermediate Steps
    offset = 1;                 %[]     Index in c corresponding to current step
    
end