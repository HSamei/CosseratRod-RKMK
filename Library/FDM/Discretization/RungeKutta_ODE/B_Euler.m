function [a,b,c,solver] = B_Euler()
% Butcher Table with Coefficients for 1st Order Implicit Backward Euler Method
%
% Written by: BD Bhai

    a = 1;      %[]     Weights for Evaluating Intermediate States
    b = 1;      %[]     Weights for Evaluating Output State
    c = 1;      %[]     Weights for Evaluating Intermediate Time Steps
    
    solver = @step_RK_I;
end