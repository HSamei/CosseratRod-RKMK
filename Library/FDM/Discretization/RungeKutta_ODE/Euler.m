function [a,b,c,solver] = Euler()
% Butcher Table with Coefficients for 1st Order Explicit Euler Method
%
% Written by: BD Bhai

    a = 0;      %[]     Weights for Evaluating Intermediate States
    b = 1;      %[]     Weights for Evaluating Output State
    c = 0;      %[]     Weights for Evaluating Intermediate Time Steps
    
    solver = @step_RK_E;
end