function [a,b,c,solver] = Imp_Mid()
% Butcher Table with Coefficients for 1st Order Implicit Midpoint Rule
% Symplectic Integrator, so whatever that means
%
% Written by: BD Bhai

    a = 1/2;    %[]     Weights for Evaluating Intermediate States
    b = 1;      %[]     Weights for Evaluating Output State
    c = 1/2;    %[]     Weights for Evaluating Intermediate Time Steps
    
    solver = @step_RK_I;
end