function [a,b,c,solver] = CrankNic()
% Butcher Table with Coefficients for 2nd Order Crank Nicholson Method
% Implicit Trapezoidal Rule that is A-Stable
%
% Written by: BD Bhai
    
    a = 1/2 .* [  0,  0   ; ...
                  1,  1  ]; %[]     Weights for Evaluating Intermediate States
              
    b = 1/2 .* [ 1; 1 ];    %[]     Weights for Evaluating Output State
    
    c = [ 0; 1 ];           %[]     Weights for Evaluating Intermediate Time Steps
    
    solver = @step_RK_I;
end