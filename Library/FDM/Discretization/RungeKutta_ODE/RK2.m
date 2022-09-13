function [a,b,c,solver] = RK2()
% Butcher Table with Coefficients for 2nd order Explicit Runge Kutta Method
%
% Written by: BD Bhai

    a = [  0,  0   ; ...
           1,  0  ];        %[]     Weights for Evaluating Intermediate States
       
    b = 1/2 .* [ 1; 1 ];    %[]     Weights for Evaluating Output State
    
    c = [ 0; 1 ];           %[]     Weights for Evaluating Intermediate Time Steps

    solver = @step_RK_E;
end