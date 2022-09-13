function [a,b,c,solver] = RK4()
% Butcher Table with Coefficients for 4th order Explicit Runge Kutta Method
%
% Written by: BD Bhai

    a = [  0,  0,  0,  0   ; ...
          .5,  0,  0,  0   ; ...
           0, .5,  0,  0   ; ...
           0,  0,  1,  0  ];        %[]     Weights for Evaluating Intermediate States
       
    b = 1/6 .* [ 1; 2; 2; 1 ];      %[]     Weights for Evaluating Output State
    
    c = 1/2 .* [ 0; 1; 1; 2 ];      %[]     Weights for Evaluating Intermediate Time Steps

    solver = @step_RK_E;
end