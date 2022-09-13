function y1 = step_RK_E_h(y0,Y_h,t0,h,Intrpl,odefcn_h,a,b,c)
% Time Steppper for the Runge Kutta Method using an Explicit Integration Scheme
% Can handle Semi-discretization with Interpolating History terms 
% 
% DETERMINE:    Spatial Derivative of Velocity and Strain Twists
% GIVEN:        Velocity, Strain Twists & Stiffness & Mass & Free Strain & Applied Loading
% 
%     - y1:         Final State
%     - y0:         Initial State
%     - Y_h:        History Terms from Semi-Discretization
%     - y_h:        Interpolation for History Terms Evaluated at a Point
%     - t0:         Initial Time (or whatever you're integrating with respect to really)
%     - h:          Temporal Step Size (or whatever you're integrating with respect to really)
%     - odefcn_h:   ODE function with arguments (time, states, state_h)
%     - odefcn:     ODE function with arguments (time, states)
%     - a:          Weights for RK Integrator
%     - b:          Weights for RK Integrator
%     - c:          Weights for RK Integrator
% 
% Written by: BD Bhai

    %% Allocate Memory 
    dim = size(y0);
    k   = zeros([ dim, length(b) ]);
    y1  = y0;
    
    %% Compute Next State 
    for i = 1: length(b)
        %% Interpolate for History Approximation 
        y_intrpl = Intrpl(Y_h, c(i));               %[]     Verify if c(i)
        odefcn = @(t,y)odefcn_h(t, y, y_intrpl);
        
        %% Compute Progress from Initial State 
        temp = zeros(dim);
        for ii = 1:i-1
            temp = temp + a(i,ii)*k(:,:,ii);
        end
        
        %% State Solution 
        k(:,:,i) = h * odefcn(t0 + c(i)*h, y0 + temp);
        y1 = y1 + b(i)*k(:,:,i);
    end
    
end