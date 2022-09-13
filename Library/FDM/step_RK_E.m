function y1 = step_RK_E(y0,t0,h,odefcn,a,b,c)
% Time Steppper for the Runge Kutta Method using an Explicit Integration Scheme 
% 
% Written by: BD Bhai

    %% Allocate Memory 
    dim = size(y0);
    k   = zeros([ dim, length(b) ]);
    y1  = y0;
    
    %% Compute Next State 
    for i = 1: length(b)
        temp = zeros(dim);
        for ii = 1:i-1
            temp = temp + a(i,ii)*k(:,:,ii);
        end
        
        %% State Solution 
        k(:,:,i) = h * odefcn(t0 + c(i)*h, y0 + temp);
        y1 = y1 + b(i)*k(:,:,i);
    end
    
end