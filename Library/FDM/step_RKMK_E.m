function y1 = step_RKMK_E(y0,t0,h,f,a,b,c)
% Time Steppper for the Runge-Kutta-Munthe-Kaas Method using an Explicit Integration Scheme 
% on a Matrix Lie Group with Matrix Exponential Mapping to Group from Algebra.
% 
% Written by: BD Bhai

    %% Allocate Memory
    dim = size(y0);
    u   = zeros([ dim, length(b) ]);
    k   = zeros([ dim, length(b) ]);
    kd  = zeros([ dim, length(b) ]);
    v   = zeros(dim);

    %% General Program for Matrix Manifold
    for i = 1: length(b)
        for ii = 1:i-1
            u(:,:,i) = u(:,:,i) + h*a(i,ii)*kd(:,:,ii);
        end
        k(:,:,i)  = f(h*c(i) + t0, expm(u(:,:,i))*y0);
        kd(:,:,i) = dexpinv(u(:,:,i), k(:,:,i));
        v = v + h* b(i)*kd(:,:,i);
    end
    y1 = y0 * expm(v);
    
end