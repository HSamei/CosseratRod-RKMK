function dei = dexpinv(u,k)
% Numerical Approximation for Inverse of Derivative of Exponential Map
% Uses bernoulli numbers and adjoint map to approximate solution
% 
% Written by: BD Bhai
    
    B    = [ 1; -1/2; 1/6; 0; -1/30; 0; 1/42; 0; -1/30; 0; 5/66 ];     %[]     Bernoulli Numbers
    temp = eye(size(u));
    dei  = zeros(size(u));
    
    for i = 1:length(B)
        dei  = dei + B(i)/factorial(i)*temp*k;
        temp = LieBracket(u,temp);
    end
    
end