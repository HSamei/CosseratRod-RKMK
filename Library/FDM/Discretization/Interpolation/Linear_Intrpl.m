function [Y_i] = Linear_Intrpl(Y,mu)
% Linear Interpolation between 2 Vectors
% 
% DETERMINE:    Interpolated State Between States Y
% GIVEN:        States, Location Between them
% 
%     - Y:      States at Start and End of Interpolation
%     - Y_i:    Interpolated States
%     - mu:     Location Between Start and End (0 = Start, 1 = End)
%
% Written by BD Bhai
    
    if length(Y(1,1,:)) ~= 2
        error('Please reconsider your life decisions')
    end
    
    Y_i = Y(:,:,1) + (Y(:,:,2) - Y(:,:,1))*mu;
    
end