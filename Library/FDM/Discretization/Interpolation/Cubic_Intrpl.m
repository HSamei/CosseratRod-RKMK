function [Y_i] = Cubic_Intrpl(Y,mu)
% Cubic Interpolation between 2 Vectors, also using  preceeding and succeeding points
% 
% DETERMINE:    Interpolated State Between States Y(2,:) and Y(3,:)
% GIVEN:        States, Location Between them
% 
%     - Y:      States at Start and End of Interpolation
%     - Y_i:    Interpolated States
%     - mu:     Location Between Start and End (0 = Start, 1 = End)
%
% Written by BD Bhai
    
    if length(Y(1,1,:)) ~= 4
        error('Please reconsider your life decisions')
    end
    
    a0 = Y(:,:,4) - Y(:,:,3) - Y(:,:,1) + Y(:,:,2);
    a1 = Y(:,:,1) - Y(:,:,2) - a0;
    a2 = Y(:,:,3) - Y(:,:,1);
    a3 = Y(:,:,2);
    
    Y_i = a0*mu^3 + a1*mu^2 + a2*mu + a3;
    
end