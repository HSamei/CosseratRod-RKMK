function [c_sd,offset_sd] = SD_Load(Discretization,h_SD)
%[]     Adopt Semi-Discretization Scheme
% 
%   c_sd:       Weights for corresponding steps
%   offset_sd:  Index corresponding to current step
% 
% Written by: BD Bhai

    if     strcmpi(Discretization,'BDF1')       
        [c_sd,offset_sd] = BDF1(h_SD);
    elseif strcmpi(Discretization,'BDF2')
        [c_sd,offset_sd] = BDF2(h_SD);
    elseif strcmpi(Discretization,'CentDiff')
        [c_sd,offset_sd] = CentDiff(h_SD);
    else
        error('Get Help')
    end
    
end