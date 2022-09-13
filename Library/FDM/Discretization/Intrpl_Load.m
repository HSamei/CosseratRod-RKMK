function [N,offset_intp,Intrpl] = Intrpl_Load(Discretization)
%[]     Adopt Interpolation Scheme
% 
%   N:              Number of nodes needed for interpolation
%   offser_intp:    Index corresponding to current step
% 
% Written by: BD Bhai

    if strcmpi(Discretization,'Cubic_Intrpl')
        N = 4;
        offset_intp = 2;
        Intrpl = @Cubic_Intrpl;
    elseif strcmpi(Discretization,'Linear_Intrpl')
        N = 2;
        offset_intp = 1;
        Intrpl = @Linear_Intrpl;
    else
        error('Get Help')
    end
    
end