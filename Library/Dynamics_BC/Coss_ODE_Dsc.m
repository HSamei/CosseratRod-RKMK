function [y_s] = Coss_ODE_Dsc(y,y_h,f_sh,Body,c0,F_dst)
% Cosserat Model for a Continuum Semi-Discretized from PDE into a Spatial ODE (!! FLEXIBLE ONLY !!)
% 
% DETERMINE:    Spatial Derivative of Velocity and Strain Twists
% GIVEN:        Velocity, Strain Twists & Stiffness & Mass & Free Strain & Applied Loading
% 
%     - y(_h):      Strain & Velocity Twists (Finite Difference Approximation using previous values)
%     - f_sh:       Strain Spatial Rate Twist (Finite Difference Approximation using previous values)
%     - BODY:       Definition of Flexible Body
%     - c_0:        FDM Coefficient at Current Time Step
%     - F_dst:     Applied Distributed Wrench on the Body (expressed in Body Frame)
%
% Written by BD Bhai
        
       %% REASSIGN VARIABLES TO CONVENTIONAL SYMBOLS 
       % Material Property Assignment
       K   = Body.Stiff;
       C   = Body.Damp;
       M   = Body.Mass;
       f_0 = Body.F_0;
       
       % State Assignment
       f   = y(1:6,:);
       eta = y(7:12,:);
       
       % State History Assignment
       f_h   = y_h(1:6,:);
       eta_h = y_h(7:12,:);
       
       %% Time Discretization (definition of c0 & y_h based on FDM used)
       f_t   = c0*f + f_h;          %[]     Local Time Discretization for history in Local Coordinates
       eta_t = c0*eta + eta_h;      %[]     Local Time Discretization for history in Local Coordinates
       
       %% Spatial Derivatives
       f_s = (K + C*c0) \ ( M*eta_t - (adj(eta)')*M*eta - C*f_sh + (adj(f)')*(K*(f - f_0) + C*f_t) + F_dst );
       eta_s = f_t + adj(eta)*f;
       
       y_s = [f_s; eta_s];

end