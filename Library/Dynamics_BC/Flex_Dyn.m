function [g,f,eta,d_eta] = Flex_Dyn(g_base, F_dst, F_base, BODY, eta_base, dt, LA_SemiDsc, LA_ODE, LG_ODE, Intrp_Fcn)
% Integrate the Cosserat PDE for a Flexible Body to return the states and transformation at its end.
%#ok<*ASGLU> 
% 
% DETERMINE:    Configuration, Acceleration @ End & Velocity, Strain over Body
% GIVEN:        Transformation, Velocity, Wrench BC @ Base, FDM Coeff, Body Definition, Applied Loads
% 
%     - g_base:     Transformation to the Base Frame of the Body wrt to Reference Frame
%     - F_dst:     Applied Distributed Wrench on the Body (expressed in Body Frame)
%     - F_base:     Twist for Body velocity of i_th CoM expressed in i_th CoM BCF
%     - BODY:       Body Object with Relevant Informaiton
%     - eta_base:   Velocity Twist of Base Frame
%     - c0:         FDM Coeff - Current Time Step
%     - c1:         FDM Coeff - Previous Time Step
%     - c2:         FDM Coeff - PrePrevious Time Step
%     - d_eta:      Body Acceleration Twists
%     - eta(_h):    Body Velocity Twists (Finite Difference Approximation using previous values)
%     - f(_h):      Body Strain Twists   (Finite Difference Approximation using previous values)
%     - eta_s:      Velocity Spatial Rate Twists
%     - f_s:        Strain Spatial Rate Twists
%     - f_sh:       Strain Spatial Rate Twist (Finite Difference Approximation using previous values)
% 
% Written by BD Bhai
    
    %% Allocate Memory
%     d_eta = zeros(6,BODY.N);    %[se(3)]    Body Acceleration Twist at Local Cross-Section
    eta = zeros(6,1,BODY.N);      %[se(3)]    Body Velocity Twist at Local Cross-Section
    f   = zeros(6,1,BODY.N);      %[se(3)]    Body Strain Twist at Local Cross-Section
    g   = zeros(4,4,BODY.N);      %[SE(3)]    Transformation to Local Cross-Section
    
    %% Apply Semi-Discretization
    [c_sd, offset_sd] = SD_Load(LA_SemiDsc, dt);    %[]     GENERALIZE somehow 
    eta_h = zeros(6, 1, BODY.N);
    f_h   = zeros(6, 1, BODY.N);
    
    for i = (1 + offset_sd):length(c_sd)            %[]     Because c_sd(offset_sd) is current step, so prev from (offset_sd + 1)
        eta_h = eta_h + c_sd(i)*BODY.eta_prev(:,:,:,i - offset_sd);
        f_h   =   f_h + c_sd(i)*BODY.f_prev(:,:,:,i - offset_sd);
    end
    
    %% Integrate Coss ODE Numerically
    % Assign Integrators and Interpolators
    [a,b,c,solver] = ODE_Load(LA_ODE);              %[]         Load Lie Algebra ODE Integrator
    solver = @step_RK_E_h;                          %[]         GENERALIZE, or just restrict to only explicit LA_ODE dsc
    
    [N_intrpl,offset_intrpl,Intrpl] = Intrpl_Load(Intrp_Fcn);
    Y_h = zeros(12,1,N_intrpl);                     %[]         Allocate Memory for Interpolation States
    
    % Initialize States
    f(:,1,1) = BODY.Stiff \ F_base + BODY.F_0;      %[se(3)]    Assign Initial Strain Twist @ Base of Body
    eta(:,1,1) = eta_base;                          %[se(3)]    Assign Initial Velocity Twist @ Base of Body
    g(:,:,1) = g_base;                              %[SE(3)]    Assign Configuration @ Base of Body wrt to Reference Frame
    ds = BODY.L / (BODY.N - 1);                     %[m]        Spatial Step Size assumed in Numerical Integration
    
    for i = 1 : BODY.N - 1
        %% 2nd Order Derivative Approx. GENERALIZE (hardcoded with Euler on top of LA_SemiDsc)
        f_sh = zeros(6,1);
        for ii = 2:length(c_sd)
            f_sh = f_sh + ( c_sd(ii)*(BODY.f_prev(:,:,i+1,ii-offset_sd) - BODY.f_prev(:,:,i,ii-offset_sd)) ) / ds;
        end
        odefcn_h = @(s,y,y_h)Coss_ODE_Dsc(y, y_h, f_sh, BODY, c_sd(offset_sd), F_dst);
        
        %% APPLY INTERPOLATED STEPPER ON LA_ODE
        % Define History States for Interpolation
        for ii = 1:N_intrpl 
            index = i + ii - offset_intrpl;                         %[]     Spatial Index (if within bounds)
            if index < 1                                            %[]     If before first node
                Y_h(:,:,ii) = [f_h(:,:,1); eta_h(:,:,1)];
            elseif index > BODY.N                                   %[]     If after last node
                Y_h(:,:,ii) = [f_h(:,:,end); eta_h(:,:,end)];
            else
                Y_h(:,:,ii) = [f_h(:,:,index); eta_h(:,:,index)];
            end
        end %[] Define Y_h for interpolation
        
        y0 = [f(:,:,i); eta(:,:,i)];                    %[]     Collection of State Variables at Current Step
        y1 = solver( y0, Y_h, (i-1)*ds, ds, Intrpl, ...
                     odefcn_h, a, b, c);                %[]     Collection of State Variables at Next Step
        
        % Define States at Next Step
        f(:,:,i+1)   = y1(1:6);                         %[se(3)]    Partition Integrated Results
        eta(:,:,i+1) = y1(7:12);                        %[se(3)]    Partition Integrated Results
        
        %% APPLY INTERPOLATED STEPPER ON LG_ODE
        [a_G,b_G,c_G,solver_G] = ODE_Load(LG_ODE);
        Str_Intrpl = @(s,g)hat(Intrpl(f(:,:,i:i+1), s/ds));     %[]   GENERALIZE (assuming linear interpolation)
        
        g(:,:,i+1) = step_RKMK_E( g(:,:,i), ds*(i-1), ds, Str_Intrpl, a_G, b_G, c_G );
    end
    
end