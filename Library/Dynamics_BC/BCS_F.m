function Error = BCS_F(InitGuess, BODY, JOINT, THETA, THETA_DOT, THETA_DDOT, F_ext, F_dst, dt, LA_SemiDsc, LA_ODE, LG_ODE, Intrp_Fcn)
% Cosserat rod Boundary Condition Solver to compute the straint and applied wrench at EE
%
% DETERMINE:    Error in the F_ext implied from assumed InitGuess
% GIVEN:        Manipulator Definition & Joint Position, Velocity, Acceleration &
%               Applied Loading & FDM Discretization Coeff.
%
%     - InitGuess:      Initial Guess for Strain at s = 0
%     - BODY:           Definition of Flexible Body
%     - JOINT:          Definition of Joint
%     - THETA:          Actuated Joint Position
%     - THETA_DOT:      Actuated Joint Velocity
%     - THETA_DDOT:     Actuated Joint Acceleration
%     - F_ext:          Applied Wrench at the EE expressed in EE BCF
%     - F_dst:          Distributed Applied Load at the EE expressed in EE BCF
%     - dt:             Time Step Size
%     - LA_SemiDsc:     Lie Algebra Semi-Discretization Scheme
%     - LA_ODE:         Lie Algebra ODE Integration Scheme
%     - LG_ODE:         Lie Group ODE Integration Scheme
%     - Intrp_Fcn:      Interpolation Function for Intermediate History Terms
% 
% Written by BD Bhai

    g_base = expm(hat(JOINT.Twist) * THETA);
    F_base = BODY.Stiff * (InitGuess - BODY.F_0);
    eta_base = JOINT.Twist .* THETA_DOT;            %[se(3)]    Velocity Twist at Joint Frame
    
    [~,f,~] = Flex_Dyn(g_base, F_dst, F_base, BODY, eta_base, dt, LA_SemiDsc, LA_ODE, LG_ODE, Intrp_Fcn);

    F_tip = BODY.Stiff * (f(:,:,end) - BODY.F_0);
    Error = F_ext - F_tip;
    
end                     %[]     FUNCTION END