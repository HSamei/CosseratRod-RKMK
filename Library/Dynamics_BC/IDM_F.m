function [F_base, g, BODY_New] = IDM_F(BODY, JOINT, THETA, THETA_DOT, THETA_DDOT, ...
                                F_ext, F_dst, dt, LA_SemiDsc, LA_ODE, LG_ODE, Intrp_Fcn, InitGuess)
% Inverse Dynamic Model for Cosserat Rod
% 
% DETERMINE:    Constraint, Actutation Forces & Body Velocities
% GIVEN:        Manipulator Definition & Joint Position, Velocity, Acceleration
% 
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
%     - InitGuess:      Initial Guess for Strain at s = 0
% 
% Written by BD Bhai
    
    %% IMPLEMENT FDM & BCS
%     c0 = 1.5/dt; c1 = -2/dt; c2 = .5/dt;                %[]     FDM Coefficients for BDF-2
    options = optimset('Display','OFF','TolFun',1e-9);
    
    [InitGuess,~,~,~] = fsolve( @(InitGuess)BCS_F(InitGuess, BODY, JOINT, THETA, THETA_DOT, ...
                                 THETA_DDOT, F_ext, F_dst, dt, LA_SemiDsc, LA_ODE, LG_ODE, Intrp_Fcn), InitGuess, options );
    
    %% RECURSIVE DEFINITION OF DYNAMICS USING EULER-POINCARE EOM
    % Initialize Base
    BODY_New = BODY;                                %[]         Initialize Structure to Save Updated History for BDF-2
    g_base = expm(hat(JOINT.Twist) * THETA);
    eta_base = JOINT.Twist .* THETA_DOT;            %[se(3)]    Velocity Twist at Joint Frame
    F_base = BODY.Stiff * (InitGuess - BODY.F_0);
    
    [g,f,eta] = Flex_Dyn(g_base, F_dst, F_base, BODY, eta_base, dt, LA_SemiDsc, LA_ODE, LG_ODE, Intrp_Fcn);
    
    %% Update History Terms (after all computations complete)
    BODY_New.f_prev(:,:,:,2)   = BODY.f_prev(:,:,:,1);
    BODY_New.f_prev(:,:,:,1)   = f;
    BODY_New.eta_prev(:,:,:,2) = BODY_New.eta_prev(:,:,:,1);
    BODY_New.eta_prev(:,:,:,1) = eta;
end