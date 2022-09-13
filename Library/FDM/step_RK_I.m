function y1 = step_RK_I(y0,t0,h,odefcn,a,b,c)
% Time Steppper for the Runge Kutta Method using an Implicit Integration Scheme 
% 
% Written by: BD Bhai

    %% Allocate Memory 
    N     = length(b);              %[]     Number of Intermediate States Evaluated
    dim   = size(y0);               %[]     Size of the States of the System
    K_smb = sym('k',[dim,N]);       %[]     Symbolic Variable for Intermediate Derivatives
    EQN   = sym(zeros([dim,N]));    %[]     Allocate Memory for Implicit Eqn
    
    %% Determine System of Eqn for K 
    for i = 1:N
        temp = zeros(dim);
        for ii = 1:N
            temp = temp + a(i,ii)*K_smb(:,:,ii);
        end
        
        % Simultaneous Equation for States     [ EQN == 0 , to use with fsolve() ] 
        EQN(:,:,i) = h * odefcn(t0 + c(i)*h, y0 + temp) - K_smb(:,:,i);
    end
    
    %% Solve System of Eqn for k 
    EQN = matlabFunction(EQN,'Vars',{K_smb});   %[]     Parameterize EQN as fcn for use with fsolve()
    IG = ones([1,1,N]) .* odefcn(t0,y0) * h;    %[]     Initial Guess to the Derivative Approx. for K_smb
    
    options = optimoptions('fsolve','Display','off');
    k = fsolve(EQN,IG,options);                 %[]     Nonlinear Optimization Solver for Root Finding
    
    %% Calculate the Next State 
    y1 = y0;                        %[]     Initialize Final State as Initial State
    for i = 1:N
        y1 = y1 + b(i)*k(:,:,i);
    end
    
end