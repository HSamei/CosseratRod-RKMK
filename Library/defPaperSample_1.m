function Robot = defPaperSample_1(N)
%#ok<*CCAT>
    
    %% DEFINE THE GEOMETRY OF THE COMPONENTS
    % BEAM PARAMETERS (Thin Rods)
    L_0 = .4;            %[m]        Free Length of the Modeled Beam
    F_0 = [1;0;0;0;0;0]; %[]         Free Static Strain under no Applied Loads 
    rho = 4.5e3;           %[kg/m3]    Density of Material
    mu  = 2e6;           %[N/m^2s]   Viscosity of Peanut Butter
    r = .01;             %[m]        Radius of Cross-Section
    E = 1e8;             %[Pa]       Youngs Modulus
    G = E/(2*(1+.3));    %[Pa]       Shear Modulus
    A = pi*r^2;          %[m2]       Cross-Sectional Area of Beam
    I = pi/4*r^4;        %[m4]       2nd Moment of Inertia of Beam
    
    J = diag([2*I,I,I]);            %[m4]       3D Moment of Inertia of Beam
    Kbt = diag([2*G*I,E*I,E*I]);    %[Nm^2]     Bending and Torsional Rigidity (Rotational)
    Kse = diag([E*A,G*A,G*A]);      %[N]        Shear and Extension Rigidity (Linear)
    Cse = diag([3*A,A,A])*mu;       %[N/s]      Shear and Extension Damping (Linear)
    Cbt = diag([2*I,I,I])*mu;       %[N/s]      Bending and Torsional Damping (Rotational)
    
    K = [Kse,zeros(3);zeros(3),Kbt];
    C = [Cse,zeros(3);zeros(3),Cbt];
    Mf = [rho*A*eye(3),zeros(3);zeros(3),rho*J];
    
    %% DEFINE THE MANIPULATOR
    % Define Bodies
    Robot = Body('Body_1','FLEXIBLE',Mf,K,C,F_0,N,L_0);
    
    % Give Previous States for Flexible Bodies
%     Robot.eta_prev  = zeros(6,Robot.N);
%     Robot.eta_pprev = zeros(6,Robot.N);
%     Robot.f_prev    = ones(1,Robot.N) .* Robot.F_0;
%     Robot.f_pprev   = ones(1,Robot.N) .* Robot.F_0;
    
end