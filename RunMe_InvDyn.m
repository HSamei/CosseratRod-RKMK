%% %%%%%%%%%%%%% Flexible Body Inverse Dynamic Simlator %%%%%%%%%%%%% %%
% #ok<*CCAT>
% Study different discretizations for geometrically integrating a 
% joint actuated dynamic Cosserat rod.
% 
% Written by BD Bhai

%% LOAD LIBRARY 
addpath(genpath(strcat(pwd,'/Library')))
tic

%% SIMULATION PARAMETERS 
% Define System
Theta       = [0];      %[rad]      Actuated Joint Position
Theta_dot   = [0];      %[rad/s]    Actuated Joint Velocity

% N           = 101;       %[]     Number of Discretized Nodes(!!! REF !!!)
% dt          = .005;      %[s]    Time Step Size (!!! REF !!!)
N           = 13;       %[]     Number of Discretized Nodes (!!! TEST !!!)
dt          = .01;     %[s]    Time Step Size (!!! TEST !!!)
Tsim        = 1.5;      %[s]    Total Duration of Simulation
Tact        = .5;       %[s]    Total Duration of Actuated Joint Trajectory
VisScale    = .5;       %[m]    Scale to Axes Dimensions if System Out of Plot View
NumStep     = 5;        %[]     Factor to shrink ds for Plotting Flexible Bodies

% Define Trajectory over time
t1 = 0:dt:Tact/2-dt;           %[s]    Time Vector 1
t2 = Tact/2:dt:Tact-dt;  %[s]    Time Vector 2
t3 = Tact:dt:Tsim;     %[s]    Time Vector 3

% Define Trajectory over time
Shape       = [30];                     %[rad/s2]   Actuated Joint Acceleration Shape
Theta_ddot1 =  Shape.* ones(size(t1));  %[rad/s2]   Joint Trajectory 1
Theta_ddot2 = -Shape.* ones(size(t2));  %[rad/s2]   Joint Trajectory 2
Theta_ddot3 = 0 .* ones(size(t3));      %[rad/s2]   Joint Trajectory 3
Theta_ddot  = [Theta_ddot1, Theta_ddot2, Theta_ddot3];  %[rad/s2]   Joint Trajectory Combined
t           = [t1,t2,t3];

F_ext       = [0;0;5;0;0;0];        %[N;Nm] Applied Wrench at the EE
F_dst       = [0;-50;0;0;0;0];       %[N;Nm] Distributed Applied Wrench over the Body
F_0         = [1;0;0;0;0;0];        %[]     Free Strain in first Flexible Body
SNAP        = [1,13,22,61,0];      %[]     TimeStep to Plot Snapshot (ascending order and padd with 0)
filename    = 'TestPaper_1_Inv.gif';%[]     Filename for Output Video

%% DEFINE BODY AND DISCRETIZATION FOR SOLUTION 
BODY       = defPaperSample_1(N);
JOINT      = Joint('Joint_1','RIGID',[0;0;0;0;0;1],0,0,0,[-pi,pi],0,null(1),BODY);

LA_SemiDsc = 'BDF2';            %[]     BDF1   or  BDF2
LA_ODE     = 'RK4';             %[]     Euler  or  RK2  or  RK4
LG_ODE     = LA_ODE;            %[]     Adopt same as LA_ODE for RKMK
Intrp_Fcn  = 'Linear_Intrpl';   %[]     Cubic_Intrpl  or  Linear_Intrpl
Animate    = 1;                 %[]     Generate Animation Video and File (1)

% Allocate memory for the history structures
[c_sd,offset_sd] = SD_Load(LA_SemiDsc,dt);
BODY.eta_prev = zeros(6, 1,BODY.N, length(c_sd) - offset_sd);
BODY.f_prev = BODY.F_0 .* ones(1,1, BODY.N, length(c_sd) - offset_sd);

[a,b,c,solver] = ODE_Load(LA_ODE);              %[]         Load Lie Algebra ODE Integrator (for theta integration maybe) Fuck it

%% ALLOCATE MEMORY FOR RESPONSE LOOP
T_H  = zeros(1,length(t));         %[N,Nm] Allocate Memory to save the Joint Angles
C_H  = zeros(1,length(t));         %[N,Nm] Allocate Memory to save the Joint Torques
g_H  = zeros(4,4,length(t));       %[N,Nm] Allocate Memory to save the Transformations
Td_H = zeros(1,length(t));         %[N,Nm] Allocate Memory to save the Joint Velocity
EE_POS = zeros(3,length(t));       %[N,Nm] Allocate Memory to save the EE Position

% Initialization of Objects for Animation
if Animate 
    % Define Figure size and View Angle
    Fig = figure(21);
    view(3), axis equal, axis(VisScale*[-1 1 -1 1 -.01 1]), grid on, Axes = gca;
    view(45,20)                    %[]     LOOK AT XY PLANE
    
    % Define Cell Array to Store Handle for each Body
    PlotObj = animatedline(Axes,'Color','r','LineStyle','none','Marker','.','MarkerSize',10);
end

%% RUN LOOP FOR TIME RESPONSE
for i = 1:length(t) 
    %% UPDATE : For next time step
    [F_base, g, BODY_New] = IDM_F(BODY, JOINT, Theta, Theta_dot, Theta_ddot(i), F_ext, F_dst, dt, LA_SemiDsc, LA_ODE, LG_ODE, Intrp_Fcn, F_0);
    
    % Save variables for post-processing
    T_H(:,i)    = Theta;
    C_H(:,i)    = F_base' * JOINT.Twist;
    Td_H(:,i)   = Theta_dot;
    g_H(:,:,i)  = g(:,:,end);
    
    % Joint Values Updated (assumes Euler integration for now)
    Theta = Theta + Theta_dot*dt;
    Theta_dot = Theta_dot + Theta_ddot(i)*dt;
    
    %% ANIMATE :
    ds = BODY.L / (BODY.N-1);                   %[m]    Size of Differential Element (used to compute f_prev)
    POS = zeros(3,(BODY.N - 1) * NumStep);            %[]     Allocate Memory for CS Positions
    g_temp = expm(hat(JOINT.Twist) * T_H(:,i));
    
    for ii = 1:(NumStep*(BODY.N - 1))
        index = ceil(ii/NumStep);                                       %[]         Determine the Local Strain Value at Step
        g_temp = g_temp* expm3(hat(BODY.f_prev(:,index)) * ds/NumStep); %[SE(3)]    Transformation to s = s + ds/NumStep
        POS(:,ii) = g_temp(1:3,4);                                      %[R3]       Save Position at s = s + ds/NumStep, index to next
    end
    
    if Animate
        clearpoints(PlotObj)
        addpoints(PlotObj,POS(1,:),POS(2,:),POS(3,:));
        drawnow;
        im = frame2im(getframe(Fig));
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,filename,'gif','Loopcount',inf);
        else
            imwrite(imind,cm,filename,'gif','DelayTime',dt,'WriteMode','append');
        end
    end
    
    %% UPDATE BODY OBJECT & INITIAL GUESS
    F_0 = BODY_New.f_prev(:,1,1); %[]    Use previous solution as initial guess
    BODY = BODY_New;              %[]    Update Body history states
    EE_POS(:,i) = POS(:,end);
end

%% POST-PROCESSING 
TOGG_PLOT = 1;
if TOGG_PLOT
    %% Save C_des 
%     if saveControl
%         C_des = C;
%         save(saveName,'C_des')  %[] Save to Output File
%     end
    
%     %% Plot Joint Actuation
%     figure()
%     plot(t,C_H,'r','LineWidth',2),hold on, grid on
%     xlabel('Time(s)')
%     ylabel('Joint Act (Nm or N)')

%     %% Plot Joint Angles
%     figure()
%     plot(t,T_H,'r','LineWidth',2), grid on
%     xlabel('Time(s)')
%     ylabel('Joint Angle')

    %% Plot EE Location
    figure()
    plot(t,EE_POS(1,:),'r','LineWidth',2),hold on, grid on
    plot(t,EE_POS(2,:),'k','LineWidth',2)
    plot(t,EE_POS(3,:),'Color',[.5,.5,.5],'LineWidth',2)

    xlabel('Time(s)')
    ylabel('EE Position')
    
end

%% DEBUG PLAYGROUND
toc