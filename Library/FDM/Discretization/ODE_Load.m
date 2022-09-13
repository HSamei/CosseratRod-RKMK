function [a,b,c,solver] = ODE_Load(Discretization)
%[]     Adopt ODE Discretization Scheme

    if     strcmpi(Discretization,'RK4')    
        [a,b,c,solver] = RK4();
    elseif strcmpi(Discretization,'RK2')
        [a,b,c,solver] = RK2();
    elseif strcmpi(Discretization,'Euler')
        [a,b,c,solver] = Euler();
    elseif strcmpi(Discretization,'B_Euler')
        [a,b,c,solver] = B_Euler();
    elseif strcmpi(Discretization,'CrankNic')
        [a,b,c,solver] = CrankNic();
	elseif strcmpi(Discretization,'Imp_Mid')
        [a,b,c,solver] = Imp_Mid();
    else 
        error('Get Help');
    end
    
end