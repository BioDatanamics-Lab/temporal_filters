function [soln, timecourse, isFired] = numerical_dscPois(ODEparams,regiemes)

    % ODE params
    C = ODEparams.C; E_L = ODEparams.E_L; g_L = ODEparams.g_L; I_app = ODEparams.I_app;
    G_ex = ODEparams.G_ex; E_ex = ODEparams.E_ex; tau_dec = ODEparams.tau_dec;
    tau_dep = ODEparams.tau_dep; tau_fac = ODEparams.tau_fac; U = ODEparams.U;
    
    % Facilitation and Depression ONLY models   
    hasDep = 1; hasFac = 1;
    if tau_dep == Inf
        % (dep=0/fac=1 triggers depression only model) => dR/dt = 0 => dep >> fac
        hasDep = 0;
    end
    if tau_fac == Inf
        % (fac=0/dep=1 triggers facilitation only model) => du/dt = 0 => fac >> dep
        hasFac = 0;
    end
    if tau_fac == 0 
       hasFac = 0; tau_fac = 1; %zero out compoentns, set = 1 so there's no division by zero
    end
    if tau_dep == 0
       hasDep = 0; tau_dep = 1; %zero out compoentns, set = 1 so there's no division by zero
    end
    
    % y(1) - R, y(2) - u, y(3) - S, y(4) - V
    f = @(t,y) [ ...
        ((1-y(1))/tau_dep) * hasDep                            ;
        ((U-y(2))/tau_fac) * hasFac                            ;
        -y(3)/tau_dec                                          ;
        (-1*g_L*(y(4)-E_L)+I_app-G_ex*y(3)*(y(4)-E_ex))/C     ];

    if g_L == 0
        Vinit = E_L;
    else
        Vinit = E_L; % + I_app/g_L; % set intial condition to leak asymptote
    end
    t(1) = 0; y(:,1) = [1 U 0 Vinit];

	% N is the number of timesteps in the actual RK4 simulation
	% h is the Delta on the RK4 timesteps.
    h = .01;
    
    [N, fire] = createPoisFiringTimes(regiemes, h);

    for i = 1:N
        % RK4 Steps
        k1 = f(t(i)      , y(:,i)             );
        k2 = f(t(i) + h/2, y(:,i) + h/2 * k1/2);
        k3 = f(t(i) + h/2, y(:,i) + h/2 * k2/2);
        k4 = f(t(i) + h  , y(:,i) + h   * k3  );
        % Update
        t(i+1) = t(i) + h;
        y(:,i+1) = y(:,i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
        % Firing Logic
        if fire(i) == 1
		    % R / u - update U first, then R with new U value
            oldR = y(1,i+1); oldu = y(2,i+1);
            y(2,i+1) = oldu + U*(1-oldu);% * hasFac;        %chance that Ca2+ channels are open increases
            y(1,i+1) = oldR - oldR * y(2,i+1);% * hasDep;   %then presynaptic resources decrease
            y(3,i+1) = y(3,i+1) + oldR * y(2,i+1);          %(old R, new u)
        end
		% tau dep/fac = 0 require different firing logic.
        if i > 1 && fire(i-1) == 1 && ODEparams.tau_dep == 0
            y(1,end) = y(1,1);
        end
        if i > 1 && fire(i-1) == 1 && ODEparams.tau_fac == 0
            y(2,end) = y(2,1);
        end
        
    end
    soln = y;
    timecourse = t;
    isFired = fire;
end


function [N, fire] = createPoisFiringTimes(regiemes, h)

    fire = [];
    for i = 1:size(regiemes,1)
        simulationSteps = (regiemes(i,2)/h);
        sim = zeros(1,simulationSteps); %train
        event_time = [];
        for j = 1:simulationSteps
            if(rand()<(regiemes(i,1)*h)/1000) %divide by 1000 to adjust units [kHz] CHECK THIS!
                sim(j) = 1;
                event_time = [event_time j*h];
            end
        end
        %{
        % check that this is exponential -- if so, then simulation is pois.
        intervals = diff(event_time); hist(intervals)
        pause
        close
        %}
        fire = [fire sim];
    end
    N = length(fire);
end

