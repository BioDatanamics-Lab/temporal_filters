function [soln, timecourse, isFired] = numerical(ODEparams,spikePer,tFin)

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

    Vinit = E_L;
    t(1) = 0; y(:,1) = [1 U 0 Vinit];
    h = .01;

    % ensures we have enough data for transients -- go out at least 50
    % spikes.
    N = ceil(tFin/h);
    if N*h/spikePer < 50
        N = ceil(spikePer*50/h);
    end
        
    fire = zeros(1, N+1);
        
    for i = 1:N
        % RK Steps
        k1 = f(t(i)      , y(:,i)             );
        k2 = f(t(i) + h/2, y(:,i) + h/2 * k1/2);
        k3 = f(t(i) + h/2, y(:,i) + h/2 * k2/2);
        k4 = f(t(i) + h  , y(:,i) + h   * k3  );
        % Update
        t(i+1) = t(i) + h;
        y(:,i+1) = y(:,i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
        % Firing Logic
        if  (t(i)/spikePer)-floor(t(i)/spikePer) < h/(spikePer)
        % R / u - update U first, then R with new U value
            oldR = y(1,i+1); oldu = y(2,i+1);
            y(2,i+1) = oldu + U*(1-oldu);% * hasFac;        %chance that Ca2+ channels are open increases
            y(1,i+1) = oldR - oldR * y(2,i+1);% * hasDep;   %then presynaptic resources decrease
            y(3,i+1) = y(3,i+1) + oldR * y(2,i+1);          %(old R, new u)

        fire(i) = 1;
        end

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
