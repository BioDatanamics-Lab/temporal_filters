function [xopt,fopt,niter,gnorm,dx] = gradDFit(delta, grad, alpha, x0)
% delta is the error between envelope and fit
% grad is the gradient wrt parameters
% alpha is step size
% x0 is init conidtions
% only difference with gradDFit is the size of minimum allowed perturbation before stopping fit

f = @(x) delta(x)*delta(x)';

% gradient descent params
% termination tolerance / maximum number of allowed iterations / minimum allowed perturbation
tol = 1e-5; maxiter = 5000000; dxmin = 1e-8;
% initialize gradient norm, optimization vector, iteration counter, perturbation
gnorm = inf; x = x0; niter = 0; dx = inf;

% gradient descent algorithm:
while and(gnorm>=tol, and(niter <= maxiter, dx >= dxmin))
    % calculate gradient:
    g = grad(x);
    gnorm = norm(g);
    % take step:
    xnew = x - .01*alpha'.*g';
    %[f(xnew)^.5 gnorm]
    % check step
    if ~isfinite(xnew)
        display(['Number of iterations: ' num2str(niter)])
        error('x is inf or NaN')
    end
    niter = niter + 1;
    dx = norm(xnew-x);
    x = xnew;
    %f(x)^.5%, pause
end

xopt = x;
fopt = f(xopt);
niter = niter - 1;
