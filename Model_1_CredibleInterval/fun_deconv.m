function   ghat  =  fun_deconv(dx,  x,  f, h)

%f: pdf
%x: time span
%h: the result of convolution, which is a cdf 
% ghat: cdf

% Define objective function
objective = @(ghat) sum((conv(f, ghat, 'same') * dx - h).^2);

% Define constraints: ghat should be non-decreasing
A = diag(ones(1,length(x))) + diag(-ones(1,length(x)-1), 1);
A = A(1:end-1, :);
b = zeros(length(x) - 1, 1);

% Lower and upper bounds for CDF
lb =  zeros(length(x), 1)- 1e-8;
ub = ones(length(x), 1) + 1e-8;
g0 = linspace(0, 1 , length(f)) ;

% Use fmincon to find the estimated CDF
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp', 'OptimalityTolerance', 1e-4, 'StepTolerance', 1e-4);
ghat = fmincon(objective, g0, A,b, [], [], lb, ub, [], options);




end