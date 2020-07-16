%f is a symbolic multivariable function to be minimized
%x0 is the starting point for direction search
function [fx, allX, x_k] = steepestDescent(f,x0)
tic
syms alp                                                      % variable for 1D search
x = sym('x', [1 length(x0)]);                                 % symbolic variables from x(1) to x(n) for n-variable problem

allX = [];                                                    % saves all x's during minimization
fx = [];                                                      % saves all functional values
criteria = [];                                                % vector that saves the bool for two successive convergence criterions

gradF = gradient(f,x);                                        % gradient of f wrt x
gradF = matlabFunction(vpa(gradF), 'vars', {x});                   % convert gradF to a function handle
f = matlabFunction(vpa(f), 'vars', {x});                           % convert f to a function handle

iter = 0;                                                     % iteration counter
x_k = x0;
criterion = 0;                                                % initial unsatisfied criterion

while 1
    iter = iter+1;
    criteria(1) = criterion;
    
    allX(iter,:) = x_k;
    fx(iter) = f(x_k);
    
    d = -gradF(x_k);
    d = d'/norm(d);                                           % normalized steepest decent direction  
    x_hat = x_k + alp*d;                                      % introduce alpha for next step in the direction of d
    f_gold(alp) = f(x_hat);
    f_gold_h = matlabFunction(vpa(f_gold));
    [~, alp_opt] = GoldenSection(0,1e-6,f_gold_h,1);         % goldensection to find alpha optimal, con_criteria = design variable(1)
    
    f_1 = double(f(x_k));
    x_k = x_k + alp_opt*d;                                    % find the next x in the direction of 
    f_2 = double(f(x_k));                                     % functional value of f at x_k+1
    
    fprintf('ITERATION %d\n', iter)
    fmt = ['x_k: [', repmat('%.16g, ', 1, numel(x_k)-1), '%.16g]\n'];
    fprintf(fmt, x_k)
    fprintf('--------------------------------------------------------------------------------------------\n');
    
    criterion = abs(f_2 - f_1) <= 1e-6 + 1e-6*abs(f_1);         % convergence criteria checker
    criteria(2) = criterion;
%     fprintf('norm of gradient: %.16f \n',norm(gradF(x_k)));
    grad_criterion = norm(gradF(x_k))<=6e-3;                  % if e_g=1e-6 is to be used here, the loop never seems to stop
    
    if criteria(1) && criteria(2) && grad_criterion
        break
    end
end
toc
end
