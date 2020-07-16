function [fx, allX, x_k] = conjugateGradient(f, x0)
tic
syms alp
x = sym('x', [1 length(x0)]);

allX = [];                                                    % saves all x's during minimization
fx = [];                                                      % saves all functional values
criteria = [];                                                % vector to save the satisfaction of two successive criterions

gradF = gradient(f,x);                                        % gradient of f wrt x
gradF = matlabFunction(vpa(gradF), 'vars', {x});              % convert gradF to a function handle
f = matlabFunction(vpa(f), 'vars', {x});                      % convert f to a function handle

iter = 0;                                                     % iteration counter
x_k = x0;
criterion = 0;

while 1
    iter = iter+1;
    criteria(1) = criterion;
    
    allX(iter,:) = x_k;
    fx(iter) = f(x_k);
    
    if iter==1
        g_new = -gradF(x_k);
        d = g_new;
    else
        g_old = g_new;
        g_new = -gradF(x_k);
        beta = (g_new'*g_new)/(g_old'*g_old);
        d = g_new + beta*d;
    end
    
    x_hat = x_k + alp*d';                                      % introduce alpha for next step in the direction of d
    f_gold(alp) = f(x_hat);
    f_gold_h = matlabFunction(f_gold);
    [~, alp_opt] = GoldenSection(0,1e-6,f_gold_h,1);           % goldensection to find alpha optimal
    
    fprintf('ITERATION %d\n', iter)
    fmt = ['x_k: [', repmat('%.16g, ', 1, numel(x_k)-1), '%.16g]\n'];
    fprintf(fmt, x_k)
    fprintf('--------------------------------------------------------------------------------------------\n');
    
    f_old = double(f(x_k));
    x_k = x_k + alp_opt*d';                                    % find the next x in the direction of d
    f_new = double(f(x_k));                                    % functional value of f at x_k+1
   
    criterion = abs(f_new - f_old) <= 1e-6 + 1e-6*abs(f_old);                  % convergence criteria checker
    criteria(2) = criterion;
    grad_criterion = norm(gradF(x_k))<=1e-6;
    
    if criteria(1) && criteria(2) && grad_criterion
        break
    end
end
toc
end
