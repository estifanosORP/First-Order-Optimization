% Requires a function handle "f", an initial guess "x1", and a step size "s"
% x_all is an array of all the values of x_min at each iteration
% x_min is the optimal value of x that minimizes the objective function f
% cr ={1,2,3} denotes the criterion that must be used for convergence checking
function [x_all, x_min] = GoldenSection(x1,s,f,cr)
r = (sqrt(5)+1)/2;
tau = 1/r;

x2 = x1 + s;

if f(x2)>f(x1)
    s = -s;
    x1=x2;
    x2=x1;
end

s = r*s;
x4 = x2+s;

x_all = [];
while 1
    while f(x4) <= f(x2)
        x1 = x2;
        x2 = x4;
        s = r*s;
        x4 = x2+s;
    end
    
    x3 = (1-tau)*x1 + tau*x4;
    
    f_bar_old = (f(x1)+f(x2)+f(x3))/3;  % to be used later for convergence checking
    if f(x2)<f(x3)
        x4 = x1;
        x1 = x3;
    else
        x1 = x2;
        x2 = x3;
    end
    x_all = [x_all, x2];
    eps_rel = 10e-16; eps_abs = 10e-4;
    f_bar = (f(x1)+f(x2)+f(x3))/3;      % new f_bar
    criterion_1 = abs(x1-x4) <= eps_rel*abs(x2) + eps_abs;
    criterion_2 = abs(f_bar - f_bar_old) <= eps_rel*abs(f(x2)) + eps_abs*10e-2;
    
    if cr == 1
        if criterion_1
            break
        end
    elseif cr == 2
        if criterion_2
            break
        end
    elseif cr == 3
        if criterion_1 && criterion_2
            break
        end
    end
end
x_min = x2;
end
