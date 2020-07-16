%% Problem 2-i : Rosenbrok's function
clear all, clc, close all
syms x1 x2

f(x1,x2) = 100*(x2-x1^2)^2 + (1-x1)^2;


%% solution by SGM
% [fx, AllX, x_opt] = steepestDescent(f, [-1.2, 1]);

%% solution by CGM
[fx, AllX, x_opt] = conjugateGradient(f, [-1.2 1]);

fprintf('x* = [%.16f %.16f]\n', x_opt)
fprintf('f(x*) = %.16f\n', fx(end))

x = AllX(:,1); y = AllX(:,2);

figure(1)
plot(fx,'linewidth', 1)
xlabel('Iteration number(k)')
ylabel('f(x_k)')
grid on

figure(2)
plot(x,y,'-r', 'linewidth', 0.5)
xlabel('x-axis')
ylabel('y-axis')
grid on

figure(3)
plot(x,y,'-r', 'linewidth', 0.5)
xlabel('x-axis')
ylabel('y-axis')
hold on
fc = fcontour(f, [min(x) max(x) min(y) max(y)]);
fc.LevelList = fx;
grid on
hold off


