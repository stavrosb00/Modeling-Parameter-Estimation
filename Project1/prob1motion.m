%% Problhma gia kinhsh mazas me ekswterikes dynameis kai elathrio + aposbesthra

clear;
clc;
close all;

tic
% c erwthma gia parametrous systhmatos kai ekswterikh dynamh eisodou u
m = 10;
b = 0.5;
k = 2.5;
u = @(t) 15 * sin(3*t) + 8; 

% mhdenikes arxikes syn8hkes systhmatos
x0(1) = 0;
x0(2) = 0;

% deigmatolhpsia ana 0.1s gia 10s diarkeia sto peirama
step = 0.1;
t_linspace = 0:step:10;


%% Lysh ODE apo to 8ewrhtiko systhma gia na bgalw tis metrhseis eisodou/eksodou gia ka8e xronikh stigmh t
[t, y_out] = ode15s(@(t, y_out)system_equations(t, y_out, m, b, k, u), t_linspace, x0);                                       
U = u(t(:));
Y = y_out(:, 1);

%% Methodos elaxistwn tetragwnwn(LSM) 
% oi poloi tou eusta8ous filtrou
lambda = [0.5, 0.5];

% to dianysma parametrwn gia to montelo dynamikou systhmatos meta thn grammikh parametropoihsh
syms m_ b_ k_; 
theta_lambda = [b_ / m_ - (lambda(1) + lambda(2)); k_ / m_ - (lambda(1) * lambda(2)); 1 / m_];

% to dianysma eisodou gia to montelo dynamikou systhmatos meta thn grammikh parametropoihsh 
% lsim(sys,u,t)
z1 = lsim(tf([-1, 0], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), y_out(:, 1), t);
z2 = lsim(tf(-1, [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), y_out(:, 1), t);
z3 = lsim(tf(1, [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U, t);                                        
zeta = [z1, z2, z3]';

% ektimhsh twn agnwstwn parametrwn
op = optimal_params(y_out, zeta);
eqns = (theta_lambda == op);
solved_params = solve(eqns, [m_, b_, k_]);
m_ = double(solved_params.m_);
b_ = double(solved_params.b_);
k_ = double(solved_params.k_);

fprintf("Least squares method calculated the following parameters of the model: m = %d, b = %d, k = %d\n", m_, b_, k_);
%% Lysh ODE apo to ektimwmeno montelo systhmatos gia na bgalw tis metrhseis eisodou/eksodou gia ka8e xronikh stigmh t
[t, y_out_estimation] = ode15s(@(t, y_out_estimation) system_equations(t, y_out_estimation, m_, b_, k_, u), t_linspace, x0);
error = zeros(length(y_out), 1); 
N = length(y_out);
for i = 1:N
    error(i) = abs((y_out(i) - y_out_estimation(i)) / y_out(i));
end

%% grafikes parastaseis ths eisodou, ths eksodou tou montelou kai tou error metaksy tous
figure(1);
plot(t, y_out(:, 1), 'Linewidth', 1);
ylabel('Output $y(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

figure(2);
plot(t, y_out_estimation(:, 1), 'Linewidth', 1);
ylabel('Estimated output $\hat{y}(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

figure(3);
plot(t, error, 'Linewidth', 1);
ylabel('Error $\big| \frac{y(t) - \hat{y}(t)}{y(t)} \big|$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');


toc
%% Eksiswseis tou systhmatos 
function dx = system_equations(t, x, m, b, k, u)
    dx(1) = x(2);
    dx(2) = (1 / m) * (-k * x(1) - b * x(2) + u(t));
    dx = dx';
end

%% Gia ka8e sthlh tou zeta briskw tis beltistes parametrous symfwna me thn norma tou MSE 
function theta = optimal_params(y_out, zeta)
    sum_denominator = 0;
    sum_numerator = 0;
    N = length(y_out);
    for i = 1:N
        sum_denominator = sum_denominator + zeta(:, i) * zeta(:, i)';
        sum_numerator = sum_numerator + zeta(:, i) * y_out(i, 1);
    end
    
    theta = sum_denominator \ sum_numerator;
end
