%% Problhma gia kyklwma RLC 

clear;
clc;
close all;

tic
% taseis eisodou twn 2 phgwn
u1 = @(t) 2 * sin(4*t);
u2 = 4;

% arxikh syn8hkh systhmatos
[x0(1), x0(2)] = v(0);

% deigmatolhpsia ana 0.00001s gia 0-50s diarkeia sto peirama gia na fanei olo to metabatiko fainomeno
step = 1e-5;
t_linspace = 0:step:50;
N = length(t_linspace);
y_out = zeros(N, 2);

%% antlhsh metrhsewn apo to arxeio gia V_C kai V_R
for i = 1:N
    [V_R, V_C] = v(t_linspace(i));
    y_out(i, 1) = V_C;
    y_out(i, 2) = V_R;
end

V_C = y_out(:, 1);

% pros8hkh 3 texnhtwn outliers kata megalyterh tash gia b erwthma
% V_C(100000) =  V_C(100000) + 50 *  V_C(100000);
% V_C(300000) =  V_C(300000) + 50 *  V_C(300000);
% V_C(500000) =  V_C(500000) + 50 *  V_C(500000);

U1 = double(u1(t_linspace))';
U2 = ones(N, 1) .* u2;
%% Methodos elaxistwn tetragwnwn(LSM) 
% oi poloi tou eusta8ous filtrou
lambda = [50, 80];

% to dianysma parametrwn gia to montelo dynamikou systhmatos meta thn grammikh parametropoihsh
z1 = lsim(tf([-1, 0], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), V_C, t_linspace');
z2 = lsim(tf([-1], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), V_C, t_linspace');
z3 = lsim(tf([1, 0], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U2, t_linspace');
z4 = lsim(tf([1], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U2, t_linspace');
z5 = lsim(tf([1, 0], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U1, t_linspace');
z6 = lsim(tf([1], [1, lambda(1) + lambda(2), lambda(1) * lambda(2)]), U1, t_linspace');
zeta = [z1, z2, z3, z4, z5, z6]';

% ektimhsh twn agnwstwn parametrwn
theta = optimal_params(y_out, zeta) + [lambda(1) + lambda(2); lambda(1) * lambda(2); 0; 0; 0; 0];
RC_inv = theta(1)
LC_inv = theta(2)


%% Lysh ODE apo to ektimwmeno montelo systhmatos gia na bgalw tis metrhseis eisodou/eksodou gia ka8e xronikh stigmh t
[t, V_C_estimation] = ode45(@(t, Vc_approx) system_equations(t, Vc_approx, RC_inv, LC_inv, u1, u2), t_linspace, x0);
sum = 0;
for i = 1: N
    sum = sum + abs((V_C(i) - V_C_estimation(i, 1)));
end
V_C_total_error = sum / N

%V_R = y_out(:, 2); den douleuei gia outliers
V_R = U1 + U2 - V_C;
V_R_estimation = U1 + U2 - V_C_estimation(i, 1);
for i = 1: N
    sum = sum + abs((V_R(i) - V_R_estimation(i)));
end
V_R_total_error = sum / N
%% grafikes parastaseis twn metrhsewn tou arxeiou kai tou ektimwmenou montelou gia ton pyknwth kai thn antistash
figure(1);
plot(t_linspace, V_C, 'Linewidth', 0.1);
ylabel('Output $y(t) = V_C(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

figure(2);
plot(t_linspace, V_C_estimation(:, 1), 'Linewidth', 0.1);
ylabel('Estimated output $\hat{V_C}(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

figure(3);
plot(t_linspace, V_R, 'Linewidth', 0.1);
ylabel('Voltage of resistance $V_R(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

figure(4);
plot(t_linspace, V_R_estimation, 'Linewidth', 0.1);
ylabel('Estimated voltage of resistance $\hat{V_R}(t)$', 'interpreter', 'latex');
xlabel('$t(sec)$', 'interpreter', 'latex');

toc
%% Eksiswseis tou systhmatos
function dx = system_equations(t, x, RC_inv, LC_inv, u1, u2)
    dx(1)  = x(2);
    dx(2) = -RC_inv * x(2) - LC_inv * x(1) + LC_inv * u2 + RC_inv * u1(t);
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
