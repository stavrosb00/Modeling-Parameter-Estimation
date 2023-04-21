%% Thema 3 : me8odo Lyapunov me meikth domh gia online ektimhsh parametrwn se deuterhs takshs systhma
clear;
clc;
close all;

%% pragmatiko systhma

% eisodos u kai parametroi methodou
u = @(t) 3.5 * sin(7.2 * t) + 2 * sin(11.7 * t);

% System matrices
a11 = -0.25;
a12 = 3;
a21 = -5;
a22 = 0;

b1 = 0.5;
b2 = 1.5;

A = [a11, a12; a21, a22];
B = [b1; b2];

gamma1 = 20;
gamma2 = 25;
theta_m = [-5, 6; -1, 1];

% arxikh syn8hkh
x0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
t_span = 0:0.01:80;


%% prosomoiwsh kai ektimhsh parametrwn 

% lysh diaforikou systhmatos gia xroniko diasthma t_span
[t, x] = ode15s(@(t, x) system_equationsV4(t, x, A, B, gamma1, gamma2, theta_m, u), t_span, x0);

y1 = x(:, 1);
y2 = x(:, 2);
y1_hat = x(:, 9);
y2_hat = x(:, 10);

a11_hat = x(:, 3);
a12_hat = x(:, 4);
a21_hat = x(:, 5);
a22_hat = x(:, 6);

b1_hat = x(:, 7);
b2_hat = x(:, 8);

% telikes ektimwmenes times parametrwn A, B
x(length(t_span), 3:8)

%% aksiologhsh basei metrikwn sfalmatos

% Percentage error
percentage_error = zeros(length(t_span), 1);
for i = 1:length(t_span)
    percentage_error(i) = sqrt(((y1(i) - y1_hat(i)) /  y1(i))^2 + ((y2(i) - y2_hat(i)) /  y2(i))^2);
end

% Mean square error
mean_square_error = zeros(length(t_span), 1);
for i = 1:length(t_span)
    mean_square_error(i) = (y1(i) - y1_hat(i)) ^ 2 + (y2(i) - y2_hat(i)) ^ 2;
end

%% grafikes parastaseis me plots

% kata poso sygklinoyn oi ektimwmenes parametroi tou montelou me to systhma
fig0 = printer_params_est_arrayA(t_span, a11_hat, a12_hat, a21_hat, a22_hat, a11, a12, a21, a22);
fig1 = printer_params_estV2(t_span, b1_hat, b2_hat, b1, b2);


% h eksodos y tou systhmatos kai tou montelou
fig2 = printer_outputs(t_span, y1, y1_hat);
fig5 = printer_outputs(t_span, y2, y2_hat);

% posostiaio sfalma
fig3 = printer_error(t_span, percentage_error, false);
                          
% meso tetragwniko sfalma
fig4 = printer_error(t_span, mean_square_error, true);

saveas(fig0, 'prob3_params_est_arrayA.png')
saveas(fig1, 'prob3_params_est.png')
saveas(fig2, 'prob3_outputsY1.png')
saveas(fig5, 'prob3_outputsY2.png')
saveas(fig3, 'prob3_errorPERC.png')
saveas(fig4, 'prob3_errorMSE.png')