%% Thema 2 i) : me8odo Lyapunov me parallhlh domh gia online ektimhsh parametrwn

clear;
clc;
close all;


%% pragmatiko systhma

% eisodos u kai parametroi methodou
u = @(t) 10 * sin(3 * t);
a = 3;
b = 0.5;

gamma1 = 4;
gamma2 = 1;

% arxikh syn8hkh
x0 = [0, 0, 0, 0]';
t_span = 0:0.01:40;

%% pros8hkh 8orybou sto shma x 

% epilogh typou 8orybou ana peirama
% 1-arxikh official
n0 = 0.5;
f = 40;
% PEIRAMATA
% 2
% n0 = 1;
% f = 40;
% 3
% n0 = 0.15;
% f = 40;
% 4
% n0 = 0.5;
% f = 20;
% 5
% n0 = 0.5;
% f = 300;

n = @(t) n0 * sin(2 * pi * f * t);


%% prosomoiwsh kai ektimhsh parametrwn 

% lysh diaforikou systhmatos gia xroniko diasthma t_span
[t, x] = ode15s(@(t, x) system_equationsV2(t, x, a, b, gamma1, gamma2, u, n), t_span, x0);

y = x(:, 1);
y_hat = x(:, 2);
theta_hat1 = x(:, 3);
theta_hat2 = x(:, 4);

% telikes ektimwmenes times parametrwn a, b
a_hat = theta_hat1(:);
b_hat = theta_hat2(:);

a_hat(length(t_span))
b_hat(length(t_span))


%% aksiologhsh basei metrikwn sfalmatos

% Percentage error
percentage_error = zeros(length(t_span), 1);
for i = 1:length(t_span)
    percentage_error(i) = abs((y(i) - y_hat(i)) / y(i));
end

% Mean square error
mean_square_error = zeros(length(t_span), 1);
for i = 1:length(t_span)
    mean_square_error(i) = (y(i) - y_hat(i)) ^ 2;
end


%% grafikes parastaseis me plots

% kata poso sygklinoyn oi ektimwmenes parametroi tou montelou me to systhma
fig1 = printer_params_est(t_span, a_hat, b_hat, a, b);

% h eksodos y tou systhmatos kai tou montelou
fig2 = printer_outputs(t_span, y, y_hat);
                                                                      
% posostiaio sfalma
fig3 = printer_error(t_span, percentage_error, false);
                          
% meso tetragwniko sfalma
fig4 = printer_error(t_span, mean_square_error, true);

saveas(fig1, 'prob2_i_params_est.png')
saveas(fig2, 'prob2_i_outputs.png')
saveas(fig3, 'prob2_i_errorPERC.png')
saveas(fig4, 'prob2_i_errorMSE.png')