%% Thema 1 : me8odo klishs gia online ektimhsh parametrwn

clear;
clc;
close all;

%% pragmatiko systhma
question = true; %gia a erwthma
% question = false; %gia b erwthma

% eisodos u kai parametroi methodou
if question
    u = @(t) 10;
    am = 5;
    gamma = 5;
else
    u = @(t) 10 * sin(3 * t);
    am = 6;
    gamma = 15; 
end

a = 3;
b = 0.5;

% arxikh syn8hkh
x0 = [0, 0, 0, 0, 0, 0]';
t_span = 0:0.01:40;

%% prosomoiwsh kai ektimhsh parametrwn 

% lysh diaforikou systhmatos gia xroniko diasthma t_span
[t, x] = ode15s(@(t, x) system_equationsV1(t, x, a, b, am, gamma, u), t_span, x0);

y = x(:, 1);
theta_hat1 = x(:, 2);
theta_hat2 = x(:, 3);
phi1 = x(:, 4);
phi2 = x(:, 5);
y_hat = x(:, 6);

% telikes ektimwmenes times parametrwn a, b
a_hat = am - theta_hat1(:);
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

if question
    saveas(fig1, 'prob1_a_params_est.png')
    saveas(fig2, 'prob1_a_outputs.png')
    saveas(fig3, 'prob1_a_errorPERC.png')
    saveas(fig4, 'prob1_a_errorMSE.png')
else
    saveas(fig1, 'prob1_b_params_est.png')
    saveas(fig2, 'prob1_b_outputs.png')
    saveas(fig3, 'prob1_b_errorPERC.png')
    saveas(fig4, 'prob1_b_errorMSE.png')
end

