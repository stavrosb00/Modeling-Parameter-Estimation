% B) erwthma: me8odo PRAGMATIKOU(ONLINE) xronou gia ton prosdiorismo 
% twn parametrwn tou montelou

clear;
clc;
close all;

dt = 0.001;
t_end = 20;
t = 0:dt:t_end;

N=t_end/dt+1;
error = cell(5,5); 
y_est = cell(5,5);
x_dyn = cell(5,5);
AIC_stats = cell(5,5);
% u = @(t)cos(4*t)+0.1*sin(6*t);
u = @(t)cos(t)+0.1*sin(2*t);
y = sys(t,u);

% u_test = @(t)cos(4*t)+0.1*sin(6*t); %+cos(16*t)+sin(20*t);
u_test = @(t)cos(t)+0.1*sin(2*t)+0.5*cos(3*t)+0.5*sin(5*t);
y_test = sys(t,u_test);
max_dim = 5;
% M = 0;
% for n = 1:max_dim
%     for m = 1:n-1
%         M = M + 1;
%     end
% end
% data_AIC = zeros(3,M); 

r = 2; %typiko gia AIC
pole = 4;
b = 0.6;
% ypologismoi kai apotelesmata(prints/plots/syntelestes aksiologhshs C erwthma)
i = 1;
fig1 = figure('Name','Online Estimation Simulation','NumberTitle','off');
fig1.WindowState = 'maximized';
% (y ,n, m, u, u_test, b, dt,t_s ,pole)
for n = 1:max_dim
    for m = 1:n-1
        y_est{n,m} = least_squares_calc_params(y, n, m, u, u_test, b, dt, t, pole);
        error{n,m} = y_test - y_est{n,m};
        AIC_stats{n,m} = calculateAIC(n, m, error{n,m}, N, r);
        e_max = max(abs(error{n,m}));
        
        fprintf('Maximum absolute error for n = %d and m = %d is: %d .\n',n,m,e_max)
        
        subplot(4,3,i)
        if abs(e_max) <= 5*10^(-1) % threshold
            plot(t,error{n,m},'LineWidth',0.8,'Color','magenta');
    
        else
            plot(t,error{n,m},'LineWidth',0.8);
        end
        title(['Error for n = ' num2str(n) ' and m = ' num2str(m) ' AIC = ' num2str(AIC_stats{n,m})],'Interpreter','Latex');
        ylabel('$$ y - \hat{y} $$','Interpreter','Latex');
        xlabel('Time (sec)','Interpreter','Latex');
        
        i = i + 1;
    end
end

saveas(fig1, 'onlineEstimation.png')

%% boh8htikes synarthseis
function y_est = least_squares_calc_params(y ,n, m, u, u_test, b, dt,t_s ,pole)
%b~ syntelesths para8yrou apwleias mnhmhs
global A_y B_y C_y D_y; 
global A_u B_u C_u D_u; 
global P_size;
Q0 = eye(n+m+1);
P0 = inv(Q0);
prev = 1;
numerator = zeros(1,m+1); 
denominator = zeros(1,n+1); 
theta_0 = zeros(1, n+m+1);
total_size = n+m+1;
% ypologismos tou eusta8ous filtrou L
for i = 1:n
    L = conv([1 pole],prev);
    prev = L;
end
% ypologismos twn mhdenikwn arxikwn syn8hkwn katastashs x0
for i = 1:n+m+1
    x0(i) = 0;
end
temp = reshape(P0, 1, []);
x0 = [x0 temp];
x0 = [x0 theta_0];
% symbolopoiw se morfh ss tis filtrarismenes eksodous kai eisodous gia na tis xeiristw sthn ODE
H_1 = zeros(n);
for i = 1:n
    for j = 1:n
        if (i == j)
              H_1(i,j) = -1;
        end
    end
end
[A_y, B_y, C_y, D_y] = tf2ss(H_1,L);

H_2 = zeros(m+1);
for c_i = 1:m+1
    for c_j = 1:m+1
        if (c_i == c_j)
            H_2(c_i,c_j) = 1;
        end
    end
end
[A_u, B_u, C_u, D_u] = tf2ss(H_2,L);
A_u = A_u(1:m+1,1:m+1);
C_u = C_u(1:m+1,1:m+1);
% lysh diaforikou systhmatos gia xroniko diasthma t_s
[t, x] = ode45(@(t,x) system_equations(t, x, y, n, m, u, b, dt, t_s),t_s,x0);

% antlhsh twn xrhsimwn plhroforiwn
theta = x(end,total_size+P_size+1:end);
for i = 1:n
    theta(i) = theta(i) + L(i+1);
end
% ektimwmeno y
for i = 1:m+1
    numerator(i) = theta(n+i);
end
denominator(1) = 1;
for i = 1:n
    denominator(i+1) = theta(i);
end
g_transfer = tf(numerator,denominator);
y_est = lsim(g_transfer,u_test(t_s),t_s);
end

function system = system_equations(t, x, y, n, m, u, beta, dt, t_s) 
    global A_y B_y C_y D_y; 
    global A_u B_u C_u D_u; 
    global P_size theta_size;
    P_size = (n + m + 1) * (n + m + 1);
    theta_size = (n + m + 1);
    total_size = (n + m + 1);
    y_now = interp1(t_s,y,t);
    for i = 1:n       %n-taksh shmatos eksodou
        dxydt(i) = A_y(i,:) * x(1:n) + B_y(i) * y_now; %ka8e seira tou A_y gia ta prwta n stoixeia tou x
    end
    phiy = C_y * x(1:n) + D_y * y_now;
    for i = 1:m+1       %m-taksh shmatos eisodou
        dxudt(i) = A_u(i,:) * x(n+1:total_size) + B_u(i) * u(t); %ka8e seira tou A_u gia ta stoixeia apo n+1 mexri ta epomena m+1 tou x
    end
    phiu = C_u * x(n+1:total_size) + D_u * u(t);

    phi = [phiy; phiu];
    P = reshape(x(total_size+1:total_size+1+P_size-1),total_size,[])';
    dP = beta * P - P * (phi * phi') * P ;
    dP = reshape(dP,1, []);
    theta = x(total_size+1+P_size:total_size+1+P_size+theta_size-1); 
    dtheta_dt = P * (y_now - phi' * theta) * phi;
    %% h eksodos tou systhmatos eksiswsewn gia thn ODE
    for i = 1:n
        system(i,1) = dxydt(i);
    end
    for i = 1:m+1
        system(n+i,1) = dxudt(i);
    end
    for i = 1:P_size
        system(total_size+i,1) = dP(i);
    end
    for i = 1:theta_size
        system(total_size+P_size+i,1) = dtheta_dt(i);
    end
end

function AIC_val = calculateAIC(n, m, error, N, r)
    %ypologismos Akaine Information Criterion 
    temp=0;
    for i=1:floor(N)
    temp=temp+error(i)^2;
    end
    I=temp/N;
    r=2;
    k=m+n+1; 
    AIC_val = N*log(I) + r*k;
end