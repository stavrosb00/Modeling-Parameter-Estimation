% B) erwthma: me8odo MH PRAGMATIKOU(OFFLINE) xronou gia ton prosdiorismo 
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
AIC_stats = cell(5,5);
u = @(t)cos(t)+0.1*sin(2*t);
y = sys(t,u);

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
% ypologismoi kai apotelesmata(prints/plots/syntelestes aksiologhshs C erwthma)
i = 1;
fig1 = figure('Name','Offline Estimation Simulation','NumberTitle','off');
fig1.WindowState = 'maximized';

for n = 1:max_dim
    for m = 1:n-1
        y_est{n,m} = least_squares_calc_params(y,n,m,u,u_test, t,pole);
        error{n,m} = y_test - y_est{n,m};
        AIC_stats{n,m} = calculateAIC(n, m, error{n,m}, N, r);
%         data_AIC(i - 1, M) = [AIC_stats{n,m} ; n ; m]; 
        e_max = max(abs(error{n,m}));
        fprintf('Maximum absolute error for n = %d and m = %d is: %d .\n',n,m,e_max)
        
        subplot(4,3,i)
        if abs(e_max) <= 10^(-2) % threshold
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

saveas(fig1, 'offlineEstimation.png')
% data_AIC
% data_AIC = sort(data_AIC, )
% minAIC = 
% errtest = find(~cellfun(@isempty,error));
% errxd = error(errtest);

%% boh8htikes synarthseis
function y_est = least_squares_calc_params(y,n,m,u,u_test, t,pole)
    prev = 1;
    filtered_delta_output = {}; 
    filtered_delta_input = {}; 
    zeta = zeros(length(t),n+m+1);
    numerator = zeros(1,m+1); 
    denominator = zeros(1,n+1);
    % ypologismos tou eusta8ous filtrou L
    for i = 1:n
        filter_L = conv([1 pole],prev);
        prev = filter_L;
    end

    for i = m:-1:0
        prev = 1;
        for j = i:-1:1
            current = conv([1 0], prev);
            prev = current;
        end
        filtered_delta_input{m-i+1} = tf(current, filter_L);
    end
    filtered_delta_input{m+1} = tf(1, filter_L);

    for i = n-1:-1:0
        prev = 1;
        for j = i:-1:1
            current = conv([1 0],prev);
            prev = current;
        end
        filtered_delta_output{n-1-i+1} = tf(current, filter_L);
    end
    filtered_delta_output{n-1+1} = tf(1, filter_L);

    for i = 1:n
        zeta(:,i) = lsim(-filtered_delta_output{i}, y, t);
    end
    for i = 1:m+1
        zeta(:,n+i) = lsim(filtered_delta_input{i}, u(t), t);
    end
    % eksiswsh elaxistwn tetragwnwn
    theta =(y' * zeta) / ((zeta' * zeta));
    for i = 1:n
        theta(i) = theta(i) + filter_L(i+1);
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
    y_est = lsim(g_transfer,u_test(t),t);
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