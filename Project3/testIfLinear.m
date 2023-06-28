% A) erwthma: dokimh tou systhmatos gia arxh epallhlias kai omoiogeneias
% gia na prokypsei an einai grammiko

clear;
clc;
close all;

t = 0 : 0.001 : 20;

%% dokimes gia arxh epallhlias 
% 1o test
u1 = @(t) cos(t);
u2 = @(t) 2 * sin(2 * t);
u12 = @(t) cos(t) + 2 * sin(2 * t);
% u12 = @(t) u1 + u2;

y1 = sys(t, u1);
y2 = sys(t, u2);
sum12 = y1 + y2;
y12 = sys(t, u12);

e1 = y12 - sum12;

fig1 = figure('Name','1st Linear criterion test 1','NumberTitle','off');
subplot(3,1,1)
plot(t,y12,'LineWidth',0.8);
title('Total output','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('$$ y_{12} $$','Interpreter','Latex');
subplot(3,1,2)
plot(t,sum12,'LineWidth',0.8,'Color','red');
leg1 = legend('$$y_{1}+y_{2}$$');
set(leg1,'Interpreter','latex');
title('Sum of single outputs','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
subplot(3,1,3)
plot(t,e1,'LineWidth',0.8,'Color','red');
leg1 = legend('$$e_{1}$$');
set(leg1,'Interpreter','latex');
title('Error $$e_{1}$$' ,'Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');

saveas(fig1, 'linear_1st_crit_1.png')

% 2o test
u3 = @(t) 5 * cos(4 * t);
u4 = @(t) 7 * sin(6 * t);
u34 = @(t) 5 * cos(4 * t) + 7 * sin(6 * t);

y3 = sys(t, u3);
y4 = sys(t, u4);
sum34 = y3 + y4;
y34 = sys(t, u34);

e2 = y34 - sum34; 

fig2 = figure('Name','1st Linear criterion test 2','NumberTitle','off');
subplot(3,1,1)
plot(t,y34,'LineWidth',0.8);
title('Total output','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
ylabel('$$ y_{34} $$','Interpreter','Latex');
subplot(3,1,2)
plot(t,sum34,'LineWidth',0.8,'Color','red');
leg1 = legend('$$y_{3}+y_{4}$$');
set(leg1,'Interpreter','latex');
title('Sum of single outputs','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
subplot(3,1,3)
plot(t,e2,'LineWidth',0.8,'Color','red');
leg1 = legend('$$e_{2}$$');
set(leg1,'Interpreter','latex');
title('Error $$e_{2}$$' ,'Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');

saveas(fig2, 'linear_1st_crit_2.png')

%% dokimh gia arxh omoiogeneias


%u_tot = @(t) u1(t) + u2(t) + u3(t) + u4(t);
u_tot = @(t) 0.5 * cos(2 * t) + 3 * sin(5 * t);


u_hom50 = @(t) 50* (0.5 * cos(2 * t) + 3 * sin(5 * t));
u_hom3 = @(t) 3* (0.5 * cos(2 * t) + 3 * sin(5 * t));
y_tot = sys(t, u_tot);
y_hom50 = sys(t, u_hom50);
y_hom3 = sys(t, u_hom3);

e_hom50 = y_hom50 - 50 * y_tot;
e_hom3 = y_hom3 - 3 * y_tot;
% 1o test x3
fig3 = figure('Name','2nd Linear criterion test 1','NumberTitle','off');
subplot(3,1,1)
plot(t, 3 * y_tot,'LineWidth',0.8);
title('$$y_{total}$$ output x3','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
subplot(3,1,2)
plot(t,y_hom3,'LineWidth',0.8,'Color','red');
title('$$y_{hom3}$$','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
subplot(3,1,3)
plot(t,e_hom3,'LineWidth',0.8);
title('$$e_{hom3}$$','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');

saveas(fig3, 'linear_2nd_crit_1.png')

% 2o test x50
fig4 = figure('Name','2nd Linear criterion test 2','NumberTitle','off');
subplot(3,1,1)
plot(t, 50 * y_tot,'LineWidth',0.8);
title('$$y_{total}$$ output x50','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
subplot(3,1,2)
plot(t,y_hom50,'LineWidth',0.8,'Color','red');
title('$$y_{hom50}$$','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');
subplot(3,1,3)
plot(t,e_hom50,'LineWidth',0.8);
title('$$e_{hom50}$$','Interpreter','Latex');
xlabel('Time (sec)','Interpreter','Latex');

saveas(fig4, 'linear_2nd_crit_2.png')


