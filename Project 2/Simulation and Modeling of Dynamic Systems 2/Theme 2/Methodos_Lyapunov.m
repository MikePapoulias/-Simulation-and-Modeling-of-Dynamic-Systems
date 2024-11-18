format long


% Pragmatikes times twn parametrwn a,b:
a=3; b=0.5;
% Xroniko diastima t:
t = 0:0.1:20;

% times parametrwn thoribou h0, f:
h0=0.5;
f=40;


% ========================================================================
% =                          (i) parallili                               =
% ========================================================================

u= @(t) 10*sin(3*t);
% Epiloges parametrwn g1, g2 methodou Lyapunov-parallilis domis:
g1=20; 
g2=1;

% Epilisi tis diaforikis exiswsis erwtimatos (i)
options = odeset(Refine=100);
[t,diaf_parallilos] = ode45(@(t,diaf_parallilos) odefun_parallili(t,diaf_parallilos,a,b,u,g1,g2,h0,f),t,[0 0 0 0],options);


% Parametroi tou pinaka diaf_parallilos
x = diaf_parallilos(:,1);
x_hat = diaf_parallilos(:,2);
A_hat = diaf_parallilos(:,3);
B_hat = diaf_parallilos(:,4);
error = x - x_hat;


% Zitoumenes ektypwseis parallilis domis
fprintf("Ektimiseis twn a kai b parallilis domis\n")
a_hat = A_hat(length(A_hat)) 
b_hat = B_hat(length(B_hat))



% Boithitikoi pinakes A, B
A = zeros(1,length(t));
B = zeros(1,length(t));
counter=1;
for i=0:0.1:(length(t)-1)/10
    A(counter) = 3;
    B(counter) = 0.5;
    counter = counter+1;
end




% Apothikeusi grafimatwn
% ------------- parallili domi ----------------
figure("Name", sprintf("Real-estimate_output x, parallili"))
plot(t, x)
hold on
plot(t, x_hat)
hold on
plot(t, error, 'green')
xlabel("t (sec)")
ylabel("x,   x_h_a_t")
title('Πραγματική και Εκτιμώμενη Έξοδος x, parallili')
legend('x', 'x_h_a_t', 'error = x - x_h_a_t')
saveas(gcf, 'Real-estimate_output x, parallili')

figure("Name", sprintf("Error, parallili"))
plot(t, error, 'green')
xlabel("t (sec)")
ylabel("error = x - x_h_a_t")
title('Error, parallili')
legend('error = x - x_h_a_t')
saveas(gcf, 'Error, parallili')

figure("Name", sprintf("Real-estimate_parameter a, parallili"))
% a=3 stathero kathe xroniki stigmi
plot(t, A)
hold on
plot(t, A_hat)
xlabel("t (sec)")
ylabel("a,   a_h_a_t")
title('Πραγματική και Εκτιμώμενη Παράμετρος a, parallili')
legend('a', 'a_h_a_t')
saveas(gcf, 'Real-estimate_parameter a, parallili')

figure("Name", sprintf("Real-estimate_parameter b, parallili"))
% b=0.5 stathero kathe xroniki stigmi
plot(t, B)
hold on
plot(t, B_hat)
xlabel("t (sec)")
ylabel("b,   b_h_a_t")
title('Πραγματική και Εκτιμώμενη Παράμετρος b, parallili')
legend('b', 'b_h_a_t')
saveas(gcf, 'Real-estimate_parameter b, parallili')







% ========================================================================
% =                            (ii) meikti                               =
% ========================================================================

u= @(t) 10*sin(3*t);
% Epiloges parametrwn g1, g2,thita_m methodou Lyapunov-meiktis domis:
g1=20; 
g2=1;
thita_m = 4;

% Epilisi tis diaforikis exiswsis erwtimatos (ii)
options = odeset(Refine=100);
[t,diaf_meiktos] = ode45(@(t,diaf_meiktos) odefun_meikti(t,diaf_meiktos,a,b,u,g1,g2,thita_m,h0,f),t,[0 0 0 0],options);


% Parametroi tou pinaka diaf_meiktos
x = diaf_meiktos(:,1);
x_hat = diaf_meiktos(:,2);
A_hat = diaf_meiktos(:,3);
B_hat = diaf_meiktos(:,4);
error = x - x_hat;


% Zitoumenes ektypwseis meiktis domis
fprintf("\nEktimiseis twn a kai b meiktis domis\n")
a_hat = A_hat(length(A_hat)) 
b_hat = B_hat(length(B_hat))



% Apothikeusi grafimatwn
% -------------- meikti domi -----------------
figure("Name", sprintf("Real-estimate_output x, meikti"))
plot(t, x)
hold on
plot(t, x_hat)
hold on
plot(t, error, 'green')
xlabel("t (sec)")
ylabel("x,   x_h_a_t")
title('Πραγματική και Εκτιμώμενη Έξοδος x, meikti')
legend('x', 'x_h_a_t', 'error = x - x_h_a_t')
saveas(gcf, 'Real-estimate_output x, meikti')

figure("Name", sprintf("Error, meikti"))
plot(t, error, 'green')
xlabel("t (sec)")
ylabel("error = x - x_h_a_t")
title('Error, meikti')
legend('error = x - x_h_a_t')
saveas(gcf, 'Error, meikti')

figure("Name", sprintf("Real-estimate_parameter a, meikti"))
% a=3 stathero kathe xroniki stigmi
plot(t, A)
hold on
plot(t, A_hat)
xlabel("t (sec)")
ylabel("a,   a_h_a_t")
title('Πραγματική και Εκτιμώμενη Παράμετρος a, meikti')
legend('a', 'a_h_a_t')
saveas(gcf, 'Real-estimate_parameter a, meikti')

figure("Name", sprintf("Real-estimate_parameter b, meikti"))
% b=0.5 stathero kathe xroniki stigmi
plot(t, B)
hold on
plot(t, B_hat)
xlabel("t (sec)")
ylabel("b,   b_h_a_t")
title('Πραγματική και Εκτιμώμενη Παράμετρος b, meikti')
legend('b', 'b_h_a_t')
saveas(gcf, 'Real-estimate_parameter b, meikti')










% function odefun 
function dy = odefun_parallili(t,diaf_parallilos,a,b,u,g1,g2,h0,f)
% =============================
% diaf_parallilos:
% (1) --> x
% (2) --> x_hat
% (3) --> a_hat
% (4) --> b_hat
% =============================
% thoribos:
noise = h0*sin(2*pi*f*t);
% sinoliko metroumeno sfalma:
error = diaf_parallilos(1)-diaf_parallilos(2)+noise;
% differential equations:
dx = -a*diaf_parallilos(1) + b*u(t);
dx_hat = -diaf_parallilos(3)*diaf_parallilos(2) + diaf_parallilos(4)*u(t);
da_hat = -g1*error*diaf_parallilos(2);
db_hat = g2*error*u(t);

dy = [dx; dx_hat; da_hat; db_hat];
end






% function odefun 
function dy = odefun_meikti(t,diaf_meiktos,a,b,u,g1,g2,thita_m,h0,f)
% =============================
% diaf_meiktos:
% (1) --> x
% (2) --> x_hat
% (3) --> a_hat
% (4) --> b_hat
% =============================
% thoribos:
noise = h0*sin(2*pi*f*t);
% sinoliko metroumeno sfalma:
error = diaf_meiktos(1)-diaf_meiktos(2)+noise;
% differential equations:
dx = -a*diaf_meiktos(1) + b*u(t);
dx_hat = -diaf_meiktos(3)*(diaf_meiktos(1)+h0*sin(2*pi*f*t)) + diaf_meiktos(4)*u(t) + thita_m*error;
da_hat = -g1*error*(diaf_meiktos(1)+h0*sin(2*pi*f*t));
db_hat = g2*error*u(t);

dy = [dx; dx_hat; da_hat; db_hat];
end

