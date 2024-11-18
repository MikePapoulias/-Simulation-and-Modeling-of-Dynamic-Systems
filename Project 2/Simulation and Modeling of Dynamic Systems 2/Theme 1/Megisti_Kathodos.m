format long


% Pragmatikes times twn parametrwn a,b:
a=3; b=0.5;
% Xroniko diastima t:
t = 0:0.1:10;



% ========================================================================
% =                         (a)    erotima                               =
% ========================================================================

u= @(t) 10;
% Bima g methodou klisis:
g=50;
% pole = o polos tou eustathous filtrou
p=3;

% Epilisi tis diaforikis exiswsis, u=10
options = odeset(Refine=100);
[t,diaf_pinakas] = ode45(@(t,diaf_pinakas) odefun(t,diaf_pinakas,a,b,u,g,p),t,[0 0 0 0 0],options);



% Parametroi tou pinaka diaf_pinakas, u=10
x = diaf_pinakas(:,1);
thita_1_hat = diaf_pinakas(:,2);
thita_2_hat = diaf_pinakas(:,3);
f1 = diaf_pinakas(:,4);
f2 = diaf_pinakas(:,5);
x_hat = thita_1_hat.*f1 + thita_2_hat.*f2;
error = x-x_hat;


% Zitoumenes ektypwseis, u=10
fprintf("Ektimiseis twn a kai b gia eisodo u=10:\n")
a_hat = p-thita_1_hat(length(thita_1_hat))
b_hat = thita_2_hat(length(thita_2_hat))


% Boithitikoi pinakes A, B, P_pole
A = zeros(1,length(t));
B = zeros(1,length(t));
P_pole = zeros(1,length(t));
counter=1;
for i=0:0.1:(length(t)-1)/10
    A(counter) = 3;
    B(counter) = 0.5;
    P_pole(counter) = p;
    counter = counter+1;
end



% Apothikeusi grafimatwn
% ---------------- u=10 -------------------
figure("Name", sprintf("Real-estimate_output x, u=10"))
plot(t, x)
hold on
plot(t, x_hat, 'r')
hold on
plot(t, error, 'green')
xlabel("t (sec)")
ylabel("x,   x_h_a_t")
title('Πραγματική και Εκτιμώμενη Έξοδος x, u=10')
legend('x', 'x_h_a_t', 'error = x - x_h_a_t')
saveas(gcf, 'Real-estimate_output x, u=10')

figure("Name", sprintf("Error, u=10"))
plot(t, x - (thita_1_hat.*f1 + thita_2_hat.*f2), 'green')
xlabel("t (sec)")
ylabel("error = x - x_h_a_t")
title('Error, u=10')
legend('error = x - x_h_a_t')
saveas(gcf, 'Error, u=10')

figure("Name", sprintf("Real-estimate_parameter a, u=10"))
% a=3 stathero kathe xroniki stigmi
plot(t, A)
hold on
plot(t, P_pole-thita_1_hat', 'r')
xlabel("t (sec)")
ylabel("a,   a_h_a_t")
title('Πραγματική και Εκτιμώμενη Παράμετρος a, u=10')
legend('a', 'a_h_a_t')
saveas(gcf, 'Real-estimate_parameter a, u=10')

figure("Name", sprintf("Real-estimate_parameter b, u=10"))
% b=0.5 stathero kathe xroniki stigmi
plot(t, B)
hold on
plot(t, thita_2_hat, 'r')
xlabel("t (sec)")
ylabel("b,   b_h_a_t")
title('Πραγματική και Εκτιμώμενη Παράμετρος b, u=10')
legend('b', 'b_h_a_t')
saveas(gcf, 'Real-estimate_parameter b, u=10')




% ========================================================================
% =                         (b)    erotima                               =
% ========================================================================

u= @(t) 10*sin(3*t);
% Bima g methodou klisis:
g=50;
% pole = o polos tou eustathous filtrou
p=2;

% Epilisi tis diaforikis exiswsis, u=10*sin(3*t)
options = odeset(Refine=100);
[t,diaf_pinakas_b] = ode45(@(t,diaf_pinakas_b) odefun(t,diaf_pinakas_b,a,b,u,g,p),t,[0 0 0 0 0],options);


% Parametroi tou pinaka diaf_pinakas, u=10*sin(3*t)
x = diaf_pinakas_b(:,1);
thita_1_hat = diaf_pinakas_b(:,2);
thita_2_hat = diaf_pinakas_b(:,3);
f1 = diaf_pinakas_b(:,4);
f2 = diaf_pinakas_b(:,5);
x_hat = thita_1_hat.*f1 + thita_2_hat.*f2;
error = x-x_hat;


% Zitoumenes ektypwseis, u=10*sin(3*t)
fprintf("\n\nEktimiseis twn a kai b gia eisodo u=10*sin(3*t):\n")
a_hat = p-thita_1_hat(length(thita_1_hat))
b_hat = thita_2_hat(length(thita_2_hat))

% Boithitikos pinakss P_pole_b
P_pole_b = zeros(1,101);
counter=1;
for i=0:0.1:10
    P_pole_b(counter) = p;
    counter = counter+1;
end


% Apothikeusi grafimatwn
% ---------------- u=10*sin(3*t) -------------------
figure("Name", sprintf("Real-estimate_output x, u=10*sin(3*t)"))
plot(t, x)
hold on
plot(t, x_hat, 'r')
hold on
plot(t, error, 'green')
xlabel("t (sec)")
ylabel("x,   x_h_a_t")
title('Πραγματική και Εκτιμώμενη Έξοδος x, u=10*sin(3*t)')
legend('x', 'x_h_a_t', 'error = x - x_h_a_t')
saveas(gcf, 'Real-estimate_output x, u=10sin(3t)')

figure("Name", sprintf("Error, u=10*sin(3*t)"))
plot(t, x - (thita_1_hat.*f1 + thita_2_hat.*f2), 'green')
xlabel("t (sec)")
ylabel("error = x - x_h_a_t")
title('Error, u=10*sin(3*t)')
legend('error = x - x_h_a_t')
saveas(gcf, 'Error, u=10sin(3t)')

figure("Name", sprintf("Real-estimate_parameter a, u=10*sin(3*t)"))
% a=3 stathero kathe xroniki stigmi
plot(t, A)
hold on
plot(t, P_pole_b-thita_1_hat', 'r')
xlabel("t (sec)")
ylabel("a,   a_h_a_t")
title('Πραγματική και Εκτιμώμενη Παράμετρος a, u=10*sin(3*t)')
legend('a', 'a_h_a_t')
saveas(gcf, 'Real-estimate_parameter a, u=10sin(3t)')

figure("Name", sprintf("Real-estimate_parameter b, u=10*sin(3*t)"))
% b=0.5 stathero kathe xroniki stigmi
plot(t, B)
hold on
plot(t, thita_2_hat, 'r')
xlabel("t (sec)")
ylabel("b,   b_h_a_t")
title('Πραγματική και Εκτιμώμενη Παράμετρος b, u=10*sin(3*t)')
legend('b', 'b_h_a_t')
saveas(gcf, 'Real-estimate_parameter b, u=10sin(3t)')









% function odefun 
function dy = odefun(t,diaf_pinakas,a,b,u,g,p)
% =============================
% diaf_pinakas:
% (1) --> x
% (2) --> thita_1_hat
% (3) --> thita_2_hat
% (4) --> f1 tou F pinaka
% (5) --> f2 tou F pinaka
% =============================
error = diaf_pinakas(1) - (diaf_pinakas(2)*diaf_pinakas(4) + diaf_pinakas(3)*diaf_pinakas(5));
% differential equations:
dx = -a*diaf_pinakas(1) + b*u(t);
dthita_1_hat = g*error*diaf_pinakas(4); 
dthita_2_hat = g*error*diaf_pinakas(5); 
df1 = -p*diaf_pinakas(4)+diaf_pinakas(1);
df2 = -p*diaf_pinakas(5)+u(t);

dy = [dx; dthita_1_hat; dthita_2_hat; df1; df2];
end




