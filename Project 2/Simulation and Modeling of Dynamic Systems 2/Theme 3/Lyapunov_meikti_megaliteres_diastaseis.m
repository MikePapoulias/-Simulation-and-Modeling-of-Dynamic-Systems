format long


% Pragmatikes times twn parametrwn a11,a12,a21,a22,b1,b2:
a11=-0.25;  a12=3;  a21=-5;  a22=0;
b1=0.5;     b2=1.5;

% Xroniko diastima t:
t = 0:0.1:20;

% Eisodos
u= @(t) 3.5*sin(7.2*t) + 2*sin(11.7*t);

% Epiloges pinakwn G1, G2 kai parametrou thita_m
G1=[100 0; 0 100]; 
G2=[100 0; 0 100];
thita_m = 4;

%arxikes_sinthikes:
x_0 = [0 ; 0];
x_hat_0 = [0 ; 0];
A_hat_0 = [0 0 ; 0 0];
B_hat_0 = [0 ; 0];


% Epilisi tis diaforikis exiswsis
options = odeset(Refine=100);
[t,diaf] = ode45(@(t,diaf) odefun(t,diaf,a11,a12,a21,a22,b1,b2,G1,G2,u,thita_m), t ,[x_0 x_hat_0 A_hat_0 B_hat_0], options);


% Parametroi tou pinaka diaf
x1 = diaf(:,1);
x2 = diaf(:,2);
x1_hat = diaf(:,3);
x2_hat = diaf(:,4);
a11_hat = diaf(:,5);
a12_hat = diaf(:,6);
a21_hat = diaf(:,7);
a22_hat = diaf(:,8);
b1_hat = diaf(:,9);
b2_hat = diaf(:,10);


% Zitoumenes ektypwseis
fprintf("Ektimiseis twn A kai B:\n")
A_hat = [a11_hat(length(a11_hat)) a12_hat(length(a12_hat)) ; a21_hat(length(a21_hat)) a22_hat(length(a22_hat))]
B_hat = [b1_hat(length(b1_hat)) ; b2_hat(length(b2_hat))]


% Boithitikoi pinakes A11, A12, A21, A22, B1, B2
A11 = zeros(1,length(t));
A12 = zeros(1,length(t));
A21 = zeros(1,length(t));
A22 = zeros(1,length(t));
B1 = zeros(1,length(t));
B2 = zeros(1,length(t));
counter=1;
for i=0:0.1:(length(t)-1)/10
    A11(counter) = a11;
    A12(counter) = a12;
    A21(counter) = a21;
    A22(counter) = a22;
    B1(counter) = b1;
    B2(counter) = b2;
    counter = counter+1;
end




% Apothikeusi grafimatwn
figure("Name", sprintf("Real-estimate_output x1"))
plot(t, x1)
hold on
plot(t, x1_hat)
hold on
plot(t, x1-x1_hat, 'green')
xlabel("t (sec)")
ylabel("x1,   x1_h_a_t")
title('Πραγματική και Εκτιμώμενη Έξοδος x1')
legend('x1', 'x1_h_a_t', 'error = x1 - x1_h_a_t')
saveas(gcf, 'Real-estimate_output x1')

figure("Name", sprintf("Real-estimate_output x2"))
plot(t, x2)
hold on
plot(t, x2_hat)
hold on
plot(t, x2-x2_hat, 'green')
xlabel("t (sec)")
ylabel("x2,   x2_h_a_t")
title('Πραγματική και Εκτιμώμενη Έξοδος x2')
legend('x2', 'x2_h_a_t', 'error = x2 - x2_h_a_t')
saveas(gcf, 'Real-estimate_output x2')

% ================== pinakas A_hat ===================
figure("Name", sprintf("Real-estimate_table A"))
subplot(2,2,1)
plot(t, A11)
hold on
plot(t, a11_hat)
xlabel("t (sec)")
ylabel("a11,   a11_h_a_t")
title('Πραγματική/Εκτιμώμενη a11')
legend('a11', 'a11_h_a_t')

subplot(2,2,2)
plot(t, A12)
hold on
plot(t, a12_hat)
xlabel("t (sec)")
ylabel("a12,   a12_h_a_t")
title('Πραγματική/Εκτιμώμενη a12')
legend('a12', 'a12_h_a_t')

subplot(2,2,3)
plot(t, A21)
hold on
plot(t, a21_hat)
xlabel("t (sec)")
ylabel("a21,   a21_h_a_t")
title('Πραγματική/Εκτιμώμενη a21')
legend('a21', 'a21_h_a_t')

subplot(2,2,4)
plot(t, A22)
hold on
plot(t, a22_hat)
xlabel("t (sec)")
ylabel("a22,   a22_h_a_t")
title('Πραγματική/Εκτιμώμενη a22')
legend('a22', 'a22_h_a_t')

saveas(gcf, 'Real-estimate_table A')

% ================== pinakas B_hat ===================
figure("Name", sprintf("Real-estimate_table B"))
subplot(2,1,1)
plot(t, B1)
hold on
plot(t, b1_hat)
xlabel("t (sec)")
ylabel("b1,   b1_h_a_t")
title('Πραγματική/Εκτιμώμενη b1')
legend('b1', 'b1_h_a_t')

subplot(2,1,2)
plot(t, B2)
hold on
plot(t, b2_hat)
xlabel("t (sec)")
ylabel("b2,   b2_h_a_t")
title('Πραγματική/Εκτιμώμενη b2')
legend('b2', 'b2_h_a_t')

saveas(gcf, 'Real-estimate_table B')









% function odefun 
function dy = odefun(t,diaf,a11,a12,a21,a22,b1,b2,G1,G2,u,thita_m)
% =============================
% diaf_meiktos:
% (1) --> x1      (2) --> x2    
% (3) --> x1_hat  (4) --> x2_hat  
% (5) --> a11_hat (6) --> a12_hat (7) --> a21_hat (8) --> a22_hat
% (9) --> b1_hat  (10) -> b2_hat
% =============================
A = [a11 a12; a21 a22];
B = [b1 ; b2];

x = [diaf(1); diaf(2)];
x_hat = [diaf(3); diaf(4)];
A_hat = [diaf(5) diaf(6); diaf(7) diaf(8)];
B_hat = [diaf(9); diaf(10)];

% sfalma e:
e = x - x_hat;
% differential equations:
dx = A*x + B*u(t);
dx_hat = A_hat*x + B_hat*u(t) + thita_m*e;
dA_hat = G1*e*x';
dB_hat = G2*e*u(t)';

dy = [dx(1); dx(2); dx_hat(1); dx_hat(2); dA_hat(1,1); dA_hat(1,2); dA_hat(2,1); dA_hat(2,2); dB_hat(1); dB_hat(2)];
end

