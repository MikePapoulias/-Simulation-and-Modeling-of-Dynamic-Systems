format long

% a = o polos tou eustathous filtrou
a=2;


% Epilisi tis diagorikis exiswsis
options = odeset(Refine=100);
[t,y_pragm] = ode45(@odefun,[0:0.1:10],[0; 0],options);



% Y_vector, U_vector
U_vector = zeros(1,101);
Y_vector = y_pragm(:,1);
counter = 1;
for i=0:0.1:10
    U_vector(counter) = 15*sin(3*i) + 8;
    counter = counter+1;
end


% Dimiourgia twn stilwn toy pinaka F
t =  0:0.1:10;  % 101 points

sys_1 = tf(-[1 0],[1 2*a a^2]);
f1 = lsim(sys_1,Y_vector,t);

sys_2 = tf(-1, [1 2*a a^2]);
f2 = lsim(sys_2,Y_vector,t);

sys_3 = tf(1, [1 2*a a^2]);
f3 = lsim(sys_3,U_vector,t);


% Dimiourgia tou F pinaka
F = [f1, f2, f3];


% Klisi tis methodou elaxistwn tetragwnwn 
thita = methodos_elaxistwn_tetragwnwn(Y_vector, F)


% Dimiourgia Y_vector_estimation pinaka ektimisis
Y_vector_estimation = thita * F';



% Dimiourgia pinaka time gia ton prosdiorismo tou x axona
time = zeros(1,101);
counter=1;
for i = 0:0.1:10
    time(counter) = i;
    counter = counter+1;
end

% Apothikeusi grafimatwn
figure("Name", sprintf("Real-estimate_output"))
plot(time, Y_vector)
hold on
plot(time, Y_vector_estimation, 'o')
hold on
plot(time, Y_vector_estimation-Y_vector', 'green')
xlabel("t (sec)")
ylabel("y,   y_h_a_t")
title('Πραγματική και Εκτιμώμενη Έξοδος')
legend('y', 'y_h_a_t', 'error = y_h_a_t - y')
saveas(gcf, 'Real-estimate_output')

figure("Name", sprintf("Error"))
plot(time, Y_vector_estimation-Y_vector', 'green')
xlabel("t (sec)")
ylabel("error = y_h_a_t - y")
title('Σφάλμα Εξόδου')
saveas(gcf, 'Error')







% function odefun 
function dy = odefun(t,x)
m=10;
b=0.5;
k=2.5;
u = 15*sin(3*t) + 8;
% differential equation:
dy = [x(2);-(b/m)*x(2)-(k/m)*x(1)+(1/m)*u];
end






% function methodos elaxistwn tetragwnwn
function [thita] = methodos_elaxistwn_tetragwnwn(Y, F)
    thita = Y'*F*inv(F'*F);
end


