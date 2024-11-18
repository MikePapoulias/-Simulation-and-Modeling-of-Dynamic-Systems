format long

% a = o polos tou eustathous filtrou
a=100;


% 1000001 time_slots
time = zeros(1,1000001);
counter=1;
for i = 0:0.00001:10
    time(counter) = i;
    counter = counter+1;
end

% Klisi tou arxeiou v.p gia antlisi twn metrisewn exodou
[VR, VC] = v(time);


% VR_vector, VC_vector, U1_vector, U2_vector
VR_vector = VR;
VC_vector = VC;
U1_vector = zeros(1,1000001);
U2_vector = zeros(1,1000001);
counter = 1;
for i = 0:0.00001:10
    U1_vector(counter) = 2*sin(4*i);
    U2_vector(counter) = 4;
    counter = counter+1;
end


% Dimiourgia twn stilwn toy pinaka F_VR
t = 0:0.00001:10;  % 1000001 points

sys_VR_1 = tf(-[1 0],[1 2*a a^2]);
f1_VR = lsim(sys_VR_1,VR_vector,t);

sys_VR_2 = tf(-1, [1 2*a a^2]);
f2_VR = lsim(sys_VR_2,VR_vector,t);

sys_VR_3 = tf(1, [1 2*a a^2]);
f3_VR = lsim(sys_VR_3,U1_vector,t);

sys_VR_4 = tf([1 0 0], [1 2*a a^2]);
f4_VR = lsim(sys_VR_4,U1_vector,t);

sys_VR_5 = tf([1 0 0], [1 2*a a^2]);
f5_VR = lsim(sys_VR_5,U2_vector,t);

% Dimiourgia twn stilwn toy pinaka F_VC
t = 0:0.00001:10;  % 1000001 points

sys_VC_1 = tf(-[1 0],[1 2*a a^2]);
f1_VC = lsim(sys_VC_1,VC_vector,t);

sys_VC_2 = tf(-1, [1 2*a a^2]);
f2_VC = lsim(sys_VC_2,VC_vector,t);

sys_VC_3 = tf(1, [1 2*a a^2]);
f3_VC = lsim(sys_VC_3,U2_vector,t);

sys_VC_4 = tf([1 0], [1 2*a a^2]);
f4_VC = lsim(sys_VC_4,U1_vector,t);

sys_VC_5 = tf(-[1 0], [1 2*a a^2]);
f5_VC = lsim(sys_VC_5,U2_vector,t);


% Dimiourgia tou F_VR pinaka
F_VR = [f1_VR, f2_VR, f3_VR, f4_VR, f5_VR];

% Dimiourgia tou F_VC pinaka
F_VC = [f1_VC, f2_VC, f3_VC, f4_VC, f5_VC];



% Klisi tis methodou elaxistwn tetragwnwn gia to thita_VR 
thita_VR = methodos_elaxistwn_tetragwnwn(VR_vector, F_VR);
fprintf("thita_VR:")
vpa(thita_VR)

% Klisi tis methodou elaxistwn tetragwnwn gia to thita_VC 
thita_VC = methodos_elaxistwn_tetragwnwn(VC_vector, F_VC);
fprintf("thita_VC:")
vpa(thita_VC)


% Dimiourgia VR pinaka ektimisis
VR_vector_estimation = thita_VR * F_VR';

% Dimiourgia VC pinaka ektimisis
VC_vector_estimation = thita_VC * F_VC';



%----------------------------------------------------------------------------
% Dimiourgia VR_noise pinaka
VR_noise = VR;
VR_noise(200000) = VR_noise(200000) + 1500;
VR_noise(450000) = VR_noise(450000) + 820;
VR_noise(800000) = VR_noise(800000) + 1250;

% Dimiourgia VC_noise pinaka
VC_noise = VC;
VC_noise(200000) = VC_noise(200000) + 1500;
VC_noise(450000) = VC_noise(450000) + 820;
VC_noise(800000) = VC_noise(800000) + 1250;

%----------------------------------------------------------------------------
% Dimiourgia twn stilwn toy pinaka F_VR_noise
t = 0:0.00001:10;  % 1000001 points

sys_VR_1_noise = tf(-[1 0],[1 2*a a^2]);
f1_VR_noise = lsim(sys_VR_1_noise,VR_noise,t);

sys_VR_2_noise = tf(-1, [1 2*a a^2]);
f2_VR_noise = lsim(sys_VR_2_noise,VR_noise,t);

sys_VR_3_noise = tf(1, [1 2*a a^2]);
f3_VR_noise = lsim(sys_VR_3_noise,U1_vector,t);

sys_VR_4_noise = tf([1 0 0], [1 2*a a^2]);
f4_VR_noise = lsim(sys_VR_4_noise,U1_vector,t);

sys_VR_5_noise = tf([1 0 0], [1 2*a a^2]);
f5_VR_noise = lsim(sys_VR_5_noise,U2_vector,t);

% Dimiourgia tou F_VR_noise pinaka
F_VR_noise = [f1_VR_noise, f2_VR_noise, f3_VR_noise, f4_VR_noise, f5_VR_noise];

% Klisi tis methodou elaxistwn tetragwnwn gia to thita_VR_noise 
thita_VR_noise = methodos_elaxistwn_tetragwnwn(VR_noise, F_VR_noise);
fprintf("\n\nthita_VR_noise:")
vpa(thita_VR_noise)

%----------------------------------------------------------------------------------
% Dimiourgia twn stilwn toy pinaka F_VC_noise
t = 0:0.00001:10;  % 1000001 points

sys_VC_1_noise = tf(-[1 0],[1 2*a a^2]);
f1_VC_noise = lsim(sys_VC_1_noise,VC_noise,t);

sys_VC_2_noise = tf(-1, [1 2*a a^2]);
f2_VC_noise = lsim(sys_VC_2_noise,VC_noise,t);

sys_VC_3_noise = tf(1, [1 2*a a^2]);
f3_VC_noise = lsim(sys_VC_3_noise,U2_vector,t);

sys_VC_4_noise = tf([1 0], [1 2*a a^2]);
f4_VC_noise = lsim(sys_VC_4_noise,U1_vector,t);

sys_VC_5_noise = tf(-[1 0], [1 2*a a^2]);
f5_VC_noise = lsim(sys_VC_5_noise,U2_vector,t);


% Dimiourgia tou F_VC_noise pinaka
F_VC_noise = [f1_VC_noise, f2_VC_noise, f3_VC_noise, f4_VC_noise, f5_VC_noise];

% Klisi tis methodou elaxistwn tetragwnwn gia to thita_VC_noise 
thita_VC_noise = methodos_elaxistwn_tetragwnwn(VC_noise, F_VC_noise);
fprintf("\n\nthita_VC_noise:")
vpa(thita_VC_noise)

%----------------------------------------------------------------------------------

% Apothikeusi grafimatwn
figure("Name", sprintf("VR Real-estimate_output"))
plot(time, VR_vector)
hold on
plot(time, VR_vector_estimation)
xlabel("t (sec)")
ylabel("VR,   VR_h_a_t")
title('VR Πραγματική και Εκτιμώμενη Έξοδος')
legend('VR', 'VR_h_a_t')
saveas(gcf, 'VR Real-estimate_output')

figure("Name", sprintf("VC Real-estimate_output"))
plot(time, VC_vector)
hold on
plot(time, VC_vector_estimation)
xlabel("t (sec)")
ylabel("VC,   VC_h_a_t")
title('VC Πραγματική και Εκτιμώμενη Έξοδος')
legend('VC', 'VC_h_a_t')
saveas(gcf, 'VC Real-estimate_output')


figure("Name", sprintf("VR_Error"))
plot(time, VR_vector_estimation-VR_vector, 'green')
axis([0 10 -0.00001 0.00001])
xlabel("t (time)")
ylabel("error = VR_h_a_t - VR")
title('VR Σφάλμα Εξόδου')
saveas(gcf, 'VR_Error')

figure("Name", sprintf("VC_Error"))
plot(time, VC_vector_estimation-VC_vector, 'green')
xlabel("t (time)")
ylabel("error = VC_h_a_t - VC")
title('VC Σφάλμα Εξόδου')
saveas(gcf, 'VC_Error')


figure("Name", sprintf("VR_noise"))
plot(time, VR_noise)
xlabel("t (time)")
ylabel("VR_n_o_i_s_e")
title('VR Θόρυβος')
saveas(gcf, 'VR_noise')

figure("Name", sprintf("VC_noise"))
plot(time, VC_noise)
xlabel("t (time)")
ylabel("VC_n_o_i_s_e")
title('VC Θόρυβος')
saveas(gcf, 'VC_noise')




% function methodos elaxistwn tetragwnwn
function [thita] = methodos_elaxistwn_tetragwnwn(Y, F)
    thita = Y*F*inv(F'*F);
end
