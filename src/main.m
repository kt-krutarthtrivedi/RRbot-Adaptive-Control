% Robot Controls - Adaptive Control for RRbot Manipulator
% Author: Krutarth Trivedi | ktrivedi@wpi.edu

clear; 
clc; 
close all;

% physical parameters of the robot
m1 = 1; r1 = 0.45; l1 = 1; I1 = 0.084;
m2 = 1; r2 = 0.45; l2 = 1; I2 = 0.084;
g = 9.81;

% define symbols
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot 'real'

%--- Initial Setup - Cubic Trajectory Generation ------%
syms t;
timeVec = [1; t; t^2; t^3];
timeVec_dot = diff(timeVec);
timeVec_ddot = diff(timeVec_dot);

t0 = 0; 
tf = 10;
theta1_t0 = pi;
theta1_tf = 0;
theta2_t0 = pi/2;
theta2_tf = 0;
theta1_dot_t0 = 0;
theta1_dot_tf = 0;
theta2_dot_t0 = 0;
theta2_dot_tf = 0;

timeMat = [1, t0, t0^2, t0^3;
           0, 1, 2*t0, 3*t0^2;
           1, tf, tf^2, tf^3;
           0, 1, 2*tf, 3*tf^2];

theta1ConfigVec = [theta1_t0; theta1_dot_t0; theta1_tf; theta1_dot_tf];
theta2ConfigVec = [theta2_t0; theta2_dot_t0; theta2_tf; theta2_dot_tf];

theta1CoefficientVec = pinv(timeMat)*theta1ConfigVec;
theta2CoefficientVec = pinv(timeMat)*theta2ConfigVec;

q1 = theta1CoefficientVec' * timeVec;
q2 = theta2CoefficientVec' * timeVec;

q1_dot = theta1CoefficientVec' * timeVec_dot;
q2_dot = theta2CoefficientVec' * timeVec_dot;

q1_ddot = theta1CoefficientVec' * timeVec_ddot;
q2_ddot = theta2CoefficientVec' * timeVec_ddot;

q_desired = [q1;q2];
qdot_desired = [q1_dot;q2_dot];
qddot_desired = [q1_ddot;q2_ddot];
Q_desired = [q_desired; qdot_desired];

%--------- Manipulator Equation Form ---------%
alpha_1 = I1 + I2 + m1*r1^2 + m2*(l1^2 + r2^2);
alpha_2 = m2*l1*r2;
alpha_3 = I2 + m2*r2^2;
alpha_4 = m1*r1 + m2*l1;
alpha_5 = m2*r2;

Mmat = [alpha_1+2*alpha_2* cos(theta2), alpha_3+alpha_2* cos(theta2); alpha_3+alpha_2* cos(theta2), alpha_3];
Cmat = [-alpha_2* sin(theta2)*theta2_dot, -alpha_2* sin(theta2)*(theta1_dot+theta2_dot); alpha_2* sin(theta2)*theta1_dot,0];
Gmat = [-alpha_4*g*sin(theta1)-alpha_5*g*sin(theta1+theta2); -alpha_5*g*sin(theta1+theta2)];

%--------- Linear Parametric Form -----------%
Y = [theta1_ddot, ...
    cos(theta2)*(2*theta1_ddot + theta2_ddot) - 2* sin(theta2)*theta1_dot*theta2_dot - sin(theta2)*theta2_dot^2, ...
    theta2_ddot, ...
    -sin(theta1)*g, ...
    -sin(theta1 + theta2)*g;
    0, ...
    sin(theta2)*theta1_dot^2 + cos(theta2)*theta1_ddot, ...
    theta1_ddot + theta2_ddot, ...
    0, ...
    -sin(theta1+theta2)*g];

alpha = [alpha_1;
        alpha_2;
        alpha_3;
        alpha_4;
        alpha_5];

%----------- Verification --------------%
Tau_MEF = Mmat*[theta1_ddot; theta2_ddot] + Cmat*[theta1_dot; theta2_dot] + Gmat;
Tau_LPE = Y * alpha;

Tau_MEF = subs(Tau_MEF, [theta1, theta2, theta1_dot, theta2_dot, theta1_ddot, theta2_ddot], [1,1,1,1,1,1]);
Tau_LPE = subs(Tau_LPE, [theta1, theta2, theta1_dot, theta2_dot, theta1_ddot, theta2_ddot], [1,1,1,1,1,1]);
res = isequal(Tau_MEF, Tau_LPE); %retuns 1 if equal, 0 otherwise.
if(res)
    disp("Both the forms are equivalent to each other.");
end

% ------- Adaptive Inverse Dynamics Control --------%
A = [0,0,1,0;
    0,0,0,1;
    0,0,0,0;
    0,0,0,0];

B = [0,0;
    0,0;
    1,0;
    0,1];

desiredEigenValues = [-2,-2,-3,-3];
K = place(A, B, desiredEigenValues);

Acl = A - B*K;
Q = eye(4);     
P = lyap(Acl',Q);
gamma = eye(5)*0.03;

%--------------- Virtual control input ----------------------%
jointValues = [theta1; theta2; theta1_dot; theta2_dot];
error = jointValues - Q_desired;
V = -K*error + qddot_desired;

%-------------- The overall control law ---------------%
Tau = Mmat*V + Cmat*[theta1_dot; theta2_dot] + Gmat;

% -------- Simulation and Plotting ---------%

% simulate the system for 10 sec for given control inputs using ODE45
T = 10;
alpha_hat0 = 0.75*alpha;
y0 = [deg2rad(200), deg2rad(125),0 ,0, ...
        alpha_hat0(1), alpha_hat0(2), alpha_hat0(3), alpha_hat0(4), alpha_hat0(5)];
[t,y] = ode45(@ode_rrbot, 0:0.1:T, y0);

%reconstruct the desired trajectories and control input
q1_reconstruct = [];
q2_reconstruct = [];
q1_dot_reconstruct = [];
q2_dot_reconstruct= [];
q1_ddot_reconstruct = [];
q2_ddot_reconstruct= [];
t1 = [];
t2 = [];

for i = 1:size(t)
    q1_reconstruct(i) = double(subs(q1, t(i,:)));
    q2_reconstruct(i) = double(subs(q2, t(i,:)));
    q1_dot_reconstruct(i) = double(subs(q1_dot, t(i,:)));
    q2_dot_reconstruct(i) = double(subs(q2_dot, t(i,:)));
    q1_ddot_reconstruct(i) = double(subs(q1_ddot, t(i,:)));
    q2_ddot_reconstruct(i) = double(subs(q2_ddot, t(i,:)));

    theta1 = y(i,1);
    theta2 = y(i,2);
    theta1_dot = y(i,3);
    theta2_dot = y(i,4);

    jointValues = [theta1; theta2; theta1_dot; theta2_dot];
    Q_desired = [q1_reconstruct(i);q2_reconstruct(i);q1_dot_reconstruct(i);q2_dot_reconstruct(i)];
    error = jointValues - Q_desired;
    qddot_desired = [q1_ddot_reconstruct(i);q2_ddot_reconstruct(i)];
    V = -K*error + qddot_desired;

    alpha1 = y(i,5);
    alpha2 = y(i,6);
    alpha3 = y(i,7);
    alpha4 = y(i,8);
    alpha5 = y(i,9);

    Mmat_hat = [alpha_1+2*alpha_2* cos(theta2), alpha_3+alpha_2* cos(theta2); alpha_3+alpha_2* cos(theta2), alpha_3];
    Cmat_hat = [-alpha_2* sin(theta2)*theta2_dot, -alpha_2* sin(theta2)*(theta1_dot+theta2_dot); alpha_2* sin(theta2)*theta1_dot,0];
    Gmat_hat = [-alpha_4*g*sin(theta1)-alpha_5*g*sin(theta1+theta2); -alpha_5*g*sin(theta1+theta2)];

    Tau = Mmat_hat*V + Cmat_hat*[y(i,3);y(i,4)] + Gmat_hat;
    t1(i) = Tau(1);
    t2(i) = Tau(2);
end

% plot the trajectories
figure('Name','Trajectories', 'NumberTitle','off');
subplot(3,2,1)
plot(t,(y(:,1)),'b');
title('theta1')
xlabel('T');
ylabel('rad');
hold on;
plot(t,q1_reconstruct,'r');
legend({'Current','Desired'})

subplot(3,2,2)
plot(t,(y(:,2)),'b')
title('theta2')
xlabel('T');
ylabel('rad');
hold on;
plot(t,q2_reconstruct,'r');
legend({'Current','Desired'})

subplot(3,2,3)
plot(t,(y(:,3)),'b')
title('theta1-dot')
xlabel('T');
ylabel('rad/s');
hold on;
plot(t,q1_dot_reconstruct,'r');
legend({'Current','Desired'})

subplot(3,2,4);
plot(t,(y(:,4)),'b');
title('theta2-dot')
xlabel('T');
ylabel('rad/s');
hold on;
plot(t,q2_dot_reconstruct,'r');
legend({'Current','Desired'})

subplot(3,2,5);
plot(t,t1,'b');
title('t1')
xlabel('s');
ylabel('Nm');

subplot(3,2,6);
plot(t,t2,'b');
title('t2')
xlabel('s');
ylabel('Nm');

figure('Name','Alpha Hat', 'NumberTitle','off');
subplot(3,2,1)
plot(t,(y(:,5)));
title('alpha1')
xlabel('T');
ylabel('alpha1');

subplot(3,2,2)
plot(t,(y(:,6)));
title('alpha2')
xlabel('T');
ylabel('alpha2');

subplot(3,2,3)
plot(t,(y(:,7)));
title('alpha3')
xlabel('T');
ylabel('alpha3');

subplot(3,2,4)
plot(t,(y(:,8)));
title('alpha4')
xlabel('T');
ylabel('alpha4');

subplot(3,2,5)
plot(t,(y(:,9)));
title('alpha5')
xlabel('T');
ylabel('alpha5');

fprintf('Eigenvalues are selected such that: \n\n');
fprintf('Torque1: %.3f < u1 < %.3f \n\n',min(t1),max(t1));
fprintf('Torque2: %.3f < u2 < %.3f \n\n', min(t2),max(t2));