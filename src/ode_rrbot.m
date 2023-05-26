% Robot Controls - Adaptive Control for RRbot Manipulator
% Author: Krutarth Trivedi | ktrivedi@wpi.edu

function dX = ode_rrbot(t,X)

% physical parameters of the robot
m1 = 1; r1 = 0.45; l1 = 1; I1 = 0.084;
m2 = 1; r2 = 0.45; l2 = 1; I2 = 0.084;
g = 9.81;

%--------------- Get the current state values ---------------%
dX = zeros(9,1);
X = num2cell(X);
[theta1, theta2, theta1_dot, theta2_dot,...
    alpha_1_hat, alpha_2_hat, alpha_3_hat, alpha_4_hat, alpha_5_hat] = deal(X{:});
jointValues = [theta1; theta2; theta1_dot; theta2_dot];

%--- Generate cubic polynomial trajectories for both the joints ------%
 
q1 = (pi*t^3)/500 - (3*pi*t^2)/100 - (6189958033024885*t)/10141204801825835211973625643008 + pi;
q2 = (pi*t^3)/1000 - (3*pi*t^2)/200 - (6189958033024885*t)/20282409603651670423947251286016 + pi/2;
q1_dot = (3*pi*t^2)/500 - (3*pi*t)/50 - 6189958033024885/10141204801825835211973625643008;
q2_dot = (3*pi*t^2)/1000 - (3*pi*t)/100 - 6189958033024885/20282409603651670423947251286016;
q1_ddot = (3*pi*t)/250 - (3*pi)/50; 
q2_ddot = (3*pi*t)/500 - (3*pi)/100;

q_desired = [q1;q2];
qdot_desired = [q1_dot;q2_dot];
qddot_desired = [q1_ddot;q2_ddot];
Q_desired = [q_desired; qdot_desired];

% ------- Adaptive Inverse Dynamics Control --------%
B = [0,0;
    0,0;
    1,0;
    0,1];

K = [6.0000,0, 5.0000, 0;
     0, 6.0000, 0, 5.0000];

P = [1.1167,        0,    0.0833,         0;
         0,    1.1167,         0,    0.0833;
    0.0833,         0,    0.1167,         0;
         0,    0.0833,         0,    0.1167];

gamma = eye(5)*0.03;

error = jointValues - Q_desired;

%--------------- Virtual control input ----------------------%
V = -K*error + qddot_desired;

%--------- Manipulator Equation Form ---------%
alpha_1 = alpha_1_hat;
alpha_2 = alpha_2_hat;
alpha_3 = alpha_3_hat;
alpha_4 = alpha_4_hat;
alpha_5 = alpha_5_hat;

Mmat_hat = [alpha_1+2*alpha_2* cos(theta2), alpha_3+alpha_2* cos(theta2); alpha_3+alpha_2* cos(theta2), alpha_3];
Cmat_hat = [-alpha_2* sin(theta2)*theta2_dot, -alpha_2* sin(theta2)*(theta1_dot+theta2_dot); alpha_2* sin(theta2)*theta1_dot,0];
Gmat_hat = [-alpha_4*g*sin(theta1)-alpha_5*g*sin(theta1+theta2); -alpha_5*g*sin(theta1+theta2)];

%-------------- The overall control law ---------------%
Tau = Mmat_hat*V + Cmat_hat*[theta1_dot; theta2_dot] + Gmat_hat;
t1 = Tau(1);
t2 = Tau(2);

%------- Integrating the system dynamics ---------%
dX(1) = theta1_dot;
dX(2) = theta2_dot;
dX(3) = (I2*t1 - I2*t2 + m2*r2^2*t1*cos(theta1 + theta2)^2 - m2*r2^2*t2*cos(theta1 + theta2)^2 + m2*r2^2*t1*sin(theta1 + theta2)^2 - m2*r2^2*t2*sin(theta1 + theta2)^2 + I2*g*l1*m2*sin(theta1) + I2*g*m1*r1*sin(theta1) + l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot^2 + l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta2_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta2_dot^2 - l1*m2*r2*t2*sin(theta1)*sin(theta1 + theta2) + g*l1*m2^2*r2^2*sin(theta1)*cos(theta1 + theta2)^2 - l1*m2*r2*t2*cos(theta1)*cos(theta1 + theta2) - l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta2_dot^2 - I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta2_dot^2 + g*m1*m2*r1*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*m1*m2*r1*r2^2*sin(theta1)*sin(theta1 + theta2)^2 - g*l1*m2^2*r2^2*cos(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2) + 2*l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot*theta2_dot - 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot*theta2_dot + l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot^2 + l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta2_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot^2 - l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta2_dot^2 - l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + 2*l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot + 2*I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot*theta2_dot)/(I1*I2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*sin(theta1 + theta2)^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2));
dX(4) = (I1*t2 - I2*t1 + I2*t2 - m2*r2^2*t1*cos(theta1 + theta2)^2 + m2*r2^2*t2*cos(theta1 + theta2)^2 - m2*r2^2*t1*sin(theta1 + theta2)^2 + m2*r2^2*t2*sin(theta1 + theta2)^2 + l1^2*m2*t2*cos(theta1)^2 + m1*r1^2*t2*cos(theta1)^2 + l1^2*m2*t2*sin(theta1)^2 + m1*r1^2*t2*sin(theta1)^2 - I2*g*l1*m2*sin(theta1) - I2*g*m1*r1*sin(theta1) + I1*g*m2*r2*sin(theta1 + theta2) - l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot^2 - l1^3*m2^2*r2*cos(theta1)^3*sin(theta1 + theta2)*theta1_dot^2 + l1^3*m2^2*r2*sin(theta1)^3*cos(theta1 + theta2)*theta1_dot^2 - l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta2_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta2_dot^2 - l1*m2*r2*t1*sin(theta1)*sin(theta1 + theta2) + 2*l1*m2*r2*t2*sin(theta1)*sin(theta1 + theta2) - g*l1*m2^2*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*l1^2*m2^2*r2*cos(theta1)^2*sin(theta1 + theta2) - l1*m2*r2*t1*cos(theta1)*cos(theta1 + theta2) + 2*l1*m2*r2*t2*cos(theta1)*cos(theta1 + theta2) + 2*l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta2_dot^2 - g*l1^2*m2^2*r2*cos(theta1)*sin(theta1)*cos(theta1 + theta2) - I1*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 + I1*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot^2 + I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta2_dot^2 + I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta2_dot^2 - g*m1*m2*r1*r2^2*sin(theta1)*cos(theta1 + theta2)^2 + g*m1*m2*r1^2*r2*cos(theta1)^2*sin(theta1 + theta2) - g*m1*m2*r1*r2^2*sin(theta1)*sin(theta1 + theta2)^2 + g*m1*m2*r1^2*r2*sin(theta1)^2*sin(theta1 + theta2) + l1^3*m2^2*r2*cos(theta1)^2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 + g*l1*m2^2*r2^2*cos(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2) - l1^3*m2^2*r2*cos(theta1)*sin(theta1)^2*sin(theta1 + theta2)*theta1_dot^2 - 2*l1*m2^2*r2^3*cos(theta1)*sin(theta1 + theta2)^3*theta1_dot*theta2_dot + 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)^3*theta1_dot*theta2_dot - l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot^2 - l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta2_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot^2 + l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta2_dot^2 + 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot^2 + l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta2_dot^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot^2 - l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta2_dot^2 - 2*l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot^2 - l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta2_dot^2 - 2*l1*m2^2*r2^3*cos(theta1)*cos(theta1 + theta2)^2*sin(theta1 + theta2)*theta1_dot*theta2_dot + 2*l1*m2^2*r2^3*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot + 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)^2*theta1_dot*theta2_dot - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*sin(theta1 + theta2)^2*theta1_dot*theta2_dot - g*l1*m1*m2*r1*r2*sin(theta1)^2*sin(theta1 + theta2) - 2*l1^2*m2^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot*theta2_dot - l1*m1*m2*r1^2*r2*cos(theta1)^3*sin(theta1 + theta2)*theta1_dot^2 + l1*m1*m2*r1^2*r2*sin(theta1)^3*cos(theta1 + theta2)*theta1_dot^2 + 2*l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)*sin(theta1 + theta2)*theta1_dot*theta2_dot - 2*I2*l1*m2*r2*cos(theta1)*sin(theta1 + theta2)*theta1_dot*theta2_dot + 2*I2*l1*m2*r2*sin(theta1)*cos(theta1 + theta2)*theta1_dot*theta2_dot - g*l1*m1*m2*r1*r2*cos(theta1)*sin(theta1)*cos(theta1 + theta2) + l1*m1*m2*r1^2*r2*cos(theta1)^2*sin(theta1)*cos(theta1 + theta2)*theta1_dot^2 - l1*m1*m2*r1^2*r2*cos(theta1)*sin(theta1)^2*sin(theta1 + theta2)*theta1_dot^2)/(I1*I2 + I2*l1^2*m2*cos(theta1)^2 + I2*m1*r1^2*cos(theta1)^2 + I2*l1^2*m2*sin(theta1)^2 + I2*m1*r1^2*sin(theta1)^2 + I1*m2*r2^2*cos(theta1 + theta2)^2 + I1*m2*r2^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + l1^2*m2^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*cos(theta1)^2*sin(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*cos(theta1 + theta2)^2 + m1*m2*r1^2*r2^2*sin(theta1)^2*sin(theta1 + theta2)^2 - 2*l1^2*m2^2*r2^2*cos(theta1)*sin(theta1)*cos(theta1 + theta2)*sin(theta1 + theta2));

%------- Integrating the adaption law -------------%
theta1_ddot = dX(3);
theta2_ddot = dX(4);

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

phi = Mmat_hat \ Y;
alpha_hat_dot = -gamma \(phi'*B'*P*error);

dX(5) = alpha_hat_dot(1);
dX(6) = alpha_hat_dot(2);
dX(7) = alpha_hat_dot(3);
dX(8) = alpha_hat_dot(4);
dX(9) = alpha_hat_dot(5);

end