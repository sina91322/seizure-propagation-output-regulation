clc
clear
close all

%% =========================================================
% Nominal parameters
% ==========================================================
a = 100;
b = 50;

A1 = 7;
B1 = 22;

A2 = 3.25;
B2 = 22;

C1 = 135;
C2 = 108;
C3 = 33.75;
C4 = 33.75;

Ad = 3;
ad = 30;

k12 = 4000;
k21 = 1;

%% =========================================================
% Nominal regulator solution
% ==========================================================
pi_nom = [ ...
     0.000700000000
    -2.569358224
     2.524581232
     0
     0
     0
     0.01093013439
     4.345459079
     3.041036007
     0
     0
     0
     0.001000000000
     0.03363118274
     0
     0 ];

c_nom = -145.8164630;

%% =========================================================
% Sigmoid and derivative
% ==========================================================
S = @(v) 5 ./ (1 + exp(0.56*(6-v)));

dS = @(v) ...
    0.56 .* S(v) .* (1 - S(v)/5);

%% =========================================================
% Sigmoid slopes at nominal operating point
% ==========================================================
q1 = dS(pi_nom(2)-pi_nom(3));
q2 = dS(C1*pi_nom(1));
q3 = dS(C3*pi_nom(1));

q4 = dS(pi_nom(8)-pi_nom(9));
q5 = dS(C1*pi_nom(7));
q6 = dS(C3*pi_nom(7));

%% =========================================================
% Nominal Jacobian
% ==========================================================
A_nom = zeros(16);

% Kinematic equations
A_nom(1,4) = 1;
A_nom(2,5) = 1;
A_nom(3,6) = 1;

A_nom(7,10) = 1;
A_nom(8,11) = 1;
A_nom(9,12) = 1;

A_nom(13,15) = 1;
A_nom(14,16) = 1;

% Column 1
A_nom(4,1) = -a^2;
A_nom(4,2) =  A1*a*q1;
A_nom(4,3) = -A1*a*q1;
A_nom(4,4) = -2*a;

A_nom(5,1)  = A1*a*C2*C1*q2;
A_nom(5,2)  = -a^2;
A_nom(5,5)  = -2*a;
A_nom(5,14) = A1*a*k21;

A_nom(6,1) = B1*b*C4*C3*q3;
A_nom(6,3) = -b^2;
A_nom(6,6) = -2*b;

% Column 2
A_nom(10,7)  = -a^2;
A_nom(10,8)  =  A2*a*q4;
A_nom(10,9)  = -A2*a*q4;
A_nom(10,10) = -2*a;

A_nom(11,7)  = A2*a*C2*C1*q5;
A_nom(11,8)  = -a^2;
A_nom(11,11) = -2*a;
A_nom(11,13) = A2*a*k12;

A_nom(12,7)  = B2*b*C4*C3*q6;
A_nom(12,9)  = -b^2;
A_nom(12,12) = -2*b;

% Inter-column pathways
A_nom(15,2)  =  Ad*ad*q1;
A_nom(15,3)  = -Ad*ad*q1;
A_nom(15,13) = -ad^2;
A_nom(15,15) = -2*ad;

A_nom(16,8)  =  Ad*ad*q4;
A_nom(16,9)  = -Ad*ad*q4;
A_nom(16,14) = -ad^2;
A_nom(16,16) = -2*ad;

%% =========================================================
% Input matrix
% ==========================================================
B = zeros(16,1);
B(5) = A1*a;

%% =========================================================
% Nominal LQR design
% ==========================================================
Q = eye(16);
R = 100;

[K_lqr,~,~] = lqr(A_nom,B,Q,R);

% MATLAB convention:
% u = -K_lqr*x
%
% Our convention:
% u = c_nom + K_nom*(x - pi_nom)

K_nom = -K_lqr;

%% =========================================================
% Verification
% ==========================================================
Acl = A_nom + B*K_nom;

eig_cl = eig(Acl);
max_real_cl = max(real(eig_cl));

disp('K_nom =')
disp(K_nom)

disp('Closed-loop eigenvalues =')
disp(eig_cl)

disp('Maximum real part =')
disp(max_real_cl)

if max_real_cl < 0
    disp('Nominal closed-loop linear model is stable.')
else
    warning('Nominal closed-loop linear model is NOT stable.')
end

%% =========================================================
% Save controller gain
% ==========================================================
save('K_nominal.mat','K_nom','A_nom','B','pi_nom','c_nom');