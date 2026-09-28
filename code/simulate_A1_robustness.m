clc
clear
close all

%% =========================================================

% ==========================================================

%% =========================================================
% Load controller and observer gains
% ==========================================================
load('K_nominal.mat','K_nom');

load('observer_A1_gains.mat', ...
    'Lx_A1','Lw_A1');

%% Tuned parameter-estimation gain
Lw_test = 0.05*Lw_A1;

fprintf('\n=========================================\n')
fprintf(' SYMBOLIC A1 ROBUSTNESS TEST\n')
fprintf('=========================================\n')

fprintf('Original Lw = %.6f\n',Lw_A1);
fprintf('Used Lw     = %.6f\n',Lw_test);

%% =========================================================
% Simulation settings
% ==========================================================
T  = 20;
dt = 0.001;

tspan = 0:dt:T;

t_ctrl   = 5;
t_region = 5.7;

%% =========================================================
% True A1 values
% ==========================================================
A1_values = [6.5 7.0 7.5];

nA = length(A1_values);

%% =========================================================
% Same cortical-input realization for all A1 values
% ==========================================================
rng(1);

P1 = normrnd(90,30,1,length(tspan));
P2 = normrnd(90,30,1,length(tspan));

fprintf('\nCortical inputs:\n')
fprintf('mean(P1) = %.6f, std(P1) = %.6f\n',mean(P1),std(P1));
fprintf('mean(P2) = %.6f, std(P2) = %.6f\n',mean(P2),std(P2));
%% =========================================================
% Performance storage
% ==========================================================
RMS_x13_error       = zeros(nA,1);
Final_x13_error     = zeros(nA,1);

A1hat_final         = zeros(nA,1);
Final_A1_error      = zeros(nA,1);
RMS_A1_error        = zeros(nA,1);

RMS_output_est_err  = zeros(nA,1);

RMS_state_est_err   = zeros(nA,1);
Final_state_est_err = zeros(nA,1);

MaxControl          = zeros(nA,1);
RMS_Control         = zeros(nA,1);

%% =========================================================
% Trajectory storage
% ==========================================================
t_save              = cell(nA,1);
x13_save            = cell(nA,1);
A1hat_save          = cell(nA,1);
u_save              = cell(nA,1);
output_est_err_save = cell(nA,1);

%% =========================================================
% Solver
% ==========================================================
opts = odeset( ...
    'RelTol',1e-3, ...
    'AbsTol',1e-5);

%% =========================================================
% Loop over true A1 values
% ==========================================================
for k = 1:nA

    A1_true = A1_values(k);

    fprintf('\n=========================================\n')
    fprintf(' Running A1 = %.2f\n',A1_true);
    fprintf('=========================================\n')

    %% =====================================================
    % Initial conditions
    % ======================================================
    z0 = zeros(33,1);

    % True system
    z0(1:16) = 0;

    % Observer
    z0(17:32) = 1e-3;

    % Same initial A1 estimate for every case
    z0(33) = 6.9;

    %% =====================================================
    % Simulation
    % ======================================================
 [t,z] = ode45( ...
    @(t,z) symbolic_robust_closed_loop( ...
        t,z,K_nom,Lx_A1,Lw_test, ...
        P1,P2,A1_true,t_ctrl,dt), ...
    tspan,z0,opts);
    %% =====================================================
    % Extract states
    % ======================================================
    x      = z(:,1:16);
    xhat   = z(:,17:32);
    A1_hat = z(:,33);

    %% Measured cortical output
    y1 = ...
        x(:,2)-x(:,3);

    %% Estimated cortical output
    y1_hat = ...
        xhat(:,2)-xhat(:,3);

    %% Propagation state
    x13 = ...
        x(:,13);

    %% =====================================================
    % Reconstruct regulator and control
    % ======================================================
    N = length(t);

    u_hist    = zeros(N,1);
    pi13_hist = zeros(N,1);

    for i = 1:N

        [pi_hat,c_hat] = ...
            symbolic_regulator_map_A1(A1_hat(i));

        pi13_hist(i) = ...
            pi_hat(13);

        if t(i) < t_ctrl

            u_hist(i) = 0;

        else

            u_hist(i) = ...
                c_hat + ...
                K_nom*(xhat(i,:)'-pi_hat);

        end

    end

    %% =====================================================
    % Errors
    % ======================================================
    reg_err = ...
        x13-pi13_hist;

    output_est_err = ...
        y1_hat-y1;

    state_est_err = ...
        vecnorm(xhat-x,2,2);

    A1_est_err = ...
        A1_hat-A1_true;

    %% =====================================================
    % Evaluation region
    % ======================================================
    idx_region = ...
        t >= t_region;

    %% =====================================================
    % Performance metrics
    % ======================================================

    %% Regulation
    RMS_x13_error(k) = ...
        sqrt(mean(reg_err(idx_region).^2));

    Final_x13_error(k) = ...
        abs(reg_err(end));

    %% Parameter estimation
    A1hat_final(k) = ...
        A1_hat(end);

    Final_A1_error(k) = ...
        abs(A1_hat(end)-A1_true);

    RMS_A1_error(k) = ...
        sqrt(mean(A1_est_err(idx_region).^2));

    %% Output observer
    RMS_output_est_err(k) = ...
        sqrt(mean(output_est_err(idx_region).^2));

    %% State observer
    RMS_state_est_err(k) = ...
        sqrt(mean(state_est_err(idx_region).^2));

    Final_state_est_err(k) = ...
        state_est_err(end);

    %% Control effort
    MaxControl(k) = ...
        max(abs(u_hist));

    RMS_Control(k) = ...
        sqrt(mean(u_hist(idx_region).^2));

    %% =====================================================
    % Save trajectories
    % ======================================================
    t_save{k} = ...
        t(:);

    x13_save{k} = ...
        x13(:);

    A1hat_save{k} = ...
        A1_hat(:);

    u_save{k} = ...
        u_hist(:);

    output_est_err_save{k} = ...
        output_est_err(:);

    %% =====================================================
    % Check trajectory lengths
    % ======================================================
    fprintf('Stored samples         = %d\n', ...
        length(t_save{k}));

    fprintf('Stored observer errors = %d\n', ...
        length(output_est_err_save{k}));

    %% =====================================================
    % Report
    % ======================================================
    fprintf('RMS x13 error t>=5.7 = %.12e\n', ...
        RMS_x13_error(k));

    fprintf('Final x13 error       = %.12e\n', ...
        Final_x13_error(k));

    fprintf('Final A1_hat          = %.12f\n', ...
        A1hat_final(k));

    fprintf('Final |A1hat-A1|      = %.12e\n', ...
        Final_A1_error(k));

    fprintf('RMS A1 error          = %.12e\n', ...
        RMS_A1_error(k));

    fprintf('RMS output est error  = %.12e\n', ...
        RMS_output_est_err(k));

    fprintf('RMS state est error   = %.12e\n', ...
        RMS_state_est_err(k));

    fprintf('Final state est error = %.12e\n', ...
        Final_state_est_err(k));

    fprintf('Maximum |u|           = %.12e\n', ...
        MaxControl(k));

end

%% =========================================================
% Summary table
% ==========================================================
Results = table( ...
    A1_values(:), ...
    RMS_x13_error, ...
    Final_x13_error, ...
    A1hat_final, ...
    Final_A1_error, ...
    RMS_A1_error, ...
    RMS_output_est_err, ...
    RMS_state_est_err, ...
    Final_state_est_err, ...
    MaxControl, ...
    RMS_Control, ...
    'VariableNames',{ ...
    'A1_true', ...
    'RMS_x13_error', ...
    'Final_x13_error', ...
    'A1hat_final', ...
    'Final_A1_error', ...
    'RMS_A1_error', ...
    'RMS_output_est_error', ...
    'RMS_state_est_error', ...
    'Final_state_est_error', ...
    'MaxControl', ...
    'RMS_Control'});

fprintf('\n\n=========================================\n')
fprintf(' SYMBOLIC A1 ROBUSTNESS SUMMARY\n')
fprintf('=========================================\n')

disp(Results)

%% =========================================================
% Verify stored trajectory dimensions
% ==========================================================
fprintf('\n=========================================\n')
fprintf(' SAVED TRAJECTORY DIMENSIONS\n')
fprintf('=========================================\n')

for k = 1:nA

    fprintf( ...
        'A1 = %.1f : t = %d , x13 = %d , err = %d\n', ...
        A1_values(k), ...
        length(t_save{k}), ...
        length(x13_save{k}), ...
        length(output_est_err_save{k}));

end

%% =========================================================
% Save results AFTER trajectories have been generated
% ==========================================================
save( ...
    'symbolic_A1_robustness_results.mat', ...
    'A1_values', ...
    't_save', ...
    'x13_save', ...
    'output_est_err_save', ...
    'A1hat_save', ...
    'u_save', ...
    'RMS_x13_error', ...
    'Final_x13_error', ...
    'A1hat_final', ...
    'Final_A1_error', ...
    'RMS_A1_error', ...
    'RMS_output_est_err', ...
    'RMS_state_est_err', ...
    'Final_state_est_err', ...
    'MaxControl', ...
    'RMS_Control');

fprintf('\nRobustness results saved successfully.\n');
fprintf('File: symbolic_A1_robustness_results.mat\n');

%% =========================================================
% CLOSED-LOOP SYSTEM
% ==========================================================
function dz = symbolic_robust_closed_loop( ...
    t,z,K,Lx,Lw,P1,P2,A1_true,t_ctrl,dt)

%% States
x      = z(1:16);
xhat   = z(17:32);
A1_hat = z(33);

%% =========================================================
% Parameters
% ==========================================================
B1 = 22;

A2 = 3.25;
B2 = 22;

a = 100;
b = 50;

C1 = 135;
C2 = 108;
C3 = 33.75;
C4 = 33.75;

Ad = 3;
ad = 30;

k12 = 4000;
k21 = 1;

%% =========================================================
% Cortical input
% ==========================================================
idxP = ceil(t/dt + 0.001);

idxP = max(1,min(idxP,length(P1)));

p1 = P1(idxP);
p2 = P2(idxP);

%% =========================================================
% Measured output
% ==========================================================
y = ...
    x(2)-x(3);

yhat = ...
    xhat(2)-xhat(3);

ytilde = ...
    yhat-y;

%% =========================================================
% Semi-analytical regulator map
% ==========================================================
[pi_hat,c_hat] = ...
    symbolic_regulator_map_A1(A1_hat);

%% =========================================================
% Controller
% ==========================================================
if t < t_ctrl

    u = 0;

else

    u = ...
        c_hat + ...
        K*(xhat-pi_hat);

end

%% =========================================================
% True nonlinear system
% ==========================================================
dx = symbolic_robust_jansen_rhs( ...
    x,u,p1,p2, ...
    A1_true,B1,A2,B2, ...
    a,b,C1,C2,C3,C4, ...
    Ad,ad,k12,k21);

%% =========================================================
% Observer
% ==========================================================
dxhat = symbolic_robust_jansen_rhs( ...
    xhat,u,p1,p2, ...
    A1_hat,B1,A2,B2, ...
    a,b,C1,C2,C3,C4, ...
    Ad,ad,k12,k21);

dxhat = ...
    dxhat + Lx*ytilde;

%% =========================================================
% Parameter observer
% ==========================================================
dA1_hat = ...
    Lw*ytilde;

%% Complete augmented system
dz = [ ...
    dx
    dxhat
    dA1_hat ];

end

%% =========================================================
% NONLINEAR TWO-COLUMN JANSEN-RIT MODEL
% ==========================================================
function dx = symbolic_robust_jansen_rhs( ...
    x,u,p1,p2, ...
    A1,B1,A2,B2, ...
    a,b,C1,C2,C3,C4, ...
    Ad,ad,k12,k21)

dx = zeros(16,1);

%% =========================================================
% Column 1
% ==========================================================
dx(1) = x(4);
dx(2) = x(5);
dx(3) = x(6);

dx(4) = ...
    A1*a*SJR(x(2)-x(3)) ...
    -2*a*x(4) ...
    -a^2*x(1);

dx(5) = ...
    A1*a*( ...
        p1 ...
        +u ...
        +C2*SJR(C1*x(1)) ...
        +k21*x(14)) ...
    -2*a*x(5) ...
    -a^2*x(2);

dx(6) = ...
    B1*b*C4*SJR(C3*x(1)) ...
    -2*b*x(6) ...
    -b^2*x(3);

%% =========================================================
% Column 2
% ==========================================================
dx(7) = x(10);
dx(8) = x(11);
dx(9) = x(12);

dx(10) = ...
    A2*a*SJR(x(8)-x(9)) ...
    -2*a*x(10) ...
    -a^2*x(7);

dx(11) = ...
    A2*a*( ...
        p2 ...
        +C2*SJR(C1*x(7)) ...
        +k12*x(13)) ...
    -2*a*x(11) ...
    -a^2*x(8);

dx(12) = ...
    B2*b*C4*SJR(C3*x(7)) ...
    -2*b*x(12) ...
    -b^2*x(9);

%% =========================================================
% Inter-column pathway 1 -> 2
% ==========================================================
dx(13) = x(15);

dx(15) = ...
    Ad*ad*SJR(x(2)-x(3)) ...
    -2*ad*x(15) ...
    -ad^2*x(13);

%% =========================================================
% Inter-column pathway 2 -> 1
% ==========================================================
dx(14) = x(16);

dx(16) = ...
    Ad*ad*SJR(x(8)-x(9)) ...
    -2*ad*x(16) ...
    -ad^2*x(14);

end

%% =========================================================
% ORIGINAL JANSEN-RIT SIGMOID
% ==========================================================
function r = SJR(v)

r = ...
    5 ./ ...
    (1 + exp(0.56*(6-v)));

end