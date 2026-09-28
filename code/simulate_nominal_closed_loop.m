clc
clear
close all

%% =========================================================
% SYMBOLIC / SEMI-ANALYTICAL
% STOCHASTIC NOMINAL OUTPUT-FEEDBACK TEST
%
% Publication version
%
% True parameter:
%       A1 = 7
%
% Regulator map:
%       symbolic_regulator_map_A1
%
% Cortical input:
%       p(t) ~ N(90,30^2)
%
% Controller ON at t = 5 s
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
fprintf(' SYMBOLIC STOCHASTIC NOMINAL TEST\n')
fprintf('=========================================\n')

fprintf('Original Lw = %.6f\n',Lw_A1);
fprintf('Used Lw     = %.6f\n',Lw_test);

%% =========================================================
% Simulation settings
% ==========================================================
T  = 20;
dt = 0.001;

tspan = 0:dt:T;

A1_true = 7;

t_ctrl   = 5;
t_region = 5.7;

%% =========================================================
% Cortical input
% Mean = 90
% Standard deviation = 30
% ==========================================================
rng(1);

P1 = normrnd(90,30,1,length(tspan));
P2 = normrnd(90,30,1,length(tspan));

% fprintf('\nCortical input:\n')
% fprintf('mean(p) = %.6f\n',mean(p1));
% fprintf('std(p)  = %.6f\n',std(p2));

%% =========================================================
% Initial conditions
% ==========================================================
z0 = zeros(33,1);

% True system states
z0(1:16) = 0;

% Observer states
z0(17:32) = 1e-3;

% Initial A1 estimate
z0(33) = 6.9;

%% =========================================================
% Solver
% ==========================================================
opts = odeset( ...
    'RelTol',1e-3, ...
    'AbsTol',1e-5);

%% =========================================================
% Simulation
% ==========================================================
[t,z] = ode45( ...
    @(t,z) symbolic_closed_loop_nominal( ...
        t,z,K_nom,Lx_A1,Lw_test, ...
        P1,P2,A1_true,t_ctrl,dt), ...
    tspan,z0,opts);

%% =========================================================
% Extract states
% ==========================================================
x      = z(:,1:16);
xhat   = z(:,17:32);
A1_hat = z(:,33);

%% Cortical outputs
y1 = x(:,2)-x(:,3);
y2 = x(:,8)-x(:,9);

y1_hat = xhat(:,2)-xhat(:,3);

%% Propagation state
x13 = x(:,13);

%% =========================================================
% Reconstruct regulator reference and control
% ==========================================================
N = length(t);

u_hist    = zeros(N,1);
pi13_hist = zeros(N,1);
c_hist    = zeros(N,1);

for i = 1:N

    [pi_hat,c_hat] = ...
        symbolic_regulator_map_A1(A1_hat(i));

    pi13_hist(i) = pi_hat(13);
    c_hist(i)    = c_hat;

    if t(i) < t_ctrl

        u_hist(i) = 0;

    else

        u_hist(i) = ...
            c_hat + ...
            K_nom*(xhat(i,:)'-pi_hat);

    end

end

%% =========================================================
% Errors
% ==========================================================
reg_err = x13-pi13_hist;

output_est_err = y1_hat-y1;

state_est_err = ...
    vecnorm(xhat-x,2,2);

A1_est_err = ...
    A1_hat-A1_true;

%% =========================================================
% Regions
% ==========================================================
idx_before = t < t_ctrl;
idx_after  = t >= t_ctrl;
idx_region = t >= t_region;

%% =========================================================
% Regulation metrics
% ==========================================================
RMS_x13_before = ...
    sqrt(mean(reg_err(idx_before).^2));

RMS_x13_after = ...
    sqrt(mean(reg_err(idx_after).^2));

RMS_x13_region = ...
    sqrt(mean(reg_err(idx_region).^2));

Final_x13_error = ...
    abs(reg_err(end));

%% =========================================================
% Observer metrics
% ==========================================================
RMS_output_est = ...
    sqrt(mean(output_est_err(idx_region).^2));

Final_output_est = ...
    abs(output_est_err(end));

RMS_state_est = ...
    sqrt(mean(state_est_err(idx_region).^2));

Final_state_est = ...
    state_est_err(end);

%% =========================================================
% Parameter-estimation metrics
% ==========================================================
Final_A1hat = ...
    A1_hat(end);

Final_A1error = ...
    abs(A1_hat(end)-A1_true);

RMS_A1error = ...
    sqrt(mean(A1_est_err(idx_region).^2));

%% =========================================================
% Control metrics
% ==========================================================
MaxControl = ...
    max(abs(u_hist));

RMS_Control = ...
    sqrt(mean(u_hist(idx_region).^2));

FinalControl = ...
    u_hist(end);

%% =========================================================
% Print results
% ==========================================================
fprintf('\n=========================================\n')
fprintf(' SYMBOLIC STOCHASTIC PERFORMANCE\n')
fprintf('=========================================\n')

fprintf('\n--- Regulation ---\n')

fprintf('Target x13                 = %.12e\n', ...
    pi13_hist(end));

fprintf('RMS x13 error before ctrl  = %.12e\n', ...
    RMS_x13_before);

fprintf('RMS x13 error after ctrl   = %.12e\n', ...
    RMS_x13_after);

fprintf('RMS x13 error t>=5.7      = %.12e\n', ...
    RMS_x13_region);

fprintf('Final x13                  = %.12e\n', ...
    x13(end));

fprintf('Final |x13-pi13|           = %.12e\n', ...
    Final_x13_error);

fprintf('\n--- Parameter estimation ---\n')

fprintf('Initial A1_hat             = %.12f\n', ...
    A1_hat(1));

fprintf('Final A1_hat               = %.12f\n', ...
    Final_A1hat);

fprintf('Final |A1hat-A1|           = %.12e\n', ...
    Final_A1error);

fprintf('RMS A1 estimation error    = %.12e\n', ...
    RMS_A1error);

fprintf('\n--- Observer ---\n')

fprintf('RMS output estimation err  = %.12e\n', ...
    RMS_output_est);

fprintf('Final output estimation err= %.12e\n', ...
    Final_output_est);

fprintf('RMS state estimation err   = %.12e\n', ...
    RMS_state_est);

fprintf('Final state estimation err = %.12e\n', ...
    Final_state_est);

fprintf('\n--- Control ---\n')

fprintf('Maximum |u|                = %.12e\n', ...
    MaxControl);

fprintf('RMS u t>=5.7              = %.12e\n', ...
    RMS_Control);

fprintf('Final u                    = %.12e\n', ...
    FinalControl);

%% =========================================================
% PUBLICATION FIGURE
% Three-panel closed-loop response
% ==========================================================

fig = figure( ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[2 2 17.5 15]);

tl = tiledlayout(3,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');

%% =========================================================
% (a) Epileptic-column output
% ==========================================================
ax1 = nexttile;

plot(t,y1,'LineWidth',1.0);

hold on

xline( ...
    t_ctrl, ...
    '--', ...
    'LineWidth',1);

ylabel( ...
    '$y_1$ (mV)', ...
    'Interpreter','latex');

text( ...
    0.015,0.90,'(a)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

grid on
box on

set(ax1, ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'LineWidth',0.8);

%% =========================================================
% (b) Second-column output
% ==========================================================
ax2 = nexttile;

plot(t,y2,'LineWidth',1.0);

hold on

xline( ...
    t_ctrl, ...
    '--', ...
    'LineWidth',1);

ylabel( ...
    '$y_2$ (mV)', ...
    'Interpreter','latex');

text( ...
    0.015,0.90,'(b)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

grid on
box on

set(ax2, ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'LineWidth',0.8);

%% =========================================================
% (c) Propagation state
% ==========================================================
ax3 = nexttile;

plot( ...
    t,x13, ...
    'LineWidth',1.1);

hold on

yline( ...
    1e-3, ...
    '--', ...
    'LineWidth',1);

xline( ...
    t_ctrl, ...
    '--', ...
    'LineWidth',1);

xlabel('Time (s)');

ylabel( ...
    '$x_{13}$', ...
    'Interpreter','latex');

legend( ...
    '$x_{13}(t)$', ...
    '$\varepsilon=10^{-3}$', ...
    'Interpreter','latex', ...
    'Location','best', ...
    'Box','off');

text( ...
    0.015,0.90,'(c)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',10);

grid on
box on

set(ax3, ...
    'FontName','Times New Roman', ...
    'FontSize',9, ...
    'LineWidth',0.8);

%% =========================================================
% Link time axes
% ==========================================================
linkaxes([ax1 ax2 ax3],'x');

xlim([0 T]);

%% =========================================================
% Controller activation label
% Only once to keep figure clean
% ==========================================================
yl = ylim(ax1);

text( ...
    ax1, ...
    t_ctrl+0.15, ...
    yl(2)-0.10*(yl(2)-yl(1)), ...
    'Controller ON', ...
    'FontName','Times New Roman', ...
    'FontSize',8);

%% =========================================================
% Export publication figure
% ==========================================================
exportgraphics( ...
    fig, ...
    'symbolic_nominal_closed_loop.pdf', ...
    'ContentType','vector');

exportgraphics( ...
    fig, ...
    'symbolic_nominal_closed_loop.png', ...
    'Resolution',600);

fprintf('\nPublication figure saved as:\n');
fprintf('symbolic_nominal_closed_loop.pdf\n');
fprintf('symbolic_nominal_closed_loop.png\n');

%% =========================================================
% CLOSED-LOOP SYSTEM
% ==========================================================
function dz = symbolic_closed_loop_nominal( ...
    t,z,K,Lx,Lw,P1,P2,A1_true,t_ctrl,dt)

%% States
x      = z(1:16);
xhat   = z(17:32);
A1_hat = z(33);

%% =========================================================
% Fixed parameters
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
idxP = ...
    ceil(t/dt + 0.001);

idxP = ...
    max(1,min(idxP,length(P1)));

p1 = P1(idxP);
p2 = P2(idxP);

%% =========================================================
% Measurement
% ==========================================================
y = x(2)-x(3);

yhat = xhat(2)-xhat(3);

ytilde = yhat-y;

%% =========================================================
% SEMI-ANALYTICAL REGULATOR MAP
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
dx = symbolic_jansen_rhs( ...
    x,u,p1,p2, ...
    A1_true,B1,A2,B2, ...
    a,b,C1,C2,C3,C4, ...
    Ad,ad,k12,k21);

%% =========================================================
% Observer
% ==========================================================
dxhat = symbolic_jansen_rhs( ...
    xhat,u,p1,p2, ...
    A1_hat,B1,A2,B2, ...
    a,b,C1,C2,C3,C4, ...
    Ad,ad,k12,k21);

dxhat = ...
    dxhat + Lx*ytilde;

%% Parameter observer
dA1_hat = ...
    Lw*ytilde;

%% Complete augmented dynamics
dz = [ ...
    dx
    dxhat
    dA1_hat ];

end

%% =========================================================
% TWO-COLUMN NONLINEAR JANSEN-RIT MODEL
% ==========================================================
function dx = symbolic_jansen_rhs( ...
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
% Pathway 1 -> 2
% ==========================================================
dx(13) = x(15);

dx(15) = ...
    Ad*ad*SJR(x(2)-x(3)) ...
    -2*ad*x(15) ...
    -ad^2*x(13);

%% =========================================================
% Pathway 2 -> 1
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