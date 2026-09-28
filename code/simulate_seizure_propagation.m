clc;
clear;
close all;

rng(1);   % Reproducibility

%% =========================================================
% Simulation settings
% ==========================================================
T  = 10;          % [s]
dt = 0.001;       % [s]
time = 0:dt:T;

x0 = zeros(16,1);

options = odeset( ...
    'RelTol',1e-1, ...
    'AbsTol',1e-5);

%% =========================================================
% Independent stochastic inputs for the two columns
% ==========================================================
mu_p    = 90;
sigma_p = 30;

p1 = mu_p + sigma_p*randn(size(time));
p2 = mu_p + sigma_p*randn(size(time));

% Continuous functions used by ode45
p1_fun = @(t) interp1(time,p1,t,'previous','extrap');
p2_fun = @(t) interp1(time,p2,t,'previous','extrap');


%% =========================================================
% CASE 1: Uncoupled columns
% ==========================================================
k12 = 0;
k21 = 0;

% Warm-up simulation
[~,x_temp] = ode45( ...
    @(t,x) jansen_double(t,x,k12,k21,p1_fun,p2_fun), ...
    time,x0,options);

x0_unc = x_temp(end,:)';

% Simulation used for plotting
[t_unc,x_unc] = ode45( ...
    @(t,x) jansen_double(t,x,k12,k21,p1_fun,p2_fun), ...
    time,x0_unc,options);

% Cortical outputs
y1_unc = x_unc(:,2) - x_unc(:,3);
y2_unc = x_unc(:,8) - x_unc(:,9);


%% =========================================================
% CASE 2: Coupled columns
% ==========================================================
k12 = 4000;
k21 = 1;

% Warm-up simulation
[~,x_temp] = ode45( ...
    @(t,x) jansen_double(t,x,k12,k21,p1_fun,p2_fun), ...
    time,x0,options);

x0_coup = x_temp(end,:)';

% Simulation used for plotting
[t_coup,x_coup] = ode45( ...
    @(t,x) jansen_double(t,x,k12,k21,p1_fun,p2_fun), ...
    time,x0_coup,options);

% Cortical outputs
y1_coup = x_coup(:,2) - x_coup(:,3);
y2_coup = x_coup(:,8) - x_coup(:,9);


%% =========================================================
% Figure
% ==========================================================
figure('Color','w','Position',[100 100 1000 600]);

tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');

% (a) Column 1 -- uncoupled
nexttile;
plot(t_unc,y1_unc,'LineWidth',1);
xlim([0 T]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('(a) Column 1 -- uncoupled');
set(gca,'FontSize',11);
box on;

% (b) Column 1 -- coupled
nexttile;
plot(t_coup,y1_coup,'LineWidth',1);
xlim([0 T]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('(b) Column 1 -- coupled');
set(gca,'FontSize',11);
box on;

% (c) Column 2 -- uncoupled
nexttile;
plot(t_unc,y2_unc,'LineWidth',1);
xlim([0 T]);
xlabel('Time (s)');
ylabel('Amplitude (mV)');
title('(c) Column 2 -- uncoupled');
set(gca,'FontSize',11);
box on;

% (d) Column 2 -- coupled
nexttile;

plot(t_coup,y2_coup,'LineWidth',1);

xlim([0 T]);

xlabel('Time (s)');
ylabel('Amplitude (mV)');

title('(d) Column 2 -- coupled');

set(gca,'FontSize',11);
box on;


%% =========================================================
% Export as vector PDF
% ==========================================================

filename = 'G:\mypapers\sim.paper.dbs\Fig3.pdf';

exportgraphics(gcf, filename, ...
    'ContentType','vector');

fprintf('Figure saved to:\n%s\n', filename);


%% =========================================================
% Coupled two-column Jansen--Rit model
% ==========================================================
function dx = jansen_double(t,x,k12,k21,p1_fun,p2_fun)

%% Neural mass parameters
a = 100;
b = 50;

% Column 1: epileptic
A1 = 7;
B1 = 22;

% Column 2: normal
A2 = 3.25;
B2 = 22;

%% Intracolumn connectivity
C1 = 135;
C2 = 108;
C3 = 33.75;
C4 = 33.75;

%% Inter-column synaptic dynamics
Ad = 3;
ad = 30;

%% Independent external inputs
p1 = p1_fun(t);
p2 = p2_fun(t);

%% State equations
dx = zeros(16,1);

% =========================================================
% COLUMN 1
% ==========================================================
dx(1) = x(4);

dx(4) = ...
    A1*a*S(x(2)-x(3)) ...
    - 2*a*x(4) ...
    - a^2*x(1);


dx(2) = x(5);

dx(5) = ...
    A1*a*( ...
        p1 ...
        + C2*S(C1*x(1)) ...
        + k21*x(14)) ...
    - 2*a*x(5) ...
    - a^2*x(2);


dx(3) = x(6);

dx(6) = ...
    B1*b*C4*S(C3*x(1)) ...
    - 2*b*x(6) ...
    - b^2*x(3);


% =========================================================
% COLUMN 2
% ==========================================================
dx(7) = x(10);

dx(10) = ...
    A2*a*S(x(8)-x(9)) ...
    - 2*a*x(10) ...
    - a^2*x(7);


dx(8) = x(11);

dx(11) = ...
    A2*a*( ...
        p2 ...
        + C2*S(C1*x(7)) ...
        + k12*x(13)) ...
    - 2*a*x(11) ...
    - a^2*x(8);


dx(9) = x(12);

dx(12) = ...
    B2*b*C4*S(C3*x(7)) ...
    - 2*b*x(12) ...
    - b^2*x(9);


% =========================================================
% INTER-COLUMN PATHWAY
% Column 1 --> Column 2
% ==========================================================
dx(13) = x(15);

dx(15) = ...
    Ad*ad*S(x(2)-x(3)) ...
    - 2*ad*x(15) ...
    - ad^2*x(13);


% =========================================================
% INTER-COLUMN PATHWAY
% Column 2 --> Column 1
% ==========================================================
dx(14) = x(16);

dx(16) = ...
    Ad*ad*S(x(8)-x(9)) ...
    - 2*ad*x(16) ...
    - ad^2*x(14);

end


%% =========================================================
% Sigmoid
% ==========================================================
function r = S(v)

r = 5 ./ (1 + exp(0.56*(6-v)));

end