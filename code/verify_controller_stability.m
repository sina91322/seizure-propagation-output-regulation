clc
clear
close all

%% =========================================================
% Robust nonlinear q-envelope
% Obtained from A1 sweep with Lw = 0.05
% ==========================================================

%% =========================================================
% Robust nonlinear q-envelope
% Post-transient regulated region: t >= 5.7 s
% Independent stochastic cortical inputs
% ==========================================================

qmin = [ ...
    0.002633014238 ...
    0.093460662170 ...
    0.091487606035 ...
    0.086082578200 ...
    0.136894504953 ...
    0.100793057632 ];

qmax = [ ...
    0.010144013427 ...
    0.097945203362 ...
    0.092569140126 ...
    0.229600033962 ...
    0.228721321624 ...
    0.115581138122 ];
%% A1 uncertainty
A1_min = 6.5;
A1_max = 7.5;

%% =========================================================
% Load fixed controller gain
% ==========================================================
load('K_nominal.mat','K_nom');
fprintf('\nFinal controller gain K:\n');
disp(K_nom);

fprintf('Size of K = %d x %d\n', ...
    size(K_nom,1),size(K_nom,2));
%% =========================================================
% Fixed model parameters
% ==========================================================
a  = 100;
b  = 50;

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

nx = 16;

%% =========================================================
% Build 128 vertices:
% 6 q uncertainties + A1 uncertainty
% ==========================================================
Nvert = 128;

Acl = cell(Nvert,1);

unstableVertices = [];

worstEig = -Inf;
worstIndex = 0;
worstA1 = NaN;
worstQ = [];

for i = 1:Nvert

    bits = dec2bin(i-1,7)-'0';

    %% q corner
    qbits = bits(1:6);

    q = qmin + qbits.*(qmax-qmin);

    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q5 = q(5);
    q6 = q(6);

    %% A1 corner
    if bits(7) == 0
        A1 = A1_min;
    else
        A1 = A1_max;
    end

    %% =====================================================
    % Linearized / slope-parameterized A matrix
    % ======================================================
    A = zeros(nx);

    %% Kinematics
    A(1,4) = 1;
    A(2,5) = 1;
    A(3,6) = 1;

    A(7,10) = 1;
    A(8,11) = 1;
    A(9,12) = 1;

    A(13,15) = 1;
    A(14,16) = 1;

    %% Column 1
    A(4,1) = -a^2;
    A(4,2) =  A1*a*q1;
    A(4,3) = -A1*a*q1;
    A(4,4) = -2*a;

    A(5,1)  = A1*a*C2*C1*q2;
    A(5,2)  = -a^2;
    A(5,5)  = -2*a;
    A(5,14) = A1*a*k21;

    A(6,1) = B1*b*C4*C3*q3;
    A(6,3) = -b^2;
    A(6,6) = -2*b;

    %% Column 2
    A(10,7)  = -a^2;
    A(10,8)  =  A2*a*q4;
    A(10,9)  = -A2*a*q4;
    A(10,10) = -2*a;

    A(11,7)  = A2*a*C2*C1*q5;
    A(11,8)  = -a^2;
    A(11,11) = -2*a;
    A(11,13) = A2*a*k12;

    A(12,7)  = B2*b*C4*C3*q6;
    A(12,9)  = -b^2;
    A(12,12) = -2*b;

    %% Inter-column pathways
    A(15,2)  =  Ad*ad*q1;
    A(15,3)  = -Ad*ad*q1;
    A(15,13) = -ad^2;
    A(15,15) = -2*ad;

    A(16,8)  =  Ad*ad*q4;
    A(16,9)  = -Ad*ad*q4;
    A(16,14) = -ad^2;
    A(16,16) = -2*ad;

    %% =====================================================
    % Input matrix depends on A1
    % ======================================================
    B = zeros(nx,1);

    B(5) = A1*a;

    %% Closed-loop matrix
    Acl_i = A + B*K_nom;

    Acl{i} = Acl_i;

    %% =====================================================
    % First check: individual Hurwitz stability
    % ======================================================
    eigCL = eig(Acl_i);

    emax = max(real(eigCL));

    if emax > worstEig
        worstEig   = emax;
        worstIndex = i;
        worstA1    = A1;
        worstQ     = q;
    end

    if emax >= 0
        unstableVertices(end+1) = i;
    end

end

%% =========================================================
% Individual vertex report
% ==========================================================
fprintf('\n=========================================\n')
fprintf(' CONTROLLER: ROBUST ENVELOPE CHECK\n')
fprintf('=========================================\n')

fprintf('Number of vertices = %d\n',Nvert);

fprintf('Worst closed-loop eigenvalue = %.12e\n', ...
    worstEig);

fprintf('Worst vertex index = %d\n', ...
    worstIndex);

fprintf('A1 at worst vertex = %.6f\n', ...
    worstA1);

fprintf('Number of unstable vertices = %d\n', ...
    length(unstableVertices));

fprintf('\nq at worst vertex:\n');

fprintf(['q1 = %.12f\n' ...
         'q2 = %.12f\n' ...
         'q3 = %.12f\n' ...
         'q4 = %.12f\n' ...
         'q5 = %.12f\n' ...
         'q6 = %.12f\n'], ...
    worstQ(1),worstQ(2),worstQ(3), ...
    worstQ(4),worstQ(5),worstQ(6));

if ~isempty(unstableVertices)

    fprintf('\nRESULT:\n')
    fprintf('Some robust-envelope vertices are unstable.\n')

    fprintf('Unstable vertices:\n')
    disp(unstableVertices)

    fprintf('Common-P test is therefore skipped.\n')

    return
end

%% =========================================================
% Common quadratic Lyapunov test
% ==========================================================
P = sdpvar(nx,nx,'symmetric');

epsP   = 1e-6;
epsLMI = 1e-6;

Constraints = [ ...
    P >= epsP*eye(nx)];

for i = 1:Nvert

    M = ...
        Acl{i}'*P + ...
        P*Acl{i};

    Constraints = [ ...
        Constraints, ...
        M <= -epsLMI*eye(nx)];

end

%% =========================================================
% Solve
% ==========================================================
ops = sdpsettings( ...
    'solver','sdpt3', ...
    'verbose',1);

diagnostics = optimize( ...
    Constraints,[],ops);

fprintf('\n=========================================\n')
fprintf(' CONTROLLER COMMON-P: ROBUST ENVELOPE\n')
fprintf('=========================================\n')

disp(diagnostics.info)

%% =========================================================
% Numerical verification
% ==========================================================
if diagnostics.problem == 0

    Pval = value(P);

    minEigP = ...
        min(real(eig(Pval)));

    worstLyap = -Inf;
    worstLyapIndex = 0;

    for i = 1:Nvert

        M = ...
            Acl{i}'*Pval + ...
            Pval*Acl{i};

        emax = ...
            max(real(eig(M)));

        if emax > worstLyap

            worstLyap = emax;
            worstLyapIndex = i;

        end

    end

    fprintf('min eig(P) = %.12e\n', ...
        minEigP);

    fprintf('Worst Lyapunov eigenvalue = %.12e\n', ...
        worstLyap);

    fprintf('Worst Lyapunov vertex = %d\n', ...
        worstLyapIndex);

    if minEigP > 0 && worstLyap < 0

        fprintf('\nRESULT:\n')
        fprintf(['Common quadratic controller stability is ' ...
                 'verified over the ROBUST nonlinear ' ...
                 'q-A1 envelope.\n']);

    else

        fprintf('\nRESULT:\n')
        fprintf('Numerical certificate verification failed.\n');

    end

else

    fprintf('\nRESULT:\n')
    fprintf(['No common quadratic controller certificate ' ...
             'was found over the robust envelope.\n']);

end