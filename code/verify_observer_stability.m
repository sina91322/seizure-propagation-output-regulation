clc
clear
close all

%% =========================================================
% FINAL OBSERVER COMMON-P TEST
%
% Uncertainties:
%   q1,...,q6
%   A1
%   E4 = df4/dA1
%   E5 = df5/dA1
%
% Total vertices:
%   2^9 = 512
%
% Final semi-analytical regulator map:
%   symbolic_regulator_map_A1
% ==========================================================

%% =========================================================
% Robust nonlinear q-envelope
% Valid for regulated operating region t >= 5.7 s
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

%% =========================================================
% A1 uncertainty
% ==========================================================
A1_min = 6.5;
A1_max = 7.5;

%% =========================================================
% Load FINAL observer gains
% ==========================================================
load('observer_A1_gains.mat', ...
    'Lx_A1','Lw_A1');

Lw_test = 0.05*Lw_A1;

L_aug = [ ...
    Lx_A1
    Lw_test ];

fprintf('\n=========================================\n')
fprintf(' FINAL OBSERVER COMMON-P TEST\n')
fprintf('=========================================\n')

fprintf('Original Lw = %.12f\n',Lw_A1);
fprintf('Used Lw     = %.12f\n',Lw_test);

%% =========================================================
% Dimensions
% ==========================================================
nx = 16;
na = 17;

%% =========================================================
% Fixed model parameters
% ==========================================================
a = 100;
b = 50;

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

p1 = 90;

%% =========================================================
% Output matrix
%
% y = x2 - x3
% ==========================================================
C = zeros(1,nx);

C(2) = 1;
C(3) = -1;

C_aug = [C 0];

%% =========================================================
% Jansen-Rit sigmoid
% ==========================================================
Sfun = @(v) ...
    5./(1 + exp(0.56*(6-v)));

%% =========================================================
% Determine numerical envelope of E_A1(A1)
%
% E_A1 = df/dA1 evaluated on regulator manifold
% ==========================================================
NAgrid = 10001;

A1_grid = linspace( ...
    A1_min, ...
    A1_max, ...
    NAgrid);

E4_grid = zeros(NAgrid,1);
E5_grid = zeros(NAgrid,1);

for kk = 1:NAgrid

    A1g = A1_grid(kk);

    [pi_A1,c_A1] = ...
        symbolic_regulator_map_A1(A1g);

    E4_grid(kk) = ...
        a*Sfun( ...
        pi_A1(2)-pi_A1(3));

    E5_grid(kk) = ...
        a*( ...
        p1 ...
        + c_A1 ...
        + C2*Sfun(C1*pi_A1(1)) ...
        + k21*pi_A1(14));

end

%% =========================================================
% Conservative numerical bounds
% ==========================================================
E4_min = min(E4_grid);
E4_max = max(E4_grid);

E5_min = min(E5_grid);
E5_max = max(E5_grid);

%% Add a small numerical safety margin
margin = 1e-6;

E4_span = max(1,abs(E4_max-E4_min));
E5_span = max(1,abs(E5_max-E5_min));

E4_min = E4_min - margin*E4_span;
E4_max = E4_max + margin*E4_span;

E5_min = E5_min - margin*E5_span;
E5_max = E5_max + margin*E5_span;

fprintf('\n=========================================\n')
fprintf(' E_A1 ENVELOPE\n')
fprintf('=========================================\n')

fprintf('E4 min = %.12e\n',E4_min);
fprintf('E4 max = %.12e\n',E4_max);

fprintf('E5 min = %.12e\n',E5_min);
fprintf('E5 max = %.12e\n',E5_max);

%% =========================================================
% Build 512 vertices
%
% bits:
%   1...6 -> q1...q6
%   7     -> A1
%   8     -> E4
%   9     -> E5
% ==========================================================
Nvert = 512;

Aaug = cell(Nvert,1);

unstableVertices = [];

worstEig = -Inf;
worstIndex = 0;

worstA1 = NaN;
worstQ = NaN(1,6);
worstE4 = NaN;
worstE5 = NaN;

for i = 1:Nvert

    bits = dec2bin(i-1,9)-'0';

    %% =====================================================
    % q vertex
    % ======================================================
    qbits = bits(1:6);

    q = ...
        qmin + ...
        qbits.*(qmax-qmin);

    q1 = q(1);
    q2 = q(2);
    q3 = q(3);
    q4 = q(4);
    q5 = q(5);
    q6 = q(6);

    %% =====================================================
    % A1 vertex
    % ======================================================
    if bits(7) == 0
        A1 = A1_min;
    else
        A1 = A1_max;
    end

    %% =====================================================
    % E4 vertex
    % ======================================================
    if bits(8) == 0
        E4 = E4_min;
    else
        E4 = E4_max;
    end

    %% =====================================================
    % E5 vertex
    % ======================================================
    if bits(9) == 0
        E5 = E5_min;
    else
        E5 = E5_max;
    end

    %% =====================================================
    % State Jacobian A(q,A1)
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

    %% =====================================================
    % Column 1
    % ======================================================
    A(4,1) = -a^2;

    A(4,2) = ...
        A1*a*q1;

    A(4,3) = ...
        -A1*a*q1;

    A(4,4) = ...
        -2*a;

    A(5,1) = ...
        A1*a*C2*C1*q2;

    A(5,2) = ...
        -a^2;

    A(5,5) = ...
        -2*a;

    A(5,14) = ...
        A1*a*k21;

    A(6,1) = ...
        B1*b*C4*C3*q3;

    A(6,3) = ...
        -b^2;

    A(6,6) = ...
        -2*b;

    %% =====================================================
    % Column 2
    % ======================================================
    A(10,7) = ...
        -a^2;

    A(10,8) = ...
        A2*a*q4;

    A(10,9) = ...
        -A2*a*q4;

    A(10,10) = ...
        -2*a;

    A(11,7) = ...
        A2*a*C2*C1*q5;

    A(11,8) = ...
        -a^2;

    A(11,11) = ...
        -2*a;

    A(11,13) = ...
        A2*a*k12;

    A(12,7) = ...
        B2*b*C4*C3*q6;

    A(12,9) = ...
        -b^2;

    A(12,12) = ...
        -2*b;

    %% =====================================================
    % Inter-column pathways
    % ======================================================
    A(15,2) = ...
        Ad*ad*q1;

    A(15,3) = ...
        -Ad*ad*q1;

    A(15,13) = ...
        -ad^2;

    A(15,15) = ...
        -2*ad;

    A(16,8) = ...
        Ad*ad*q4;

    A(16,9) = ...
        -Ad*ad*q4;

    A(16,14) = ...
        -ad^2;

    A(16,16) = ...
        -2*ad;

    %% =====================================================
    % Uncertain parameter column
    % ======================================================
    E_A1 = zeros(nx,1);

    E_A1(4) = E4;
    E_A1(5) = E5;

    %% =====================================================
    % Augmented 17x17 model
    % ======================================================
    Aaug_i = [ ...
        A, E_A1;
        zeros(1,nx), 0 ];

    Aaug{i} = Aaug_i;

    %% =====================================================
    % Fixed-gain observer error matrix
    % ======================================================
    Aerr = ...
        Aaug_i + ...
        L_aug*C_aug;

    emax = ...
        max(real(eig(Aerr)));

    if emax > worstEig

        worstEig = emax;
        worstIndex = i;

        worstA1 = A1;
        worstQ = q;

        worstE4 = E4;
        worstE5 = E5;

    end

    if emax >= 0

        unstableVertices(end+1) = i;

    end

end

%% =========================================================
% Fixed-gain report
% ==========================================================
fprintf('\n=========================================\n')
fprintf(' OBSERVER CHECK ON 512 VERTICES\n')
fprintf('=========================================\n')

fprintf('Worst observer eigenvalue = %.12e\n', ...
    worstEig);

fprintf('Worst vertex index = %d\n', ...
    worstIndex);

fprintf('A1 at worst vertex = %.6f\n', ...
    worstA1);

fprintf('Number unstable observer vertices = %d\n', ...
    length(unstableVertices));

fprintf('\nq at worst vertex:\n');

for ii = 1:6

    fprintf('q%d = %.12f\n', ...
        ii,worstQ(ii));

end

fprintf('\nE values at worst vertex:\n');

fprintf('E4 = %.12e\n',worstE4);
fprintf('E5 = %.12e\n',worstE5);

%% =========================================================
% Stop if any vertex is unstable
% ==========================================================
if ~isempty(unstableVertices)

    fprintf('\nRESULT:\n')
    fprintf(['Some observer vertices are not Hurwitz.\n' ...
             'Common-P test is skipped.\n']);

    return

end

%% =========================================================
% Common quadratic Lyapunov matrix
% ==========================================================
Pobs = ...
    sdpvar(na,na,'symmetric');

epsP   = 1e-6;
epsLMI = 1e-6;

Constraints = [ ...
    Pobs >= ...
    epsP*eye(na)];

for i = 1:Nvert

    Aerr = ...
        Aaug{i} + ...
        L_aug*C_aug;

    M = ...
        Aerr'*Pobs + ...
        Pobs*Aerr;

    Constraints = [ ...
        Constraints, ...
        M <= -epsLMI*eye(na)];

end

%% =========================================================
% Solve LMI
% ==========================================================
ops = ...
    sdpsettings( ...
    'solver','sdpt3', ...
    'verbose',1);

diagnostics = ...
    optimize( ...
    Constraints, ...
    [], ...
    ops);

fprintf('\n=========================================\n')
fprintf(' OBSERVER COMMON-P: FINAL 512-VERTEX TEST\n')
fprintf('=========================================\n')

disp(diagnostics.info)

%% =========================================================
% Verification
% ==========================================================
if diagnostics.problem == 0

    Pobs_val = ...
        value(Pobs);

    minEigPobs = ...
        min(real(eig(Pobs_val)));

    worstLyap = -Inf;
    worstLyapIndex = 0;

    for i = 1:Nvert

        Aerr = ...
            Aaug{i} + ...
            L_aug*C_aug;

        M = ...
            Aerr'*Pobs_val + ...
            Pobs_val*Aerr;

        emax = ...
            max(real(eig(M)));

        if emax > worstLyap

            worstLyap = emax;
            worstLyapIndex = i;

        end

    end

    fprintf('\nmin eig(Pobs) = %.12e\n', ...
        minEigPobs);

    fprintf( ...
        'Worst observer Lyapunov eigenvalue = %.12e\n', ...
        worstLyap);

    fprintf( ...
        'Worst Lyapunov vertex = %d\n', ...
        worstLyapIndex);

    if ...
        minEigPobs > 0 && ...
        worstLyap < 0

        fprintf('\nRESULT:\n');

        fprintf([ ...
            'Common quadratic observer stability is ' ...
            'verified over the final bounded ' ...
            'q-A1-E region.\n']);

    else

        fprintf('\nRESULT:\n');

        fprintf([ ...
            'Numerical certificate verification ' ...
            'failed.\n']);

    end

else

    fprintf('\nRESULT:\n');

    fprintf([ ...
        'No common quadratic observer certificate ' ...
        'was found over the final bounded region.\n']);

end