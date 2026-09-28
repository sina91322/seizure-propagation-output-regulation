function [pi_vec,c] = symbolic_regulator_map_A1(A1)

% =========================================================
% SEMI-ANALYTICAL / SYMBOLIC REGULATOR MAP
% Derived from the nonlinear regulator equations
% without linearizing the Jansen-Rit sigmoid.
% =========================================================

a  = 100;
b  = 50;

A2 = 3.25;
B1 = 22;
B2 = 22;

C1 = 135;
C2 = 108;
C3 = 33.75;
C4 = 33.75;

Ad = 3;
ad = 30;

k12 = 4000; %#ok<NASGU>
k21 = 1;

p1 = 90;

% Regulation target
pi13 = 1e-3;

% Constant quantities obtained after symbolic reduction
% pi7 is the numerical root of the remaining scalar
% transcendental regulator equation.
pi7  = 0.01093013440;
pi8  = 4.345459082;
pi9  = 3.041036009;
pi14 = 0.03363118277;

% =========================================================
% Explicit A1-dependent regulator solution
% =========================================================

pi1 = A1*ad*pi13/(Ad*a);

pi3 = (B1*C4/b)*SJR_symbolic(C3*pi1);

s13 = ad*pi13/Ad;
delta1 = SJR_inv_symbolic(s13);

pi2 = pi3 + delta1;

% Feedforward regulator input
c = a*pi2/A1 ...
    - p1 ...
    - C2*SJR_symbolic(C1*pi1) ...
    - k21*pi14;

% =========================================================
% Complete regulator vector
% =========================================================

pi_vec = zeros(16,1);

pi_vec(1)  = pi1;
pi_vec(2)  = pi2;
pi_vec(3)  = pi3;

pi_vec(7)  = pi7;
pi_vec(8)  = pi8;
pi_vec(9)  = pi9;

pi_vec(13) = pi13;
pi_vec(14) = pi14;

end


% =========================================================
% Original nonlinear Jansen-Rit sigmoid
% =========================================================
function r = SJR_symbolic(v)

r = 5./(1 + exp(0.56*(6-v)));

end


% =========================================================
% Exact inverse sigmoid
% =========================================================
function v = SJR_inv_symbolic(s)

if any(s <= 0 | s >= 5)
    error('Inverse sigmoid requires 0 < s < 5.');
end

v = 6 - log(5./s - 1)/0.56;

end