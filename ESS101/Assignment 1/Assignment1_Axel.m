%% Assignment 1 - Part 1a)
clear; clc;

syms m1 m2 L g real
syms p1_x p1_y p1_z theta phi real
syms p1_x_dot p1_y_dot p1_z_dot theta_dot phi_dot real
syms p1_x_dotdot p1_y_dotdot p1_z_dotdot theta_dotdot phi_dotdot real
syms ux uy uz real

u = [ux; uy; uz];

p1 = [p1_x; p1_y; p1_z];
p1_dot = [p1_x_dot; p1_y_dot; p1_z_dot];
p1_dotdot = [p1_x_dotdot; p1_y_dotdot; p1_z_dotdot]; 
p2 = p1 + L * [sin(phi)*cos(theta); sin(phi)*sin(theta); cos(phi)];

q = [p1; theta; phi];
q_dot = [p1_dot; theta_dot; phi_dot];
q_dotdot = [p1_dotdot; theta_dotdot; phi_dotdot];

p2_dot = jacobian(p2, q) * q_dot;

T = 0.5*m1*(p1_dot.'*p1_dot) + 0.5*m2*(p2_dot.'*p2_dot);
V1 = g*m1*[0 0 1]*p1;
V2 = g*m2*[0 0 1]*p2;
V = V1 + V2;

Lag = T-V;

Q = [u; 0; 0];

v = q_dot;

dLdq = jacobian(Lag, q).';
dLdq_dot = jacobian(Lag, q_dot).';
ddt_dLdq_dot = jacobian(dLdq_dot, [q; q_dot]) * [q_dot; q_dotdot];

E = simplify(ddt_dLdq_dot - dLdq - Q); % = 0

M_1a = simplify(jacobian(E, q_dotdot))
b_1a = simplify(-(E-M_1a*q_dotdot))

%% Part 1b)

syms Z real
syms p2_x p2_y p2_z real
syms p2_x_dot p2_y_dot p2_z_dot real
syms p2_x_dotdot p2_y_dotdot p2_z_dotdot real

p2 = [p2_x; p2_y; p2_z];
p2_dot = [p2_x_dot; p2_y_dot; p2_z_dot];
p2_dotdot = [p2_x_dotdot; p2_y_dotdot; p2_z_dotdot];

q = [p1; p2];
q_dot = [p1_dot; p2_dot];
q_dotdot = [p1_dotdot; p2_dotdot];

e = p1 - p2;

C = 0.5 * (e.'*e - L^2);

V = g * (m1*p1_z + m2*p2_z);
T = 0.5*m1*(p1_dot.'*p1_dot) + 0.5*m2*(p2_dot.'*p2_dot);

Lag = T - V - Z.' * C;

v = q_dot;

gradLdq = jacobian(Lag, q).';
gradLdq_dot = jacobian(Lag, q_dot).';
ddt_gradLdq_dot = jacobian(gradLdq_dot, [q; q_dot]) * [q_dot; q_dotdot];

Q = [u; 0; 0; 0];

dCdq = jacobian(C, q);
ddq_dCdq = jacobian(dCdq*q_dot, q);

E = simplify(ddt_gradLdq_dot - gradLdq - Q);

M_1b = simplify(jacobian(E, q_dotdot))
b_1b = simplify(-(E-M_1b*q_dotdot))

c_dotdot = dCdq * q_dotdot + ddq_dCdq * q_dot;

%% Part 2a)

a = jacobian(C, q).';

LHS_mat = [M_1b a; a.' 0];

zero_qdd = sym(zeros(size(q_dotdot)));

old = [q_dotdot; Z];
new = [zero_qdd;  sym(0)];

E_RHS        = subs(E,        old, new);
c_dotdot_RHS = subs(c_dotdot, old, new);

x = [q_dotdot; Z];
LHS = LHS_mat*x;
RHS = -[E_RHS; c_dotdot_RHS];

a = a
c = RHS

%% Part 2b)
x_sol = LHS_mat \ c





